"""Estimate and test a plate-position correction for the validation plates.

The screen estimates a plate-position effect (a "tilt": an additive row plus column effect
for every feature) from compounds that were placed at random on 11 platemaps
(``position_correction_utils``).
The validation plates have one layout, shared by five replicate plates: both kinds of DMSO
(healthy and failing cells) in blocks of columns over all rows, and 15 compounds in
columns 5 to 8, three on each plate at positions that differ between plates.
Compounds do not average out of a tilt estimate with so few compounds, but every group of
wells (a compound, or a kind of DMSO) sits at several positions, so we fit one effect per
group together with the tilt (``fit_compound_adjusted_tilt``).
The strength of the tilt differs between plates, so every plate has its own amplitude, like
the platemaps of the screen (``shrunk_plate_amplitudes``).

The validation plates keep their own feature names (the d2 channel is ``Golgi``), and the
tilt is fit on the validation plates only.

This module holds what is specific to the validation plates.
It reuses ``position_correction_utils`` and ``position_correction_effect_utils`` wherever
it can, so the validation plates are processed like the screen: the fit has the same
format (``coef``, ``center``, ``features`` and ``amplitudes``), so
``position_correction_utils.save_fit``, ``load_fit`` and ``correct_cells`` work with it.
"""

import json
import pathlib
from collections.abc import Callable, Iterable, Sequence

import numpy as np
import pandas as pd
import position_correction_effect_utils as pce
import position_correction_utils as pcu
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HEALTHY_CELL_TYPE = "nonfailing"  # the screen calls this "healthy"
FAILING_CELL_TYPE = "failing"
N_ROWS = len(pcu.ROWS)
# the channels of the validation plates (the d2 channel is Golgi, and PM in the screen)
CHANNELS = ("Actin", "DNA", "ER", "Mito", "Golgi")


def load_platemap_plate_numbers(platemap_dir: str | pathlib.Path) -> dict[str, int]:
    """Read the plate number of every validation plate from its plate map.

    Parameters
    ----------
    platemap_dir : str | pathlib.Path
        Folder with one ``<plate barcode>_platemap.csv`` per plate, with the columns
        ``plate_barcode`` and ``plate_number``.

    Returns
    -------
    dict[str, int]
        Plate barcode mapped to its plate number.

    Raises
    ------
    FileNotFoundError
        If the folder has no plate map.
    """
    numbers = {}
    for path in sorted(pathlib.Path(platemap_dir).glob("*_platemap.csv")):
        platemap = pd.read_csv(path, usecols=["plate_barcode", "plate_number"])
        numbers[str(platemap["plate_barcode"].iloc[0])] = int(platemap["plate_number"].iloc[0])
    if not numbers:
        raise FileNotFoundError(f"No *_platemap.csv files in {platemap_dir}")
    return numbers


def discover_validation_profiles(
    profile_dir: str | pathlib.Path,
    suffix: str = "_sc_normalized.parquet",
    platemap_dir: str | pathlib.Path | None = None,
) -> pd.DataFrame:
    """Find the single-cell profile file of every validation plate.

    Looks for ``<plate><suffix>`` directly in ``profile_dir``.
    The result has the layout that ``position_correction_utils.build_well_medians``
    expects, with the plate number as the platemap number, the key of the amplitude of
    the plate in the fit.
    The plate number is the one of the plate map of the plate (``platemap_dir``), so it
    does not depend on which other plates are in the folder.
    Without ``platemap_dir`` it is the position of the plate (1 to n) in the order of the
    plate names, which is only right when the folder has all five plates.

    Parameters
    ----------
    profile_dir : str | pathlib.Path
        Folder with the single-cell profile files of the validation plates.
    suffix : str, optional
        End of the file names, by default ``"_sc_normalized.parquet"``.
    platemap_dir : str | pathlib.Path | None, optional
        Folder with the plate maps (see ``load_platemap_plate_numbers``).

    Returns
    -------
    pd.DataFrame
        One row per plate, sorted by plate number, with the columns ``plate`` (str),
        ``platemap`` (int) and ``path`` (pathlib.Path).
        The table is empty, with no columns, when no file matches.

    Raises
    ------
    ValueError
        If ``platemap_dir`` is given and a plate has no plate map.
    """
    paths = sorted(pathlib.Path(profile_dir).glob(f"*{suffix}"))
    if platemap_dir is None:
        records = [
            {"plate": path.name.removesuffix(suffix), "platemap": number, "path": path}
            for number, path in enumerate(paths, start=1)
        ]
        return pd.DataFrame(records)
    numbers = load_platemap_plate_numbers(platemap_dir)
    plates = [path.name.removesuffix(suffix) for path in paths]
    unknown = [plate for plate in plates if plate not in numbers]
    if unknown:
        raise ValueError(f"No plate map for the plates {unknown} in {platemap_dir}")
    records = [
        {"plate": plate, "platemap": numbers[plate], "path": path}
        for plate, path in zip(plates, paths, strict=True)
    ]
    return pd.DataFrame(records).sort_values("platemap", ignore_index=True)


def dropped_features_path(cache_path: str | pathlib.Path) -> pathlib.Path:
    """Path of the file that records the features that building the well cache dropped."""
    cache_path = pathlib.Path(cache_path)
    return cache_path.with_name(cache_path.stem + "_dropped_features.json")


def check_feature_selection_is_plate_independent(cache_path: str | pathlib.Path) -> None:
    """Check that the held-out plate does not decide which features the tests use.

    Building the well cache drops every feature with a missing value in any well of any
    plate. A test that holds out a plate keeps a feature that is missing on that plate
    only, so a feature that is missing on exactly one plate would make the feature set of
    that fold differ from the set of the cache. The check reads the record that the build
    wrote (``dropped_features_path``).

    Parameters
    ----------
    cache_path : str | pathlib.Path
        Parquet file of the well medians (see ``load_or_build_well_medians``).

    Raises
    ------
    FileNotFoundError
        If the cache has no record of its dropped features (it was built by an older
        version): delete the cache and build it again.
    ValueError
        If a feature is missing on exactly one plate. The message names the features.
    """
    record = dropped_features_path(cache_path)
    if not record.is_file():
        raise FileNotFoundError(
            f"{record} is missing, so we cannot tell which features the cache dropped. "
            f"Delete {cache_path} and build it again."
        )
    dropped = json.loads(record.read_text())
    one_plate = {feature: plates for feature, plates in dropped.items() if len(plates) == 1}
    if one_plate:
        raise ValueError(
            f"{len(one_plate)} features are missing on exactly one plate, so the feature set "
            f"of a test that holds that plate out would differ from the cache: {one_plate}"
        )


def load_or_build_well_medians(
    cache_path: str | pathlib.Path,
    profile_dir: str | pathlib.Path,
    n_jobs: int = 1,
    platemap_dir: str | pathlib.Path | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Load the median profile of every well, or build it from the single-cell profiles.

    Reading the five plates takes a few minutes, so the table is cached in a parquet file.
    When the cache is missing, the function summarizes every well of every plate as the
    median profile of its cells (``position_correction_utils.build_well_medians``, which
    drops the features with a missing value in any well) and saves the cache.
    Delete the cache when the single-cell profiles change.

    Parameters
    ----------
    cache_path : str | pathlib.Path
        Parquet file of the well medians.
    profile_dir : str | pathlib.Path
        Folder with the ``<plate>_sc_normalized.parquet`` files, only read when the cache
        is missing.
    n_jobs : int, optional
        Number of plates to read in parallel, by default 1.
    platemap_dir : str | pathlib.Path | None, optional
        Folder with the plate maps, which give the plate numbers (see
        ``discover_validation_profiles``). Only used when the cache is missing.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        The well table, with one row per plate and well and the columns ``plate``,
        ``platemap``, ``well``, ``treatment``, ``cell_type``, ``n`` (the number of cells in
        the well) and the features, and the list of feature names.

    Raises
    ------
    FileNotFoundError
        If the cache is missing and ``profile_dir`` has no profile files.
    """
    cache_path = pathlib.Path(cache_path)
    meta_columns = ["plate", "platemap", "well", "treatment", "cell_type", "n"]
    if cache_path.exists():
        well_medians = pd.read_parquet(cache_path)
        return well_medians, [c for c in well_medians.columns if c not in meta_columns]
    plate_table = discover_validation_profiles(
        pathlib.Path(profile_dir).resolve(strict=True), platemap_dir=platemap_dir
    )
    if plate_table.empty:
        raise FileNotFoundError(f"No *_sc_normalized.parquet files in {profile_dir}")
    well_medians, features, dropped = pcu.build_well_medians(
        plate_table, n_jobs=n_jobs, return_dropped=True
    )
    well_medians.to_parquet(cache_path, index=False)
    # record what the build dropped, so that the held-out tests can check the features
    dropped_features_path(cache_path).write_text(json.dumps(dropped, indent=1))
    return well_medians, features


def fit_compound_adjusted_tilt(
    well_medians: pd.DataFrame,
    features: Sequence[str],
    min_cells: int = 50,
    include_controls: bool = False,
    ridge: float = 1e-6,
) -> dict:
    """Fit a tilt map that accounts for the effect of each group of wells.

    The screen places every compound at one position, so it can only average compound
    effects out over many compounds.
    Here every group of wells (a compound, or with ``include_controls`` the healthy DMSO
    and the failing DMSO) sits at several positions, so we fit one effect for each group
    together with an additive row plus column effect per feature::

        profile(group, row, column) = group_effect(group)
                                      + row_effect(row) + column_effect(column)

    The group effect absorbs everything that is the same for all wells of a group, and the
    row and column effects are what differs between the positions of a group.
    The column effects have no reference column.
    When the groups sit in blocks of columns (for example the failing DMSO only in columns
    2, 4, and 10, the healthy DMSO in columns 3, 9, and 11, and the compounds in columns 5
    to 8), a group effect can absorb the offset of its block of columns, so those offsets
    cannot be estimated and an unconstrained fit splits them between group and column
    effects arbitrarily, which distorts the difference between the groups.
    A tiny ridge penalty on the column effects (``ridge``) sets what the data cannot tell
    to zero, so the group effects keep the differences between groups; the effects that the
    data do determine are unchanged.
    The fit has the same two passes as ``position_correction_utils.fit_tilt_map``
    (least squares, then a refit after winsorizing the residuals at three robust standard
    deviations), and it averages the replicate wells of a group (the same position on the
    replicate plates) first.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians of the validation plates (``load_or_build_well_medians``).
    features : Sequence[str]
        Features to model.
    min_cells : int, optional
        Minimum number of cells for a well to be used, by default 50.
    include_controls : bool, optional
        Whether the DMSO wells enter the fit, as two groups (healthy and failing cells).
        By default only the compound wells do, one group for each compound. The corrected
        layout of the validation plates has the controls of both cell types at many rows
        and columns, so the fit needs them to see rows B and G.
    ridge : float, optional
        Ridge penalty on the column effects, by default 1e-6 (small enough not to change
        the effects the data determine).

    Returns
    -------
    dict
        The keys ``coef`` and ``center`` (a tilt map in the layout of
        ``position_correction_utils``, so ``position_correction_utils.tilt`` works with it,
        with the tilt averaging to zero over the wells that the tilt map was fit on) and
        ``features``.
    """
    features = list(features)
    wells = well_medians[well_medians["n"] >= min_cells]
    if include_controls:
        group = wells["treatment"] + "|" + wells["cell_type"]
    else:
        wells = wells[~wells["treatment"].isin(pcu.CONTROL_TREATMENTS)]
        group = wells["treatment"]
    averaged = (
        wells.assign(group=group)
        .groupby(["group", "well"])[features]
        .mean()
        .reset_index()
    )
    rows, cols = pcu.well_index(averaged["well"])
    groups = sorted(averaged["group"].unique())

    group_terms = (averaged["group"].to_numpy()[:, None] == np.array(groups)).astype(float)
    # row B (index 0) is the reference row (the group terms are the intercept); the columns
    # have one term each, with a ridge penalty instead of a reference column
    row_terms = np.eye(pcu.N_ROW_TERMS)[rows][:, 1:]
    n_columns = pcu.N_COL_TERMS + 1
    column_terms = np.eye(n_columns)[cols]
    design = np.hstack([group_terms, row_terms, column_terms])
    penalty = np.hstack(
        [
            np.zeros((n_columns, design.shape[1] - n_columns)),
            np.sqrt(ridge) * np.eye(n_columns),
        ]
    )
    profiles = averaged[features].to_numpy(float)

    def solve(values: np.ndarray) -> np.ndarray:
        stacked = np.vstack([values, np.zeros((n_columns, values.shape[1]))])
        return np.linalg.lstsq(np.vstack([design, penalty]), stacked, rcond=None)[0]

    beta = solve(profiles)
    mad_to_sigma = 1.4826
    resid = profiles - design @ beta
    resid_center = np.median(resid, axis=0)
    mad = mad_to_sigma * np.median(np.abs(resid - resid_center), axis=0) + 1e-9
    winsorized = design @ beta + np.clip(
        resid, resid_center - 3 * mad, resid_center + 3 * mad
    )
    beta = solve(winsorized)

    row_effects = beta[len(groups) : len(groups) + pcu.N_ROW_TERMS - 1]  # rows C to G
    column_effects = beta[len(groups) + pcu.N_ROW_TERMS - 1 :]  # columns 2 to 11
    coef = np.vstack(
        [
            np.zeros((1, len(features))),  # row B, the reference
            row_effects,
            column_effects[1:] - column_effects[0],  # columns 3 to 11 relative to column 2
        ]
    )
    center = (pcu.design_matrix(rows, cols) @ coef).mean(axis=0)
    return {"coef": coef, "center": center, "features": features}


def with_plate_amplitudes(fit: dict, amplitudes: pd.DataFrame) -> dict:
    """Attach the amplitude of every plate to a tilt map.

    The amplitudes are keyed by the plate number, which is the platemap number of the
    plate in ``discover_validation_profiles``, so the result works with
    ``position_correction_utils.save_fit``, ``position_correction_utils.correct_cells`` and
    ``position_correction_effect_utils.correct_well_medians``.

    Parameters
    ----------
    fit : dict
        Tilt map with the keys ``coef``, ``center`` and ``features``.
    amplitudes : pd.DataFrame
        Output of ``shrunk_plate_amplitudes``, with the columns ``plate_number`` and
        ``shrunk_amplitude``.

    Returns
    -------
    dict
        Copy of ``fit`` with the key ``amplitudes`` (plate number mapped to its amplitude)
        and the key ``plate_numbers`` (plate name mapped to its plate number), so that a
        plate gets its own amplitude whatever other plates are processed with it.
    """
    return {
        **fit,
        "amplitudes": {
            int(number): float(value)
            for number, value in zip(
                amplitudes["plate_number"], amplitudes["shrunk_amplitude"], strict=True
            )
        },
        "plate_numbers": {
            str(plate): int(number)
            for plate, number in zip(
                amplitudes["plate"], amplitudes["plate_number"], strict=True
            )
        },
    }


def validate_plate_amplitudes(
    amplitudes: pd.DataFrame, expected_plates: Iterable[str]
) -> None:
    """Check that every expected plate has a finite amplitude before the fit is saved.

    Parameters
    ----------
    amplitudes : pd.DataFrame
        Output of ``shrunk_plate_amplitudes``.
    expected_plates : Iterable[str]
        Barcodes of the plates that the fit has to cover.

    Raises
    ------
    ValueError
        If a plate is missing (for example because no well of it has enough cells) or
        has an amplitude, a standard error, or a shrunk amplitude that is not finite.
        The message names the plates.
    """
    present = set(amplitudes["plate"])
    missing = [plate for plate in expected_plates if plate not in present]
    finite = np.isfinite(
        amplitudes[["amplitude", "standard_error", "shrunk_amplitude"]].to_numpy(float)
    ).all(axis=1)
    invalid = amplitudes.loc[~finite, "plate"].tolist()
    if missing or invalid:
        raise ValueError(
            f"The fit does not cover every plate. Missing plates: {missing}. "
            f"Plates without a finite amplitude: {invalid}."
        )


def save_validation_fit(fit: dict, path: str | pathlib.Path) -> None:
    """Save a fit of the validation plates with the plate that owns each amplitude.

    This is ``position_correction_utils.save_fit`` with the arrays ``plate_names`` and
    ``plate_numbers`` (``fit["plate_numbers"]``) added. ``position_correction_utils.load_fit``
    still reads the file and ignores them.
    """
    amplitudes = fit["amplitudes"]
    plate_numbers = fit["plate_numbers"]
    np.savez_compressed(
        path,
        coef=fit["coef"].astype("float32"),
        center=fit["center"].astype("float32"),
        features=np.array(fit["features"]),
        amplitude_platemaps=np.array(list(amplitudes), dtype=int),
        amplitude_values=np.array(list(amplitudes.values()), dtype=float),
        plate_names=np.array(list(plate_numbers)),
        plate_numbers=np.array(list(plate_numbers.values()), dtype=int),
    )


def load_validation_fit(path: str | pathlib.Path) -> dict:
    """Load a fit written by ``save_validation_fit``.

    Returns the output of ``position_correction_utils.load_fit`` with the key
    ``plate_numbers`` (plate name mapped to its plate number).

    Raises
    ------
    ValueError
        If the file has no plate names (it was saved by ``save_fit``): rerun
        ``1.fit_position_correction.ipynb``.
    """
    fit = pcu.load_fit(path)
    stored = np.load(path, allow_pickle=False)
    if "plate_names" not in stored.files:
        raise ValueError(
            f"{path} does not say which plate owns each amplitude. "
            "Run 4a.position_correction/1.fit_position_correction.ipynb again."
        )
    fit["plate_numbers"] = dict(
        zip(stored["plate_names"].tolist(), stored["plate_numbers"].tolist(), strict=True)
    )
    return fit


def plate_number_in_fit(fit: dict, plate: str) -> int:
    """Return the plate number that the amplitude of ``plate`` is stored under.

    Raises
    ------
    KeyError
        If the fit has no amplitude for the plate, or the amplitude is not finite.
    """
    if plate not in fit["plate_numbers"]:
        raise KeyError(
            f"The fit has no amplitude for the plate {plate}; it covers {sorted(fit['plate_numbers'])}"
        )
    number = fit["plate_numbers"][plate]
    if not np.isfinite(fit["amplitudes"].get(number, np.nan)):
        raise KeyError(f"The amplitude of the plate {plate} is not finite")
    return int(number)


def shrunk_plate_amplitudes(
    well_medians: pd.DataFrame,
    fit: dict,
    min_cells: int = 50,
    n_boot: int = 200,
    seed: int = 0,
) -> pd.DataFrame:
    """Measure the amplitude of every plate and shrink it like the screen does.

    The amplitude of a plate is how strongly the tilt shows in it (``plate_amplitudes``).
    Like the amplitude of a platemap in the screen, it is noisy, so we shrink it toward the
    average of the other plates (``position_correction_utils.shrink_amplitude``).

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians of the validation plates.
    fit : dict
        Tilt map with the keys ``coef``, ``center`` and ``features``.
    min_cells : int, optional
        Minimum number of cells for a well to be used, by default 50.
    n_boot : int, optional
        Number of resamples of the standard error, by default 200.
    seed : int, optional
        Seed of the resampling, by default 0.

    Returns
    -------
    pd.DataFrame
        The output of ``plate_amplitudes`` with the column ``shrunk_amplitude`` added.
    """
    table = plate_amplitudes(well_medians, fit, min_cells, n_boot, seed)
    # a plate whose amplitude or standard error is undefined is left out of the shrinking
    # and gets no shrunk amplitude (``validate_plate_amplitudes`` reports it)
    usable = np.isfinite(table[["amplitude", "standard_error"]].to_numpy(float)).all(axis=1)
    raw = {
        int(row.plate_number): (row.amplitude, row.standard_error)
        for row, ok in zip(table.itertuples(), usable, strict=True)
        if ok
    }
    table["shrunk_amplitude"] = [
        pcu.shrink_amplitude(raw, int(number)) if ok and len(raw) >= 3 else np.nan
        for number, ok in zip(table["plate_number"], usable, strict=True)
    ]
    return table


def _tilt_scale(profiles: np.ndarray, tilt_values: np.ndarray) -> float:
    """``position_correction_utils.amplitude``, or NaN when the tilt has no variation."""
    if not (tilt_values * tilt_values).sum() > 0:
        return float("nan")
    return pcu.amplitude(profiles, tilt_values)


def plate_amplitudes(
    well_medians: pd.DataFrame,
    fit: dict,
    min_cells: int = 50,
    n_boot: int = 200,
    seed: int = 0,
    min_valid_fraction: float = 0.9,
) -> pd.DataFrame:
    """Measure how strongly the tilt shows in each replicate plate.

    The amplitude of a plate is the scale of the tilt that best matches its wells
    (``position_correction_utils.amplitude``) after the mean of each group of wells (a
    compound, or the healthy or the failing DMSO) is removed from the profiles and from the
    tilt, so that what a group does is not counted as position.
    Its standard error comes from resampling the wells of the plate.
    The tilt was fit on all five plates, so the amplitudes are around 1, and a plate that
    differs a lot would show that the position effect is not the same on every plate.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians of the validation plates.
    fit : dict
        Tilt map with the keys ``coef``, ``center`` and ``features``.
    min_cells : int, optional
        Minimum number of cells for a well to be used, by default 50.
    n_boot : int, optional
        Number of resamples of the standard error, by default 200.
    seed : int, optional
        Seed of the resampling, by default 0.
    min_valid_fraction : float, optional
        A resample of the wells of a plate that has no variation in the tilt has no
        amplitude. The standard error is NaN when fewer than this fraction of the
        resamples have one, by default 0.9.

    Returns
    -------
    pd.DataFrame
        One row per plate with the columns ``plate``, ``plate_number`` (the ``platemap``
        column of ``well_medians``, which is the plate number of the plate map),
        ``amplitude`` and ``standard_error``. The amplitude is NaN when the tilt has no
        variation within the groups of the plate.
    """
    wells = well_medians[well_medians["n"] >= min_cells].reset_index(drop=True)
    features = list(fit["features"])
    rows, cols = pcu.well_index(wells["well"])
    profiles = wells[features].to_numpy(float)
    tilt_values = pcu.tilt(fit["coef"], fit["center"], rows, cols)
    group = (wells["treatment"] + "|" + wells["cell_type"]).to_numpy()
    rng = np.random.default_rng(seed)
    plate_numbers = wells.groupby("plate")["platemap"].first()
    records = []
    for plate in sorted(wells["plate"].unique()):
        number = int(plate_numbers[plate])
        in_plate = (wells["plate"] == plate).to_numpy()
        labels = group[in_plate]
        centered_profiles = profiles[in_plate] - pd.DataFrame(profiles[in_plate]).groupby(labels).transform("mean").to_numpy()
        centered_tilt = tilt_values[in_plate] - pd.DataFrame(tilt_values[in_plate]).groupby(labels).transform("mean").to_numpy()
        n_wells = in_plate.sum()
        boots = []
        for _ in range(n_boot):
            draw = rng.integers(0, n_wells, n_wells)
            boots.append(_tilt_scale(centered_profiles[draw], centered_tilt[draw]))
        valid = np.array(boots)[np.isfinite(boots)]
        records.append(
            {
                "plate": plate,
                "plate_number": number,
                "amplitude": _tilt_scale(centered_profiles, centered_tilt),
                "standard_error": float(np.std(valid))
                if len(valid) >= max(2, min_valid_fraction * n_boot)
                else np.nan,
            }
        )
    return pd.DataFrame(records)


def parse_feature_family(feature: str) -> dict[str, str]:
    """Split a CellProfiler feature name into compartment, measurement and channel.

    This is ``position_correction_effect_utils.parse_feature_family`` for the channels of
    the validation plates (``CHANNELS``, with ``Golgi`` where the screen has ``PM``).

    Parameters
    ----------
    feature : str
        Feature name such as ``"Cells_Texture_SumAverage_DNA_3_00_256"``.

    Returns
    -------
    dict[str, str]
        The keys ``compartment``, ``measurement`` and ``channel``; the channel is
        ``"Multiple channels"`` for features that combine channels and ``"No channel"``
        for features without a channel.
    """
    parts = feature.split("_")
    channels = []
    for part in parts[2:]:
        if part in CHANNELS and part not in channels:
            channels.append(part)
    if not channels:
        channel = pce.NO_CHANNEL
    elif len(channels) > 1:
        channel = pce.MULTI_CHANNEL
    else:
        channel = channels[0]
    return {"compartment": parts[0], "measurement": parts[1], "channel": channel}


def tilt_size_by_feature(fit: dict) -> pd.DataFrame:
    """Measure the size of the tilt of every feature and label its family.

    The size of a feature is the root mean square of its tilt over all 60 positions of the
    plate, with an amplitude of 1, in the units of the normalized features (standard
    deviations of a single cell).
    This is ``position_correction_effect_utils.tilt_size_by_feature`` with the channels of
    the validation plates.

    Parameters
    ----------
    fit : dict
        Tilt map with the keys ``coef``, ``center`` and ``features``.

    Returns
    -------
    pd.DataFrame
        One row per feature with the columns ``feature``, ``tilt_rms``, ``compartment``,
        ``measurement`` and ``channel``.
    """
    rows, cols = np.meshgrid(
        np.arange(pcu.N_ROW_TERMS), np.arange(pcu.N_COL_TERMS + 1), indexing="ij"
    )
    tilt = pcu.tilt(fit["coef"], fit["center"], rows.ravel(), cols.ravel())
    sizes = pd.DataFrame(
        {"feature": fit["features"], "tilt_rms": np.sqrt((tilt**2).mean(axis=0))}
    )
    families = pd.DataFrame([parse_feature_family(name) for name in sizes["feature"]])
    return pd.concat([sizes, families], axis=1)


def leave_one_plate_out_correction(
    wells: pd.DataFrame, features: Sequence[str], min_cells: int = 50
) -> np.ndarray:
    """Correct every plate with a tilt that was fit on the other plates.

    For each plate, the tilt is fit on the wells of the other plates (compounds and
    controls, ``fit_compound_adjusted_tilt`` with ``include_controls``) and subtracted from
    the wells of the plate, so the correction never saw the plate.

    Parameters
    ----------
    wells : pd.DataFrame
        Wells with at least ``min_cells`` cells, with the columns ``plate``, ``well``,
        ``treatment``, ``cell_type``, ``n`` and the features.
    features : Sequence[str]
        Features to correct.
    min_cells : int, optional
        Minimum cells of a well, by default 50.

    Returns
    -------
    np.ndarray
        Corrected profiles with shape ``(len(wells), len(features))``, in the order of
        ``wells``.
    """
    features = list(features)
    rows, cols = pcu.well_index(wells["well"])
    raw = wells[features].to_numpy(float)
    corrected = raw.copy()
    plate = wells["plate"].to_numpy()
    for name in sorted(set(plate)):
        test = plate == name
        fit = fit_compound_adjusted_tilt(
            wells[~test], features, min_cells, include_controls=True
        )
        corrected[test] = raw[test] - pcu.tilt(
            fit["coef"], fit["center"], rows[test], cols[test]
        )
    return corrected


def position_variance_parts(
    values: np.ndarray, wells: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    """Measure the variation within groups of wells that is position, per feature.

    The biology of a group of wells (the healthy DMSO, or the failing DMSO) is removed by
    looking only at how its own wells differ. A group sits at several positions, the same on
    every replicate plate. The variation of the position means of a group (which repeats on
    all plates) is position, plus any interaction of the group with its position. The
    variation of the replicate plates around their position mean is noise.

    A feature with no position effect still has some variation between positions by chance,
    in proportion to the noise (the degrees of freedom of the positions over those of the
    replicates), so the function subtracts that chance level. What remains is the position
    variation (the ``excess``), and the total is the variation within the groups.

    Parameters
    ----------
    values : np.ndarray
        Profiles of the compound wells with shape ``(len(wells), n_features)``.
    wells : pd.DataFrame
        The wells with the columns ``group`` (healthy or failing DMSO) and ``well``, in the order of
        ``values``.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The position variation above chance (the ``excess``) and the total variation within
        the compounds, for every feature, in squared standardized units summed over the
        wells. The share of the variation that is position is the ratio of the two.
    """
    compound = wells["group"].to_numpy()
    position = (wells["group"] + "_" + wells["well"]).to_numpy()
    frame = pd.DataFrame(values)
    compound_mean = frame.groupby(compound).transform("mean").to_numpy()
    position_mean = frame.groupby(position).transform("mean").to_numpy()
    between = np.zeros(values.shape[1])
    within = np.zeros(values.shape[1])
    chance = np.zeros(values.shape[1])
    for name in sorted(set(compound)):
        in_compound = compound == name
        n_positions = len(set(position[in_compound]))
        n_wells = int(in_compound.sum())
        noise = ((values[in_compound] - position_mean[in_compound]) ** 2).sum(axis=0)
        between += (
            (position_mean[in_compound] - compound_mean[in_compound]) ** 2
        ).sum(axis=0)
        within += noise
        if n_wells <= n_positions:
            # one well per position leaves no replicate to estimate the chance level from
            chance += np.nan
        else:
            chance += noise * (n_positions - 1) / (n_wells - n_positions)
    return between - chance, between + within


def position_variance_share(values: np.ndarray, wells: pd.DataFrame) -> np.ndarray:
    """Measure the share of the variation within compounds that is position, per feature.

    See ``position_variance_parts``. A share of 0 means no more variation between positions
    than chance.

    Parameters
    ----------
    values : np.ndarray
        Profiles of the compound wells with shape ``(len(wells), n_features)``.
    wells : pd.DataFrame
        The wells with the columns ``group`` (healthy or failing DMSO) and ``well``, in the order of
        ``values``.

    Returns
    -------
    np.ndarray
        The position share of every feature (missing for features without variation).
    """
    excess, total = position_variance_parts(values, wells)
    with np.errstate(divide="ignore", invalid="ignore"):
        return excess / total


def held_out_position_variance(
    well_medians: pd.DataFrame,
    features: Sequence[str],
    min_cells: int = 50,
    amplitudes: dict[int, float] | None = None,
) -> pd.DataFrame:
    """Measure the position variation of every feature before and after the correction.

    The measurement uses the DMSO wells, which sit at the same positions on all five
    plates, so the replicate plates give the noise (see ``position_variance_parts``): the
    healthy DMSO and the failing DMSO are two groups.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Median profile of every well (see ``load_or_build_well_medians``).
    features : Sequence[str]
        Features to measure.
    min_cells : int, optional
        Minimum cells of a well, by default 50.
    amplitudes : dict[int, float] | None, optional
        Amplitude of each plate (plate number mapped to the shrunk amplitude of the fit).
        The in-sample correction multiplies the tilt of a plate by its amplitude, as step 4b
        does. The held-out correction cannot know the amplitude of a plate that it never saw,
        and uses the tilt as it is (amplitude 1). By default, every amplitude is 1.

    Returns
    -------
    pd.DataFrame
        A table with one row per feature and the position share (``position_variance_share``) of
        the uncorrected wells (``raw``), of the wells corrected by a tilt that was fit on all
        plates (``in_sample``, which a fit on its own data always favors), and of the wells
        corrected by a tilt that was fit on the other four plates (``held_out``). The
        columns ending in ``_excess`` are the position variation above chance in absolute
        units (``position_variance_parts``).
    """
    features = list(features)
    wells = well_medians[well_medians["n"] >= min_cells].reset_index(drop=True)
    wells["group"] = wells["treatment"] + "|" + wells["cell_type"]
    rows, cols = pcu.well_index(wells["well"])
    raw = wells[features].to_numpy(float)
    fit = fit_compound_adjusted_tilt(wells, features, min_cells, include_controls=True)
    scale = (
        wells["platemap"].map(amplitudes).to_numpy(float)[:, None]
        if amplitudes is not None
        else 1.0
    )
    in_sample = raw - scale * pcu.tilt(fit["coef"], fit["center"], rows, cols)
    held_out = leave_one_plate_out_correction(wells, features, min_cells)
    is_control = (wells["treatment"] == "DMSO").to_numpy()
    controls = wells[is_control].reset_index(drop=True)
    table = {"feature": features}
    for name, values in (("raw", raw), ("in_sample", in_sample), ("held_out", held_out)):
        excess, total = position_variance_parts(values[is_control], controls)
        with np.errstate(divide="ignore", invalid="ignore"):
            table[name] = excess / total
        table[f"{name}_excess"] = excess
    return pd.DataFrame(table)


def well_score_table(
    medians: pd.DataFrame, before: np.ndarray, after: np.ndarray
) -> pd.DataFrame:
    """Combine the principal component scores of wells before and after correction.

    This is ``position_correction_effect_utils.well_score_table`` with the plate number
    (the ``platemap`` column of ``medians``, which is the number of the plate map of the
    plate) in the column ``plate_number``, because the validation plates share one layout.

    Parameters
    ----------
    medians : pd.DataFrame
        One row per well, with the columns ``plate``, ``platemap``, ``well``, ``treatment``
        and ``cell_type``.
    before, after : np.ndarray
        Wells by components scores before and after the correction, in the same principal
        components and in the row order of ``medians``.

    Returns
    -------
    pd.DataFrame
        Two rows per well (``version`` is ``"Before correction"`` or ``"After
        correction"``) with the columns of
        ``position_correction_effect_utils.well_score_table``, and ``plate_number`` in
        place of ``platemap_number``.
    """
    table = pce.well_score_table(medians, before, after)
    # the plate number of the plate map, which does not depend on the plates in the table
    table.insert(2, "plate_number", table["platemap_number"])
    return table.drop(columns="platemap_number")


def _row_probabilities(
    train_x: np.ndarray,
    train_rows: np.ndarray,
    test_x: np.ndarray,
    n_components: int,
    seed: int,
) -> np.ndarray:
    """Train a row classifier and predict the row probabilities of new wells.

    The model is the one of the screen's held-out validation
    (``3a.position_correction/1.validate_position_correction.ipynb``): standardize the
    profiles, reduce them to ``n_components`` principal components, and fit a
    class-balanced multinomial logistic regression of the plate row.

    Parameters
    ----------
    train_x : np.ndarray
        Training wells by features.
    train_rows : np.ndarray
        Row index (0 to 5, rows B to G) of each training well.
    test_x : np.ndarray
        Wells to predict, with the same features as ``train_x``.
    n_components : int
        Number of principal components.
    seed : int
        Random seed of the PCA.

    Returns
    -------
    np.ndarray
        Probabilities with shape ``(len(test_x), N_ROWS)``; a row that is absent from
        ``train_rows`` gets probability 0.
    """
    model = make_pipeline(
        StandardScaler(),
        PCA(n_components=n_components, random_state=seed),
        LogisticRegression(max_iter=2000, class_weight="balanced"),
    )
    model.fit(train_x, train_rows)
    probabilities = np.zeros((len(test_x), N_ROWS))
    probabilities[:, model.classes_] = model.predict_proba(test_x)
    return probabilities


def leave_one_plate_out_row_prediction(
    wells: pd.DataFrame,
    features: Sequence[str],
    make_offsets: Callable[
        [pd.DataFrame], Callable[[np.ndarray, np.ndarray], np.ndarray]
    ]
    | None = None,
    shuffle_seed: int | None = None,
    n_components: int = 50,
    seed: int = 0,
) -> np.ndarray:
    """Predict the plate row of wells of a plate that the correction and the classifier never saw.

    The question is the one of the screen: after a correction, can a classifier still tell
    which plate row a well came from?
    The screen holds out a platemap. The validation plates share the positions of the DMSO
    wells and differ in their compounds, so we hold out one plate at a time. For every plate,
    the function

    1. fits the correction on the other plates only (``make_offsets``), or leaves the
       profiles as they are when ``make_offsets`` is None;
    2. corrects all wells with it;
    3. trains the row classifier on the wells of the other plates;
    4. predicts the row probabilities of the wells of the held-out plate.

    Parameters
    ----------
    wells : pd.DataFrame
        Wells with the columns ``plate``, ``well`` and the features.
    features : Sequence[str]
        Features of the classifier.
    make_offsets : Callable | None, optional
        Function that takes the training wells (a data frame like ``wells``) and returns a
        function of the zero-based row and column indices of wells that gives the offsets
        to subtract, with one row per well and one column per feature.
        By default, there is no correction.
    shuffle_seed : int | None, optional
        When set, the row labels of the training wells are shuffled with this seed.
        This gives the chance-level floor, as in the screen.
    n_components : int, optional
        Number of principal components of the classifier, by default 50.
    seed : int, optional
        Random seed of the PCA, by default 0.

    Returns
    -------
    np.ndarray
        Row probabilities with shape ``(len(wells), N_ROWS)``, in the order of ``wells``.
    """
    features = list(features)
    wells = wells.reset_index(drop=True)
    rows, cols = pcu.well_index(wells["well"])
    profiles = wells[features].to_numpy(float)
    plate = wells["plate"].to_numpy()
    rng = np.random.default_rng(shuffle_seed)
    probabilities = np.zeros((len(wells), N_ROWS))
    for name in sorted(set(plate)):
        test = plate == name
        train = ~test
        corrected = profiles
        if make_offsets is not None:
            corrected = profiles - make_offsets(wells[train])(rows, cols)
        train_rows = rows[train]
        if shuffle_seed is not None:
            train_rows = rng.permutation(train_rows)
        probabilities[test] = _row_probabilities(
            corrected[train], train_rows, corrected[test], n_components, seed
        )
    return probabilities


def row_prediction_scores(rows: np.ndarray, probabilities: np.ndarray) -> dict:
    """Score row predictions as in the screen.

    Parameters
    ----------
    rows : np.ndarray
        True row index (0 to 5) of each well.
    probabilities : np.ndarray
        Row probabilities with shape ``(len(rows), N_ROWS)``.

    Returns
    -------
    dict
        ``macro_auc``, the one-vs-rest AUC averaged over the rows that the wells are in
        (the probabilities of the other rows are left out and the rest are rescaled to add
        up to 1), and ``row_b_auc``, the AUC of the row B probability for row B wells against
        all other wells (missing when no well is in row B, or all are).
        A value of 0.5 is chance level.
    """
    present = sorted(np.unique(rows))
    present_probabilities = probabilities[:, present]
    present_probabilities = present_probabilities / present_probabilities.sum(
        axis=1, keepdims=True
    )
    is_row_b = rows == 0
    if len(present) < 2:
        macro_auc = float("nan")  # one row has nothing to be told apart from
    elif len(present) == 2:
        # two rows: scikit-learn wants the binary form, and the two one-vs-rest AUCs are equal
        macro_auc = float(roc_auc_score(rows == present[1], present_probabilities[:, 1]))
    else:
        macro_auc = float(
            roc_auc_score(
                rows,
                present_probabilities,
                multi_class="ovr",
                average="macro",
                labels=present,
            )
        )
    return {
        "macro_auc": macro_auc,
        "row_b_auc": float(roc_auc_score(is_row_b, probabilities[:, 0]))
        if 0 < is_row_b.sum() < len(rows)
        else float("nan"),
    }


def control_contrast(
    well_medians: pd.DataFrame,
    features: Sequence[str],
    offsets: np.ndarray | None = None,
) -> pd.DataFrame:
    """Measure the healthy minus failing DMSO contrast of every plate.

    The healthy and the failing DMSO wells are selected by ``treatment`` and ``cell_type``
    (not by row). Both sit in all six rows, in blocks of columns (failing in columns 2, 4,
    and 10, healthy in 3, 9, and 11), so a tilt of rows and columns shifts them alike.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians of the validation plates, with the columns ``plate``, ``treatment``,
        ``cell_type`` and the features.
    features : Sequence[str]
        Features of the contrast.
    offsets : np.ndarray | None, optional
        Offsets to subtract from the wells, with one row per row of ``well_medians`` and
        one column per feature.
        By default, there is no correction.

    Returns
    -------
    pd.DataFrame
        One row per plate with the column ``norm``, the Euclidean length of the difference
        between the mean healthy DMSO and the mean failing DMSO profile.
    """
    features = list(features)
    profiles = well_medians[features].to_numpy(float)
    if offsets is not None:
        profiles = profiles - offsets
    is_dmso = (well_medians["treatment"] == "DMSO").to_numpy()
    healthy = is_dmso & (well_medians["cell_type"] == HEALTHY_CELL_TYPE).to_numpy()
    failing = is_dmso & (well_medians["cell_type"] == FAILING_CELL_TYPE).to_numpy()
    records = []
    for plate in well_medians["plate"].unique():
        in_plate = (well_medians["plate"] == plate).to_numpy()
        difference = profiles[in_plate & healthy].mean(axis=0) - profiles[
            in_plate & failing
        ].mean(axis=0)
        records.append({"plate": plate, "norm": float(np.linalg.norm(difference))})
    return pd.DataFrame(records)


def training_tilt_offsets(
    training_wells: pd.DataFrame, features: Sequence[str], min_cells: int = 50
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """Fit the tilt on the training plates only and return its offsets of well positions.

    The result is a function of the row and column indices of wells, as the
    ``make_offsets`` argument of ``leave_one_plate_out_row_prediction`` expects.
    """
    fit = fit_compound_adjusted_tilt(
        training_wells, features, min_cells, include_controls=True
    )
    return lambda rows, columns: pcu.tilt(fit["coef"], fit["center"], rows, columns)


def held_out_row_prediction_scores(
    rows: np.ndarray,
    probabilities: dict[str, np.ndarray],
    shuffled: Sequence[np.ndarray],
    subset: np.ndarray,
) -> tuple[dict, dict, pd.DataFrame]:
    """Pooled held-out scores of the raw, corrected, and shuffled classifiers on a subset of wells.

    ``probabilities`` has the keys ``"raw"`` and ``"corrected"``, and ``shuffled`` is a list
    of the probabilities of classifiers that were trained on shuffled rows.
    """
    raw = row_prediction_scores(rows[subset], probabilities["raw"][subset])
    corrected = row_prediction_scores(rows[subset], probabilities["corrected"][subset])
    floor = pd.DataFrame([row_prediction_scores(rows[subset], p[subset]) for p in shuffled])
    return raw, corrected, floor


def tilt_recovery(
    estimated: dict, true_tilt: np.ndarray, rows: np.ndarray, columns: np.ndarray
) -> tuple[float, float]:
    """Median cosine and relative size of an estimated tilt, compared with the true tilt (both centered)."""
    tilt_estimate = pcu.tilt(estimated["coef"], estimated["center"], rows, columns)
    tilt_estimate = tilt_estimate - tilt_estimate.mean(axis=0)
    centered_truth = true_tilt - true_tilt.mean(axis=0)
    cosines = np.sum(tilt_estimate * centered_truth, axis=1) / (
        np.linalg.norm(tilt_estimate, axis=1) * np.linalg.norm(centered_truth, axis=1)
    )
    return float(np.median(cosines)), float(pcu.amplitude(tilt_estimate, centered_truth))


def group_effect_recovery(
    wells: pd.DataFrame,
    features: Sequence[str],
    estimated: dict,
    true_effects: dict[str, np.ndarray],
) -> tuple[float, float]:
    """Compare the group effects after the correction with the true group effects of a simulation.

    The tilt that ``estimated`` holds is subtracted from the simulated wells, the corrected
    profiles are averaged within every group (a compound, or the DMSO of one cell type) and
    centered over the groups, and the result is compared with the true group effects
    (centered the same way).

    Parameters
    ----------
    wells : pd.DataFrame
        Simulated wells with the columns ``treatment``, ``cell_type``, ``well`` and the features.
    features : Sequence[str]
        Names of the feature columns.
    estimated : dict
        Tilt map with the keys ``coef`` and ``center``.
    true_effects : dict[str, np.ndarray]
        True effect of every group over the features, keyed ``"<treatment>|<cell_type>"``.

    Returns
    -------
    tuple[float, float]
        The median over the groups of the cosine between the corrected and the true group
        effect, and the scale of the corrected effects along the true ones (1 if the
        correction keeps the group effects at their size).
    """
    rows, cols = pcu.well_index(wells["well"])
    corrected = wells[list(features)].to_numpy(float) - pcu.tilt(
        estimated["coef"], estimated["center"], rows, cols
    )
    keys = (wells["treatment"] + "|" + wells["cell_type"]).to_numpy()
    groups = sorted(true_effects)
    corrected_means = np.vstack([corrected[keys == group].mean(axis=0) for group in groups])
    corrected_means -= corrected_means.mean(axis=0)
    truth = np.vstack([true_effects[group] for group in groups])
    truth -= truth.mean(axis=0)
    cosines = np.sum(corrected_means * truth, axis=1) / (
        np.linalg.norm(corrected_means, axis=1) * np.linalg.norm(truth, axis=1)
    )
    return float(np.median(cosines)), float((corrected_means * truth).sum() / (truth * truth).sum())


def mean_row_pattern(
    wells: pd.DataFrame, features: Sequence[str]
) -> tuple[list[str], np.ndarray]:
    """Mean profile of every plate row (rows B to G), centered over the rows.

    Parameters
    ----------
    wells : pd.DataFrame
        Well medians with the column ``well`` (such as ``"B02"``) and the features.
    features : Sequence[str]
        Names of the feature columns.

    Returns
    -------
    tuple[list[str], np.ndarray]
        The row labels and the rows by features pattern, each feature centered over the rows.

    Raises
    ------
    ValueError
        If a row has no well.
    """
    means = (
        wells.assign(row=wells["well"].str[0])
        .groupby("row")[list(features)]
        .mean()
        .reindex(list(pcu.ROWS))
    )
    if means.isna().any().any():
        raise ValueError("A plate row has no well, so the row pattern is not defined")
    return list(means.index), means.to_numpy(float) - means.to_numpy(float).mean(axis=0)


def tilt_row_pattern(fit: dict) -> np.ndarray:
    """Row effects of a tilt (rows B to G by features), each feature centered over the rows."""
    rows = fit["coef"][: pcu.N_ROW_TERMS]
    return rows - rows.mean(axis=0)


def cosine_and_slope(observed: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    """Cosine and slope between two patterns (positions by features)."""
    inner = float((observed * predicted).sum())
    return (
        inner / np.sqrt((observed**2).sum() * (predicted**2).sum()),
        inner / float((predicted**2).sum()),
    )


def distance_from_class_mean(
    values: np.ndarray, is_healthy: np.ndarray, is_failing: np.ndarray
) -> np.ndarray:
    """Distance of each DMSO well from the mean profile of its cell type (NaN for the other wells)."""
    result = np.full(len(values), np.nan)
    for mask in (is_healthy, is_failing):
        result[mask] = np.linalg.norm(values[mask] - values[mask].mean(axis=0), axis=1)
    return result


def tilt_position_terms(fit: dict) -> np.ndarray:
    """Row effects against row B and column effects against column 2: one vector over features per term."""
    coefficients = fit["coef"]
    rows = coefficients[1 : pcu.N_ROW_TERMS] - coefficients[0]
    cols = coefficients[pcu.N_ROW_TERMS :]
    return np.vstack([rows, cols])


def plate_layout_grid(
    table: pd.DataFrame,
    plate: str,
    value_column: str,
    row_labels: Sequence[str],
    column_numbers: Sequence[int],
) -> np.ndarray:
    """Values of one plate on the plate layout (rows by columns), NaN where the table has no well.

    ``table`` has the columns ``plate``, ``row`` (a label of ``row_labels``), ``column``
    (a number of ``column_numbers``) and ``value_column``.
    """
    grid = np.full((len(row_labels), len(column_numbers)), np.nan)
    row_index = {label: i for i, label in enumerate(row_labels)}
    column_index = {number: i for i, number in enumerate(column_numbers)}
    for _, well in table[table["plate"] == plate].iterrows():
        grid[row_index[well["row"]], column_index[well["column"]]] = well[value_column]
    return grid
