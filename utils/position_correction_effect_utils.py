"""Study how the plate-position correction changes the feature space.

``position_correction_utils`` estimates the plate-position effect ("tilt") and corrects
single cells with it. The correction of a well is ``amplitude * tilt(row, column)``, and
it shifts every cell of the well by the same vector. The median profile of a well after
the correction is therefore the median profile before the correction minus that vector.
The functions here use this to compare wells before and after the correction, without
reading single cells.
"""

from collections.abc import Sequence

import numpy as np
import pandas as pd
import position_correction_utils as pcu
from sklearn.decomposition import PCA

CHANNELS = ("Actin", "DNA", "ER", "Mito", "PM")
MULTI_CHANNEL = "Multiple channels"
NO_CHANNEL = "No channel"


def correct_well_medians(medians: pd.DataFrame, fit: dict) -> pd.DataFrame:
    """Apply the position correction to well median profiles.

    Parameters
    ----------
    medians : pd.DataFrame
        One row per well, with the columns ``platemap`` (number), ``well`` and the
        features of the fit.
    fit : dict
        Position correction fit (see ``position_correction_utils.load_fit``).

    Returns
    -------
    pd.DataFrame
        Copy of ``medians`` where each feature of the fit has the correction
        ``amplitude * tilt(row, column)`` of the platemap of the well subtracted.

    Raises
    ------
    KeyError
        If a platemap of the wells has no amplitude in the fit.
    """
    rows, cols = pcu.well_index(medians["well"])
    amplitudes = medians["platemap"].map(fit["amplitudes"])
    if amplitudes.isna().any():
        missing = sorted(medians.loc[amplitudes.isna(), "platemap"].unique())
        raise KeyError(f"No amplitude in the fit for platemaps {missing}")
    offsets = amplitudes.to_numpy()[:, None] * pcu.tilt(
        fit["coef"], fit["center"], rows, cols
    )
    corrected = medians.copy()
    corrected[fit["features"]] = (
        medians[fit["features"]].to_numpy(dtype=float) - offsets
    )
    return corrected


def parse_feature_family(feature: str) -> dict[str, str]:
    """Split a CellProfiler feature name into compartment, measurement and channel.

    Parameters
    ----------
    feature : str
        Feature name such as ``"Cells_Texture_SumAverage_DNA_3_00_256"``. The first
        part is the compartment and the second part is the measurement type.

    Returns
    -------
    dict[str, str]
        A dictionary with the keys:

        - "compartment": for example ``"Cells"``
        - "measurement": for example ``"Texture"``
        - "channel": the channel (``"DNA"``), ``MULTI_CHANNEL`` for features that
          combine channels (such as correlations), or ``NO_CHANNEL`` for features
          without a channel (such as area and shape)
    """
    parts = feature.split("_")
    channels = []
    for part in parts[2:]:
        if part in CHANNELS and part not in channels:
            channels.append(part)
    if not channels:
        channel = NO_CHANNEL
    elif len(channels) > 1:
        channel = MULTI_CHANNEL
    else:
        channel = channels[0]
    return {"compartment": parts[0], "measurement": parts[1], "channel": channel}

def tilt_size_by_feature(fit: dict) -> pd.DataFrame:
    """Measure the size of the tilt of every feature and label its family.

    The size of a feature is the root mean square of its tilt over all positions of the
    plate (6 rows by 10 columns), with an amplitude of 1. The features are z-scored on
    each plate, so the size is in units of the standard deviation of a single cell.

    Parameters
    ----------
    fit : dict
        Position correction fit (see ``position_correction_utils.load_fit``).

    Returns
    -------
    pd.DataFrame
        One row per feature with the columns ``feature``, ``tilt_rms``,
        ``compartment``, ``measurement`` and ``channel``.
    """
    rows, cols = np.meshgrid(
        np.arange(pcu.N_ROW_TERMS), np.arange(pcu.N_COL_TERMS + 1), indexing="ij"
    )
    tilt = pcu.tilt(fit["coef"], fit["center"], rows.ravel(), cols.ravel())
    sizes = pd.DataFrame(
        {
            "feature": fit["features"],
            "tilt_rms": np.sqrt((tilt**2).mean(axis=0)),
        }
    )
    families = pd.DataFrame([parse_feature_family(name) for name in sizes["feature"]])
    return pd.concat([sizes, families], axis=1)

def fit_pca(profiles: np.ndarray, n_components: int = 10, random_state: int = 0) -> PCA:
    """Fit a PCA on well profiles without scaling the features.

    The features are z-scored on each plate, so they are already on a similar scale.

    Parameters
    ----------
    profiles : np.ndarray
        Wells by features matrix.
    n_components : int, optional
        Number of principal components, by default 10.
    random_state : int, optional
        Seed of the PCA, by default 0.

    Returns
    -------
    PCA
        Fitted PCA. ``transform`` centers new profiles with the mean of the fitted ones.
    """
    return PCA(n_components=n_components, random_state=random_state).fit(profiles)

def variance_explained_by(values: np.ndarray, groups: Sequence) -> float:
    """Fraction of the variance of a score that the groups explain.

    This is the between-group sum of squares divided by the total sum of squares (eta
    squared) of a one-way analysis of variance.

    Parameters
    ----------
    values : np.ndarray
        Score of each well, for example its first principal component.
    groups : Sequence
        Group of each well, for example its plate row.

    Returns
    -------
    float
        Fraction between 0 and 1. It is 0 when the variance of ``values`` is 0.
    """
    values = np.asarray(values, dtype=float)
    total = ((values - values.mean()) ** 2).sum()
    if total == 0:
        return 0.0
    group_means = pd.Series(values).groupby(np.asarray(groups)).transform("mean")
    return float(((group_means.to_numpy() - values.mean()) ** 2).sum() / total)

def well_score_table(
    medians: pd.DataFrame, before: np.ndarray, after: np.ndarray
) -> pd.DataFrame:
    """Combine the principal component scores of wells before and after correction.

    Parameters
    ----------
    medians : pd.DataFrame
        One row per well, with the columns ``plate``, ``platemap`` (number), ``well``,
        ``treatment`` and ``cell_type``.
    before, after : np.ndarray
        Wells by components scores before and after the correction, in the same
        principal components and in the row order of ``medians``.

    Returns
    -------
    pd.DataFrame
        Two rows per well, one for each ``version`` (``"Before correction"`` and
        ``"After correction"``), with the columns ``plate``, ``platemap_number``,
        ``well``, ``well_row``, ``well_col``, ``treatment``, ``cell_type``,
        ``is_control`` (DMSO wells), ``version`` and one column ``PC1``, ``PC2``, ...
        for each component.
    """
    wells = medians["well"].astype(str).reset_index(drop=True)
    info = pd.DataFrame(
        {
            "plate": medians["plate"].to_numpy(),
            "platemap_number": medians["platemap"].astype(int).to_numpy(),
            "well": wells,
            "well_row": wells.str[0],
            "well_col": wells.str[1:].astype(int),
            "treatment": medians["treatment"].to_numpy(),
            "cell_type": medians["cell_type"].to_numpy(),
            "is_control": (medians["treatment"] == "DMSO").to_numpy(),
        }
    )
    components = [f"PC{i + 1}" for i in range(before.shape[1])]
    tables = [
        pd.concat([info, pd.DataFrame(scores, columns=components)], axis=1).assign(
            version=version
        )
        for version, scores in [
            ("Before correction", before),
            ("After correction", after),
        ]
    ]
    return pd.concat(tables, ignore_index=True)


def r_squared(values: np.ndarray, predictors: np.ndarray) -> float:
    """R squared of a linear model of ``values`` with an intercept and the ``predictors`` columns.

    Parameters
    ----------
    values : np.ndarray
        The response, one value for every well, with shape ``(n_wells,)``.
    predictors : np.ndarray
        The predictors, one column for each, with shape ``(n_wells, n_predictors)``. The intercept is added.

    Returns
    -------
    float
        The share of the variance of ``values`` that the model explains, or 0 when ``values``
        has no variance.
    """
    design = np.column_stack([np.ones(len(values)), predictors])
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    residual = values - design @ coefficients
    total = ((values - values.mean()) ** 2).sum()
    if total == 0:
        return 0.0  # a score without variation has nothing to explain (as variance_explained_by)
    return 1 - (residual**2).sum() / total


def adjusted_r_squared(scores: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Adjusted R squared of every column of ``scores`` on ``design`` (which has its own intercept).

    The residual degrees of freedom use the rank of the design, because the design of
    ``position_correction_utils.design_matrix`` has an indicator for every position, also
    for the positions that the wells do not occupy. The result is NaN for every column
    when there are no more wells than the rank, and for a column without variation.

    Parameters
    ----------
    scores : np.ndarray
        Scores of the wells, with shape ``(n_wells, n_columns)`` (for example principal
        component scores).
    design : np.ndarray
        Design matrix of the wells, with shape ``(n_wells, n_terms)`` and its own intercept.

    Returns
    -------
    np.ndarray
        The adjusted R squared of every column of ``scores``, with shape ``(n_columns,)``.
    """
    n_wells = design.shape[0]
    coefficients, _, rank, _ = np.linalg.lstsq(design, scores, rcond=None)
    if n_wells <= rank:
        return np.full(scores.shape[1], np.nan)
    residual = ((scores - design @ coefficients) ** 2).sum(axis=0)
    total = ((scores - scores.mean(axis=0)) ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1 - (residual / total) * (n_wells - 1) / (n_wells - rank)


def pc_variance_table(
    wells: pd.DataFrame,
    profiles: np.ndarray,
    label: str,
    n_components: int = 10,
    random_state: int = 0,
    factors: Sequence[str] = ("row", "col", "plate", "treatment"),
) -> tuple[pd.DataFrame, PCA, np.ndarray]:
    """Fraction of each principal component that the ``factors`` of ``wells`` explain.

    A factor with a single level is skipped.

    Parameters
    ----------
    wells : pd.DataFrame
        The wells, with the columns that the ``factors`` name, in the order of ``profiles``.
    profiles : np.ndarray
        The profile of every well, with shape ``(n_wells, n_features)``.
    label : str
        Name of the analysis, written to the ``analysis`` column of the table.
    n_components : int, optional
        Number of principal components, by default 10.
    random_state : int, optional
        Seed of the PCA, by default 0.
    factors : Sequence[str], optional
        Columns of ``wells`` to explain the components with, by default ``("row", "col", "plate", "treatment")``.

    Returns
    -------
    tuple[pd.DataFrame, PCA, np.ndarray]
        The table (one row per component, with the variance of the component and the fraction
        that each factor explains), the fitted PCA, and the scores of the wells.
    """
    pca = fit_pca(profiles, n_components, random_state)
    scores = pca.transform(profiles)
    table = pd.DataFrame(
        {
            "analysis": label,
            "component": [f"PC{i + 1}" for i in range(n_components)],
            "variance_explained": pca.explained_variance_ratio_,
        }
    )
    for factor in factors:
        if factor in wells and wells[factor].nunique() > 1:
            table[f"{factor}_explains"] = [
                variance_explained_by(scores[:, i], wells[factor].to_numpy())
                for i in range(n_components)
            ]
    return table, pca, scores


def explained_by_position_and_cells(
    rows: np.ndarray, log_cells: np.ndarray, scores: np.ndarray, label: str
) -> pd.DataFrame:
    """R squared of each component for plate row, log cells per well, and both.

    Parameters
    ----------
    rows : np.ndarray
        Plate row of every well (the row labels, one per well).
    log_cells : np.ndarray
        Logarithm of the number of cells of every well.
    scores : np.ndarray
        Principal component scores of the wells, with shape ``(n_wells, n_components)``.
    label : str
        Name of the analysis, written to the ``analysis`` column of the table.

    Returns
    -------
    pd.DataFrame
        One row per component, with the R squared for the row (``row``), the log number of
        cells (``log_cells``), and both (``row_and_log_cells``).
    """
    present = np.unique(rows)
    row_terms = (rows[:, None] == present[None, 1:]).astype(float)
    records = []
    for i in range(scores.shape[1]):
        records.append(
            {
                "analysis": label,
                "component": f"PC{i + 1}",
                "row": r_squared(scores[:, i], row_terms),
                "log_cells": r_squared(scores[:, i], log_cells[:, None]),
                "row_and_log_cells": r_squared(
                    scores[:, i], np.column_stack([row_terms, log_cells])
                ),
            }
        )
    return pd.DataFrame(records)


def partition_variance(
    medians: pd.DataFrame,
    fit: dict,
    nonfailing_label: str,
    n_components: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Split the variance of the principal components into position and heart status.

    The components are fit on the uncorrected wells and the corrected wells are projected
    on them. The position part of a component is the adjusted R squared of its scores on
    the row and column of the wells, and the status part is the squared correlation with
    the non-failing against failing label.

    The two parts overlap and do not add up: they are associations of the score, not
    variance that one explains and the other does not. On the validation plates the
    non-failing wells are in columns 3, 9, and 11 and no other wells are, so the status label
    is a sum of column indicators and a component that follows only the non-failing columns
    scores 1 on both. The row alone (``row_*``) is not confounded with the status there,
    because every row has the same number of wells of each kind.

    Parameters
    ----------
    medians : pd.DataFrame
        Well medians, with the columns ``well``, ``plate``, ``platemap``, ``cell_type`` and the features.
    fit : dict
        Tilt map and plate amplitudes (``position_correction_utils.load_fit``).
    nonfailing_label : str
        Value of ``cell_type`` of the non-failing cells.
    n_components : int, optional
        Number of principal components, by default 20.
    seed : int, optional
        Seed of the PCA, by default 0.

    Returns
    -------
    pd.DataFrame
        One row per stage (``before`` and ``after`` the correction) and component, with the
        percent of the total variance that the component holds and the parts that position, the
        row alone, and the heart status are associated with, and their R squared.
    """
    features = list(fit["features"])
    raw = medians[features].to_numpy(float)
    corrected = correct_well_medians(medians, fit)[features].to_numpy(float)
    pca = fit_pca(raw, n_components, seed)
    total_variance = pca.explained_variance_[0] / pca.explained_variance_ratio_[0]
    scores = {"before": pca.transform(raw), "after": pca.transform(corrected)}
    rows, cols = pcu.well_index(medians["well"])
    design = pcu.design_matrix(rows, cols)
    row_design = design[:, : pcu.N_ROW_TERMS]  # the row indicators, which act as the intercept
    status = (medians["cell_type"] == nonfailing_label).to_numpy(float)
    records = []
    for stage, values in scores.items():
        variance = values.var(axis=0, ddof=1)
        position = np.clip(adjusted_r_squared(values, design), 0, None)
        row_only = np.clip(adjusted_r_squared(values, row_design), 0, None)
        status_r2 = np.array(
            [np.corrcoef(values[:, k], status)[0, 1] ** 2 for k in range(n_components)]
        )
        for k in range(n_components):
            records.append(
                {
                    "stage": stage,
                    "component": f"PC{k + 1}",
                    "variance_pct_of_total": 100 * variance[k] / total_variance,
                    "position_pct_of_total": 100 * variance[k] * position[k] / total_variance,
                    "row_pct_of_total": 100 * variance[k] * row_only[k] / total_variance,
                    "status_pct_of_total": 100 * variance[k] * status_r2[k] / total_variance,
                    "position_r2": position[k],
                    "row_r2": row_only[k],
                    "status_r2": status_r2[k],
                }
            )
    return pd.DataFrame(records)


def nonfailing_failing_difference(
    medians: pd.DataFrame,
    fit: dict,
    nonfailing_label: str,
    failing_label: str = "failing",
) -> pd.DataFrame:
    """Distance between the non-failing and the failing DMSO profile of every plate, before and after the correction.

    The distance is the Euclidean length of the difference between the mean non-failing DMSO
    profile and the mean failing DMSO profile of a plate, over all features. A plate without a
    non-failing or without a failing DMSO well (after the cell filter of the caller) has NaN
    distances. The columns ``n_nonfailing_dmso_wells`` and ``n_failing_dmso_wells`` say which
    plates they are, and the caller should report them and count the valid plates.

    Parameters
    ----------
    medians : pd.DataFrame
        Well medians, with the columns ``plate``, ``platemap``, ``treatment``, ``cell_type`` and the features.
    fit : dict
        Tilt map and plate amplitudes (``position_correction_utils.load_fit``).
    nonfailing_label : str
        Value of ``cell_type`` of the non-failing cells.
    failing_label : str, optional
        Value of ``cell_type`` of the failing cells, by default ``"failing"``.

    Returns
    -------
    pd.DataFrame
        One row per plate, with the numbers of wells of each kind, the distance before (``raw``)
        and after (``corrected``) the correction, and ``share_that_survives``, the corrected
        distance divided by the uncorrected one.
    """
    features = list(fit["features"])
    raw = medians[features].to_numpy(float)
    corrected = correct_well_medians(medians, fit)[features].to_numpy(float)
    is_dmso = (medians["treatment"] == "DMSO").to_numpy()
    nonfailing = is_dmso & (medians["cell_type"] == nonfailing_label).to_numpy()
    failing = is_dmso & (medians["cell_type"] == failing_label).to_numpy()
    records = []
    for plate in sorted(medians["plate"].unique()):
        in_plate = (medians["plate"] == plate).to_numpy()
        h, f = in_plate & nonfailing, in_plate & failing
        n_nonfailing, n_failing = int(h.sum()), int(f.sum())
        if n_nonfailing and n_failing:
            raw_difference = raw[h].mean(0) - raw[f].mean(0)
            corrected_difference = corrected[h].mean(0) - corrected[f].mean(0)
            raw_length = np.linalg.norm(raw_difference)
            corrected_length = np.linalg.norm(corrected_difference)
        else:  # a plate without both groups has no difference (NaN, not an error)
            raw_length = corrected_length = np.nan
        records.append(
            {
                "plate": plate,
                "platemap": int(medians.loc[in_plate, "platemap"].iloc[0]),
                "n_nonfailing_dmso_wells": n_nonfailing,
                "n_failing_dmso_wells": n_failing,
                "raw": raw_length,
                "corrected": corrected_length,
            }
        )
    table = pd.DataFrame(records)
    table["share_that_survives"] = table["corrected"] / table["raw"]
    return table


def _chi_square_statistic(first: np.ndarray, second: np.ndarray) -> float:
    """Pearson chi-square statistic of the contingency table of two label arrays.

    Parameters
    ----------
    first, second : np.ndarray
        Two arrays of labels of the same length, one entry for each observation.

    Returns
    -------
    float
        The sum over the cells of the table of ``(observed - expected) ** 2 / expected``.
    """
    first_codes, first_levels = pd.factorize(first)
    second_codes, second_levels = pd.factorize(second)
    table = np.zeros((len(first_levels), len(second_levels)))
    np.add.at(table, (first_codes, second_codes), 1)
    expected = table.sum(axis=1, keepdims=True) * table.sum(axis=0, keepdims=True) / table.sum()
    return float(((table - expected) ** 2 / expected).sum())


def within_map_association_test(
    maps: Sequence,
    labels: Sequence,
    positions: Sequence,
    n_permutations: int = 2000,
    seed: int = 0,
) -> tuple[float, float]:
    """Test whether a label (such as the pathway) is associated with the position within platemaps.

    The statistic is the sum over the platemaps of the chi-square statistic of the table of
    the labels against the positions of that platemap, so an association that has a different
    direction on different platemaps does not cancel (as it can in one table of all
    platemaps pooled). The p-value comes from permuting the labels within each platemap,
    which keeps the number of wells of each label and position on every platemap and does not
    rely on the asymptotic distribution of the statistic, which is poor for sparse tables.

    Parameters
    ----------
    maps, labels, positions : Sequence
        One entry per well: its platemap, its label, and its position (such as the plate row).
    n_permutations : int, optional
        Number of permutations, by default 2000.
    seed : int, optional
        Seed of the permutations, by default 0.

    Returns
    -------
    tuple[float, float]
        The statistic and the permutation p-value, ``(1 + permutations at least as large) /
        (1 + n_permutations)``.
    """
    frame = pd.DataFrame(
        {"map": np.asarray(maps), "label": np.asarray(labels), "position": np.asarray(positions)}
    )
    groups = [
        (group["label"].to_numpy(), group["position"].to_numpy())
        for _, group in frame.groupby("map")
    ]

    def statistic(permute: bool) -> float:
        total = 0.0
        for group_labels, group_positions in groups:
            shuffled = rng.permutation(group_labels) if permute else group_labels
            total += _chi_square_statistic(shuffled, group_positions)
        return total

    rng = np.random.default_rng(seed)
    observed = statistic(False)
    at_least = sum(statistic(True) >= observed for _ in range(n_permutations))
    return observed, (1 + at_least) / (1 + n_permutations)
