"""Estimate and apply the slot A adjustment of the validation plates.

The three compounds of every validation plate rotate through the 24 compound wells in the same
pattern on every plate, so a well position always holds the compound of one *slot* (A, B, or C).
The compounds were dispensed into the slots at random, and the compounds of slot A score lower
than those of slots B and C on every plate. A position effect that is a function of the slot
cannot be told from the effect of the compound on one plate, so the position correction
(``position_correction_utils``, a row plus column tilt) does not remove it.

The slot A adjustment estimates what separates the compounds of slot A from those of slots B and C
on all plates, and removes a shrunk part of it from the cells of the compounds in slot A:

1. For every plate, the profile of each compound is the mean of the position-corrected well
   medians of its 8 wells. The effect of slot A on a plate is the profile of the compound in slot A
   minus the mean of the two other compounds, and the estimate ``vector`` is its mean over the
   plates.
2. The estimate also holds the differences that the five compounds of slot A have by chance. Their
   size is found by permuting which compound of each plate is the one of slot A (``3 ** n_plates``
   arrangements, 243 for five plates): the mean squared size of the vector (its ``energy``) is
   compared with the energies of the arrangements.
3. The share of the observed energy that random assignment does not give, ``factor = (energy -
   median of the arrangements) / energy``, is the part of the estimate that we treat as the slot
   effect. The adjustment subtracts ``factor * vector`` from every cell of a compound in a slot A well.

The adjustment is a separate step that does not change the position-corrected profiles.
"""

import itertools
import json
import pathlib
from collections.abc import Sequence

import numpy as np
import pandas as pd
import position_correction_effect_utils as pce

SLOT_WELLS = {
    "A": ("B05", "B08", "C07", "D06", "E05", "E08", "F06", "G07"),
    "B": ("B06", "C05", "C08", "D07", "E06", "F07", "G05", "G08"),
    "C": ("B07", "C06", "D05", "D08", "E07", "F05", "F08", "G06"),
}
SLOT_OF_WELL = {well: slot for slot, wells in SLOT_WELLS.items() for well in wells}
CONTROL = "DMSO"
MAX_PLATES = 8  # 3 ** 8 arrangements are enumerated


def slot_of_wells(wells: pd.Series) -> pd.Series:
    """Slot (A, B, or C) of well positions such as ``"B05"``; missing for the other wells."""
    return wells.map(SLOT_OF_WELL)


def compound_slot_table(
    well_medians: pd.DataFrame, fit: dict, features: Sequence[str], min_cells: int = 50
) -> pd.DataFrame:
    """Profile of every compound after the position correction, with its plate and slot.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians before the correction, with the columns ``plate``, ``platemap``, ``well``,
        ``treatment``, ``n`` and the features (``validation_position_utils.load_or_build_well_medians``).
    fit : dict
        Position correction fit (``position_correction_utils.load_fit``).
    features : Sequence[str]
        Features of the table, which have to be features of the fit.
    min_cells : int, optional
        Minimum number of cells of a well, by default 50.

    Returns
    -------
    pd.DataFrame
        One row per compound with the columns ``plate``, ``treatment``, ``slot`` and the features
        (the mean over the wells of the compound of the corrected well medians).

    Raises
    ------
    ValueError
        If a compound has wells in more than one slot, or a plate does not have one compound
        in each slot.
    """
    wells = well_medians[well_medians["n"] >= min_cells]
    corrected = pce.correct_well_medians(wells, fit)
    compounds = corrected[corrected["treatment"] != CONTROL].copy()
    compounds["slot"] = slot_of_wells(compounds["well"])
    if compounds["slot"].isna().any():
        raise ValueError("A compound has a well that is not in a slot")
    slots = compounds.groupby(["plate", "treatment"])["slot"].nunique()
    if (slots != 1).any():
        raise ValueError(f"Compounds with wells in several slots: {slots[slots != 1].index.tolist()}")
    table = (
        compounds.groupby(["plate", "treatment", "slot"])[list(features)].mean().reset_index()
    )
    per_plate = table.groupby("plate")["slot"].apply(lambda s: "".join(sorted(s)))
    if (per_plate != "ABC").any():
        raise ValueError(f"Plates without one compound in each slot: {per_plate[per_plate != 'ABC'].index.tolist()}")
    return table


def _plate_profiles(table: pd.DataFrame, features: Sequence[str]) -> np.ndarray:
    """Array of shape ``(n_plates, 3, n_features)`` with the compounds of slots A, B, and C."""
    plates = sorted(table["plate"].unique())
    stacked = [
        np.stack([table[(table["plate"] == plate) & (table["slot"] == slot)][list(features)].to_numpy(float)[0] for slot in "ABC"])
        for plate in plates
    ]
    return np.stack(stacked)


def slot_a_vector(profiles: np.ndarray, choice: Sequence[int] | None = None) -> np.ndarray:
    """Mean over the plates of the compound in slot A minus the mean of the other two.

    Parameters
    ----------
    profiles : np.ndarray
        Output of ``_plate_profiles``, shape ``(n_plates, 3, n_features)``.
    choice : Sequence[int], optional
        For every plate, the index (0, 1, or 2) of the compound that counts as the one of slot A;
        by default 0, the compound that is in slot A.

    Returns
    -------
    np.ndarray
        The estimate of the effect of slot A, one value for each feature.
    """
    choice = [0] * len(profiles) if choice is None else list(choice)
    differences = []
    for plate, k in zip(profiles, choice, strict=True):
        others = [i for i in range(3) if i != k]
        differences.append(plate[k] - plate[others].mean(axis=0))
    return np.mean(differences, axis=0)


def _estimate_from_profiles(profiles: np.ndarray) -> dict:
    """Vector, energy, permutation null, p-value and shrinkage factor of the plates of ``profiles``."""
    n_plates = len(profiles)
    if n_plates > MAX_PLATES:
        raise ValueError(f"{n_plates} plates are too many to enumerate the arrangements (maximum {MAX_PLATES})")
    vector = slot_a_vector(profiles)
    energy = float((vector**2).mean())
    null = np.array(
        [
            float((slot_a_vector(profiles, choice) ** 2).mean())
            for choice in itertools.product(range(3), repeat=n_plates)
        ]
    )
    median = float(np.median(null))
    return {
        "vector": vector,
        "energy": energy,
        "null_energies": null,
        "null_median": median,
        "p_value": float((null >= energy - 1e-12).mean()),
        "factor": max(0.0, (energy - median) / energy) if energy > 0 else 0.0,
    }


def estimate_slot_a_adjustment(
    well_medians: pd.DataFrame, fit: dict, features: Sequence[str] | None = None, min_cells: int = 50
) -> dict:
    """Estimate the slot A vector, its permutation null, and the shrinkage factor.

    The estimate is also made without each plate in turn (a jackknife), to show how much it depends
    on a single plate and to carry the uncertainty of the estimate into the hit calling.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Well medians before the correction (see ``compound_slot_table``).
    fit : dict
        Position correction fit (``position_correction_utils.load_fit``).
    features : Sequence[str], optional
        Features to estimate, by default the features of the fit that are in ``well_medians``.
    min_cells : int, optional
        Minimum number of cells of a well, by default 50.

    Returns
    -------
    dict
        ``features`` and ``vector`` (the full estimate), ``energy`` (its mean squared size),
        ``null_energies`` (the energy of every arrangement), ``null_median``, ``p_value`` (the share of
        arrangements with an energy at least as large), ``factor`` (the shrinkage), ``plates``,
        ``compounds`` (the table of ``compound_slot_table`` without the features), ``min_cells``, and the
        estimate without each plate: ``jackknife_vectors`` (one row for each plate, in the order of
        ``plates``), ``jackknife_factors``, ``jackknife_energies``, and ``jackknife_p_values``.
    """
    if features is None:
        features = [f for f in fit["features"] if f in well_medians.columns]
    features = list(features)
    table = compound_slot_table(well_medians, fit, features, min_cells)
    profiles = _plate_profiles(table, features)
    full = _estimate_from_profiles(profiles)
    jackknife = [_estimate_from_profiles(np.delete(profiles, i, axis=0)) for i in range(len(profiles))]
    return {
        "features": features,
        **full,
        "plates": sorted(table["plate"].unique()),
        "compounds": table[["plate", "treatment", "slot"]],
        "min_cells": min_cells,
        "jackknife_vectors": np.stack([j["vector"] for j in jackknife]),
        "jackknife_factors": np.array([j["factor"] for j in jackknife]),
        "jackknife_energies": np.array([j["energy"] for j in jackknife]),
        "jackknife_p_values": np.array([j["p_value"] for j in jackknife]),
    }


def save_slot_adjustment(adjustment: dict, path: str | pathlib.Path) -> None:
    """Save the slot A adjustment to a compressed ``.npz`` file (``load_slot_adjustment`` reads it)."""
    np.savez_compressed(
        path,
        features=np.array(adjustment["features"]),
        vector=adjustment["vector"],
        energy=adjustment["energy"],
        null_energies=adjustment["null_energies"],
        null_median=adjustment["null_median"],
        p_value=adjustment["p_value"],
        factor=adjustment["factor"],
        min_cells=adjustment["min_cells"],
        plates=np.array(adjustment["plates"]),
        compounds=json.dumps(adjustment["compounds"].to_dict(orient="records")),
        jackknife_vectors=adjustment["jackknife_vectors"],
        jackknife_factors=adjustment["jackknife_factors"],
        jackknife_energies=adjustment["jackknife_energies"],
        jackknife_p_values=adjustment["jackknife_p_values"],
        slot_wells=json.dumps({slot: list(wells) for slot, wells in SLOT_WELLS.items()}),
    )


def load_slot_adjustment(path: str | pathlib.Path) -> dict:
    """Load a slot A adjustment saved by ``save_slot_adjustment``."""
    with np.load(path, allow_pickle=False) as data:
        return {
            "features": data["features"].tolist(),
            "vector": data["vector"],
            "energy": float(data["energy"]),
            "null_energies": data["null_energies"],
            "null_median": float(data["null_median"]),
            "p_value": float(data["p_value"]),
            "factor": float(data["factor"]),
            "min_cells": int(data["min_cells"]),
            "plates": data["plates"].tolist(),
            "compounds": pd.DataFrame(json.loads(str(data["compounds"]))),
            "jackknife_vectors": data["jackknife_vectors"],
            "jackknife_factors": data["jackknife_factors"],
            "jackknife_energies": data["jackknife_energies"],
            "jackknife_p_values": data["jackknife_p_values"],
        }


def apply_slot_a_adjustment(
    profiles: pd.DataFrame, adjustment: dict, factor: float | None = None
) -> pd.DataFrame:
    """Subtract the shrunk slot A vector from the cells (or wells) of the compounds in slot A wells.

    The controls and the compounds in slots B and C are not changed, and neither is any feature
    that is not in the adjustment. The result has two more columns: ``Metadata_slot`` (the slot of the well
    of a compound, missing for the controls) and ``Metadata_slot_a_adjusted`` (whether the
    row was shifted), so that the adjusted profiles are marked as such.

    Parameters
    ----------
    profiles : pd.DataFrame
        Position-corrected single-cell or bulk profiles, with the columns ``Metadata_Well``,
        ``Metadata_treatment`` and the features.
    adjustment : dict
        Output of ``estimate_slot_a_adjustment`` or ``load_slot_adjustment``.
    factor : float, optional
        Share of the vector to subtract, by default the shrinkage factor of the adjustment
        (1 subtracts the full estimate).

    Returns
    -------
    pd.DataFrame
        Copy of ``profiles`` with the adjusted features and the two marker columns.
    """
    factor = adjustment["factor"] if factor is None else factor
    result = profiles.copy()
    is_compound = (result["Metadata_treatment"] != CONTROL).to_numpy()
    slot = slot_of_wells(result["Metadata_Well"])
    shifted = ((slot == "A").to_numpy()) & is_compound
    features = [f for f in adjustment["features"] if f in result.columns]
    vector = pd.Series(adjustment["vector"], index=adjustment["features"])[features].to_numpy()
    values = result[features].to_numpy(dtype=float)
    values[shifted] -= factor * vector
    result[features] = values
    result["Metadata_slot"] = np.where(is_compound, slot, None)
    result["Metadata_slot_a_adjusted"] = shifted
    return result
