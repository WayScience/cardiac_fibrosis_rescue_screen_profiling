#!/usr/bin/env python
# coding: utf-8

# # Apply the slot A adjustment to the profiles
# 
# This notebook writes **new** profiles in which the cells of the compounds in slot A wells are shifted by the shrunk slot A estimate of `0.estimate_slot_a_effect.ipynb` (`slot_a_adjustment.npz`).
# The position-corrected profiles of step 4b and 5 are not changed, and the new files have their own names, so that the two kinds of profiles cannot be mixed up:
# 
# | New file | From |
# |---|---|
# | `<plate>_sc_slot_a_adjusted.parquet` | `<plate>_sc_position_corrected.parquet` (all features) |
# | `<plate>_sc_slot_a_adjusted_feature_selected.parquet` | the same features as `<plate>_sc_position_corrected_feature_selected.parquet` |
# | `bulk_slot_a_adjusted_feature_selected.parquet` | `bulk_position_corrected_feature_selected.parquet` |
# 
# Every cell or well gets two more columns: `Metadata_slot` (A, B, or C for the compound wells, missing for the controls) and `Metadata_slot_a_adjusted` (true for the cells and wells that were shifted).
# Controls and the compounds in slots B and C are not changed, and the position tilt is not touched.
# 
# - **Same features.** The feature-selected files keep the features of step 4b, so that the only difference to the position-corrected profiles is the shift. Feature selection is not run again.
# - **Bulk profiles.** The median of cells that all shift by the same vector is the median shifted by that vector, so the bulk profiles are shifted in the same way, without aggregating the cells again. The notebook checks this on every plate.
# 
# **Input** (in `../data/`): the position-corrected single-cell profiles of step 4b and the bulk profiles of step 5.
# 
# **Output** (in `../data/`, which git ignores): the three kinds of files above and a copy of `slot_a_adjustment.npz` (with its jackknife estimates, which the hit calling uses), and in this folder `results/slot_a_adjustment_summary.csv`.

# In[1]:


import pathlib
import shutil
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.append("../../../utils")
import slot_adjustment_utils as sau

data_dir = pathlib.Path("../data")
single_cell_dir = data_dir / "single_cell_profiles"
bulk_dir = data_dir / "bulk_profiles"
results_dir = pathlib.Path("results")
results_dir.mkdir(exist_ok=True)

adjustment = sau.load_slot_adjustment("slot_a_adjustment.npz")
print(f"Shrinkage factor {adjustment['factor']:.3f} on {len(adjustment['features'])} features (p = {adjustment['p_value']:.3f})")

# the estimate goes next to the profiles, so that the hit calling can read it with them
shutil.copy("slot_a_adjustment.npz", data_dir / "slot_a_adjustment.npz")


# ## Bulk profiles
# 
# The adjusted bulk profiles are made here and written after the checks below have passed.

# In[2]:


bulk_path = bulk_dir / "bulk_position_corrected_feature_selected.parquet"
bulk = pd.read_parquet(bulk_path)
bulk_adjusted = sau.apply_slot_a_adjustment(bulk, adjustment)
print(f"{len(bulk_adjusted)} wells, {int(bulk_adjusted['Metadata_slot_a_adjusted'].sum())} shifted (the compound wells of slot A)")


# ## Check every plate before anything is written
# 
# For every plate we adjust the single cells with all features and with the features of step 4b, and check that
# 
# - only the cells of the compounds in slot A wells changed,
# - the median of the adjusted cells of a well equals the adjusted bulk profile of the well.
# 
# The second check is the same median computed in two ways, so the two agree to rounding (about 1e-15 for features in standard deviations of a single cell).
# The notebook stops before it writes any file if a difference is not finite or is larger than `TOLERANCE`, which means that the bulk profiles and the single-cell profiles do not come from the same preprocessing run.

# In[3]:


# 1e-6 is far above the rounding of the two computations and far below any difference that preprocessing runs would give
TOLERANCE = 1e-6


def adjust_plate(path):
    """The adjusted cells of a plate with all features and with the features of step 4b."""
    plate = path.name.replace("_sc_position_corrected.parquet", "")
    cells = pd.read_parquet(path)
    adjusted = sau.apply_slot_a_adjustment(cells, adjustment)
    selected = [
        name
        for name in pq.read_schema(single_cell_dir / f"{plate}_sc_position_corrected_feature_selected.parquet").names
        if name in adjusted.columns
    ]
    adjusted_selected = adjusted[selected + ["Metadata_slot", "Metadata_slot_a_adjusted"]]
    return plate, cells, adjusted, adjusted_selected, selected


def check_plate(plate, cells, adjusted, adjusted_selected, selected):
    """Check one plate against the rules above and return its row of the summary."""
    features = [name for name in adjustment["features"] if name in cells.columns]
    before, after = cells[features].to_numpy(), adjusted[features].to_numpy()
    # a missing value stays missing, so two missing values count as equal
    changed = ~((before == after) | (np.isnan(before) & np.isnan(after))).all(axis=1)
    shifted = adjusted["Metadata_slot_a_adjusted"].to_numpy()
    # only the shifted cells changed (a cell whose shift is exactly zero cannot be told from an unchanged one)
    if changed[~shifted].any():
        raise ValueError(f"{plate}: a cell outside slot A changed")

    # the median of the adjusted cells of a well is the adjusted bulk profile of the well
    selected_features = [name for name in selected if not name.startswith("Metadata_")]
    medians = adjusted_selected.groupby("Metadata_Well")[selected_features].median()
    expected = bulk_adjusted[bulk_adjusted["Metadata_Plate"] == plate].set_index("Metadata_Well")
    if expected.empty:
        raise ValueError(f"{plate}: the bulk profiles have no wells of this plate")
    shared_features = [name for name in selected_features if name in expected.columns]
    difference = np.abs(medians.loc[expected.index, shared_features].to_numpy() - expected[shared_features].to_numpy()).max()
    if not np.isfinite(difference) or difference > TOLERANCE:
        raise ValueError(
            f"{plate}: the median of the adjusted cells differs from the adjusted bulk profile by {difference} (tolerance {TOLERANCE})"
        )
    return {
        "plate": plate,
        "cells": len(cells),
        "cells_shifted": int(shifted.sum()),
        "wells_shifted": int(adjusted.loc[shifted, "Metadata_Well"].nunique()),
        "features_all": len(features),
        "features_selected": len(selected_features),
        "mean_abs_shift": float(np.nanmean(np.abs(before[shifted] - after[shifted]))),
        "max_difference_to_bulk": float(difference),
    }


paths = sorted(single_cell_dir.glob("*_sc_position_corrected.parquet"))
records = []
for path in paths:
    records.append(check_plate(*adjust_plate(path)))
    print(records[-1], flush=True)
print("All plates passed the checks")


# ## Write the profiles
# 
# Nothing has been written until this point.

# In[4]:


bulk_adjusted.to_parquet(bulk_dir / "bulk_slot_a_adjusted_feature_selected.parquet", index=False)
for path in paths:
    plate, cells, adjusted, adjusted_selected, selected = adjust_plate(path)
    adjusted.to_parquet(single_cell_dir / f"{plate}_sc_slot_a_adjusted.parquet", index=False)
    adjusted_selected.to_parquet(single_cell_dir / f"{plate}_sc_slot_a_adjusted_feature_selected.parquet", index=False)
    print("Wrote the adjusted profiles of", plate, flush=True)
summary = pd.DataFrame(records)
summary.to_csv(results_dir / "slot_a_adjustment_summary.csv", index=False)
summary.round(6)

