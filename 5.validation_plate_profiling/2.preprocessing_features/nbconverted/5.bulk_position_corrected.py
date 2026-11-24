#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles from position-corrected cells
# 
# We aggregate the position-corrected single-cell profiles from step 4b to the well level (median per well) and pool the wells of **all validation plates**.
# The single-cell profiles are already standardized (z-scored using all cells on each plate) and corrected for plate-position effects, so we do not normalize the bulk profiles again.
# We then perform one feature selection on the pooled profiles.
# This is the bulk step of the screen (`../../4.image_based_profiling/4.bulk_processing.ipynb`) for the validation plates.
# 
# `3.bulk_processing.ipynb` builds bulk profiles from the uncorrected cells and normalizes them to the failing DMSO controls.
# The profiles of this notebook have a different name, and the profiles of that notebook are not changed.
# 
# **Input** (per plate, from step 4b): `<plate>_sc_position_corrected.parquet`
# The notebook reads the profiles from `data/single_cell_profiles/` when that folder exists, and from the folder on bandicoot otherwise; set `VALIDATION_CORRECTED_PROFILE_DIR` to read another folder.
# 
# **Output** (in `data/bulk_profiles/`)
# 
# - `bulk_position_corrected_feature_selected.parquet`: pooled well-level profiles of all plates after feature selection

# ## Import libraries

# In[1]:


import os
import pathlib
import sys

import pandas as pd
from pycytominer import aggregate, feature_select

sys.path.append("../../utils")
import validation_position_utils as vpu


# ## Set paths and variables

# In[2]:


# Local directory to save bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles/")
output_dir.mkdir(parents=True, exist_ok=True)
output_pooled_feature_select_file = output_dir / "bulk_position_corrected_feature_selected.parquet"

# Directory with the position-corrected single-cell profiles of step 4b
# (set VALIDATION_CORRECTED_PROFILE_DIR to read another folder; VALIDATION_PROFILE_DIR is the
# folder of the normalized profiles that step 4b reads, not of its output)
local_profile_dir = pathlib.Path("./data/single_cell_profiles/")
bandicoot_profile_dir = pathlib.Path(
    "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/single_cell_profiles/"
).expanduser()
if os.environ.get("VALIDATION_CORRECTED_PROFILE_DIR"):
    input_dir = pathlib.Path(os.environ["VALIDATION_CORRECTED_PROFILE_DIR"]).expanduser()
elif any(local_profile_dir.glob("*_sc_position_corrected.parquet")):
    input_dir = local_profile_dir
else:
    input_dir = bandicoot_profile_dir
input_dir = input_dir.resolve(strict=True)
print("Reading the position-corrected profiles from", input_dir)

# settings for feature selection (the same as for the pooled bulk profiles of step 3)
# pycytominer applies the operations in this order, so we first drop features with
# missing values and blocklisted features, and then apply the variance and
# correlation filters
feature_select_kwargs = {
    "operation": [
        "drop_na_columns",
        "blocklist",
        "variance_threshold",
        "correlation_threshold",
    ],
    # drop every feature with a missing value (the pycytominer default allows 5%)
    "na_cutoff": 0,
    "blocklist_file": "./blocklist_features.txt",
    # 0.95 is less strict than the pycytominer default of 0.9, so fewer correlated
    # features are removed and more features remain
    "corr_threshold": 0.95,
    # 0.05 is the pycytominer default for the most common value of a feature
    "freq_cut": 0.05,
}

# columns that identify a well and are constant within it; used as the aggregation
# strata so that the well-level metadata are kept
well_strata = [
    "Metadata_well_row",
    "Metadata_well_column",
    "Metadata_heart_number",
    "Metadata_cell_type",
    "Metadata_heart_failure_type",
    "Metadata_treatment",
    "Metadata_Plate",
    "Metadata_Well",
]

# negative control wells (failing-cell DMSO), counted to check the pooled table
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"


# ## Set list of plates to process
# 
# We process the position-corrected profiles of every validation plate.

# In[3]:


plate_table = vpu.discover_validation_profiles(
    input_dir, suffix="_sc_position_corrected.parquet"
)
if plate_table.empty:
    raise FileNotFoundError(f"No *_sc_position_corrected.parquet files in {input_dir}")

print("Number of plates to process:", len(plate_table))
print(plate_table["plate"].to_string(index=False))


# ## Aggregate to the well level and pool all plates
# 
# For each plate, we compute the median profile of each well.
# QC-failed cells were already removed in step 2 and the profiles are already annotated, so the well-level metadata are kept by using them as the aggregation strata.
# We then concatenate the wells of all plates.

# In[4]:


well_profiles = []
for plate_info in plate_table.itertuples():
    print("Performing aggregation for", plate_info.plate, "...")
    well_profiles.append(
        aggregate(
            population_df=pd.read_parquet(plate_info.path),
            operation="median",
            strata=well_strata,
        )
    )

pooled_df = pd.concat(well_profiles, ignore_index=True)
print("Pooled well-level profiles:", pooled_df.shape)
print("Plates:", pooled_df["Metadata_Plate"].nunique())
print("Negative control wells:", len(pooled_df.query(neg_control_query)))


# ## Feature selection
# 
# We perform one feature selection on the pooled profiles of all plates, so that all wells share one set of features.

# In[5]:


print("Feature selecting the pooled profiles...")
feature_select(
    profiles=pooled_df,
    output_file=output_pooled_feature_select_file,
    output_type="parquet",
    **feature_select_kwargs,
)
print(f"Saved feature-selected profiles to {output_pooled_feature_select_file}")

