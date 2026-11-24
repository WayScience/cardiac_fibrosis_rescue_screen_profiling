#!/usr/bin/env python
# coding: utf-8

# # Apply the plate-position correction to single cells
# 
# This step applies the position correction that `4a.position_correction/1.fit_position_correction.ipynb` saved to the normalized single-cell profiles from step 2, then performs feature selection on the corrected profiles with the same settings as step 2.
# It is the step 3b of the screen (`../../4.image_based_profiling/3b.apply_position_correction.ipynb`) for the validation plates.
# 
# For every cell, we subtract `amplitude x tilt(row, column)` from the features, where the amplitude is that of the plate (the fit has one amplitude per plate, shrunk toward the other plates like the platemap amplitudes of the screen).
# We shift all cells in a well by the same amount.
# We correct controls like all other wells.
# 
# **Input** (per plate): `<plate>_sc_normalized.parquet` from step 2, and the fit `4a.position_correction/position_correction_fit.npz`.
# The notebook reads the profiles from `data/single_cell_profiles/` when that folder exists, and from the folder on bandicoot otherwise; set `VALIDATION_PROFILE_DIR` to read another folder.
# 
# **Output** (per plate, in `data/single_cell_profiles/`, with a copy on bandicoot)
# 
# - `<plate>_sc_position_corrected.parquet`: normalized profiles after correction
# - `<plate>_sc_position_corrected_feature_selected.parquet`: corrected profiles after feature selection
# 
# The profiles of step 2 are not changed.

# ## Import libraries

# In[1]:


import os
import pathlib
import pprint
import sys

import pandas as pd
from pycytominer import feature_select
from pycytominer.cyto_utils import output

sys.path.append("../../utils")
import position_correction_utils as pcu
import validation_position_utils as vpu


# ## Set paths and variables

# In[2]:


# Local directory to save single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles/")
output_dir.mkdir(parents=True, exist_ok=True)

# Directory on bandicoot with the profiles of step 2, and where we save a copy of the new profiles
# (set BANDICOOT_OUTPUT_DIR to save the copy somewhere else)
bandicoot_output_dir = pathlib.Path(
    os.environ.get(
        "BANDICOOT_OUTPUT_DIR",
        "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/single_cell_profiles/",
    )
).expanduser()

# Directory with the normalized profiles of step 2
if os.environ.get("VALIDATION_PROFILE_DIR"):
    input_dir = pathlib.Path(os.environ["VALIDATION_PROFILE_DIR"]).expanduser()
elif any(output_dir.glob("*_sc_normalized.parquet")):
    input_dir = output_dir
else:
    input_dir = bandicoot_output_dir
input_dir = input_dir.resolve(strict=True)
print("Reading the normalized profiles from", input_dir)

# position correction fit written by 4a.position_correction/1.fit_position_correction.ipynb
fit_path = pathlib.Path("./4a.position_correction/position_correction_fit.npz").resolve(
    strict=True
)

# operations to perform for feature selection (same as step 2)
feature_select_ops = [
    "blocklist",
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
]

# suffixes of the files saved for each plate
output_suffixes = ["sc_position_corrected", "sc_position_corrected_feature_selected"]


# ## Load the fit
# 
# The fit has the tilt of every feature and the amplitude of every plate.

# In[3]:


fit = pcu.load_fit(fit_path)
if not set(range(1, 6)) <= set(fit["amplitudes"]):
    raise ValueError(
        "No amplitude for every plate in the position correction fit. "
        "Rerun 4a.position_correction/1.fit_position_correction.ipynb."
    )
print("Amplitude of each plate:", fit["amplitudes"])
print("Features in fit:", len(fit["features"]))


# ## Set list of plates to process

# In[4]:


plate_info_dictionary = {
    row.plate: {"normalized_path": row.path, "plate_number": row.platemap}
    for row in vpu.discover_validation_profiles(input_dir).itertuples()
}

# View dictionary
print("Number of plates to process:", len(plate_info_dictionary))
pprint.pprint(plate_info_dictionary, indent=4)


# ## Correct the plates and select features

# In[5]:


for plate, info in plate_info_dictionary.items():
    output_corrected_file = output_dir / f"{plate}_sc_position_corrected.parquet"
    output_feature_select_file = (
        output_dir / f"{plate}_sc_position_corrected_feature_selected.parquet"
    )

    print("Applying position correction to", plate)
    normalized_df = pd.read_parquet(info["normalized_path"])

    # Step 1: subtract the tilt of every well from its cells
    corrected_df = pcu.correct_cells(normalized_df, fit, info["plate_number"])
    n_corrected = len([f for f in fit["features"] if f in normalized_df.columns])
    print(f"Corrected {n_corrected} features in {len(corrected_df)} cells")
    corrected_df.to_parquet(output_corrected_file, index=False)

    # Step 2: feature selection on the corrected profiles
    print("Performing feature selection for", plate, "...")
    feature_select(
        profiles=corrected_df,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
        output_file=str(output_feature_select_file),
        output_type="parquet",
    )
    print(f"Position correction and feature selection complete for {plate}")

print(f"All plates have been processed and saved to: {output_dir.resolve()}")


# ## Check the last plate
# 
# The corrected profiles keep the cells and the metadata of the normalized profiles.

# In[6]:


normalized_df = pd.read_parquet(info["normalized_path"])
corrected_df = pd.read_parquet(output_corrected_file)
selected_df = pd.read_parquet(output_feature_select_file)

metadata_cols = [c for c in normalized_df.columns if c.startswith("Metadata_")]
print("Normalized:", normalized_df.shape)
print("Corrected:", corrected_df.shape)
print("Corrected and feature selected:", selected_df.shape)
print(
    "Metadata unchanged:",
    normalized_df[metadata_cols].equals(corrected_df[metadata_cols]),
)

feature_cols = [f for f in fit["features"] if f in normalized_df.columns]
shift = (normalized_df[feature_cols] - corrected_df[feature_cols]).abs()
print(f"Mean absolute shift per feature: {shift.mean().mean():.3f}")
selected_df.head(2)


# ## Save a copy of the profiles to bandicoot

# In[7]:


for plate in plate_info_dictionary:
    for suffix in output_suffixes:
        file_name = f"{plate}_{suffix}.parquet"
        output(
            df=pd.read_parquet(output_dir / file_name),
            output_filename=bandicoot_output_dir / file_name,
            output_type="parquet",
        )
print(f"Saved a copy of all position-corrected profiles to: {bandicoot_output_dir}")

