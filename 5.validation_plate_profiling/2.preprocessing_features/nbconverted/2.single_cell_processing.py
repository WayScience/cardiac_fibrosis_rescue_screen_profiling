#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the plates for all samples as we only have three wells associated with the healthy controls, which is insufficient for normalization.

# ## Import libraries

# In[ ]:


import pathlib
import sys
from functools import partial

import pandas as pd

from pycytominer import annotate, feature_select, normalize

sys.path.append("../../utils")
import io_utils


# ## Set paths and variables

# In[ ]:


# Where to save the processed profiles: "local", "bandicoot", or "both"
save_to = "both"

# Directory with QC-labeled profiles
qc_labeled_dir = pathlib.Path("./data/qc_labeled_profiles/").resolve(strict=True)

# Find all plates that have a QC-labeled profile (one per validation plate)
plate_ids = sorted(
    path.name.replace("_qc_labeled.parquet", "")
    for path in qc_labeled_dir.glob("*_qc_labeled.parquet")
)
print(f"There are {len(plate_ids)} plate(s) to process:")
for plate_id in plate_ids:
    print(plate_id)

# Local directory to save single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles/")

# Directory on bandicoot to save single-cell profiles
bandicoot_output_dir = pathlib.Path(
    "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/single_cell_profiles/"
).expanduser()

# Save a profile to the location(s) set by save_to
save_profile = partial(
    io_utils.save_parquet,
    save_to=save_to,
    local_dir=output_dir,
    bandicoot_dir=bandicoot_output_dir,
)

# Directory to read saved profiles back from in the checks below
saved_dir = bandicoot_output_dir if save_to == "bandicoot" else output_dir

# Path to the platemap (all validation plates share the same layout)
platemap_path = pathlib.Path(
    "../metadata/platemaps/platemap_validation.csv"
).resolve(strict=True)

# operations to perform for feature selection
feature_select_ops = [
    "blocklist",
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
]


# ## Process data with pycytominer

# In[ ]:


# Load the platemap once since all plates share the same layout
platemap_df = pd.read_csv(platemap_path)

for plate_id in plate_ids:
    print("Performing preprocessing on", plate_id)

    profile_df = pd.read_parquet(qc_labeled_dir / f"{plate_id}_qc_labeled.parquet")

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Annotation
    print("Performing annotation for", plate_id, "...")
    annotated_df = annotate(
        profiles=profile_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
    )

    # Rename site column to have the Metadata prefix
    annotated_df = annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"})
    save_profile(annotated_df, f"{plate_id}_sc_annotated.parquet")

    # Step 2: Normalization
    print("Performing normalization for", plate_id, "...")
    normalized_df = normalize(
        profiles=annotated_df,
        method="standardize",
        samples="all",
    )
    save_profile(normalized_df, f"{plate_id}_sc_normalized.parquet")

    # Step 3: Feature selection
    print("Performing feature selection for", plate_id, "...")
    feature_selected_df = feature_select(
        profiles=normalized_df,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
    )
    save_profile(feature_selected_df, f"{plate_id}_sc_feature_selected.parquet")

    print(f"Annotation, normalization, and feature selection complete for {plate_id}")

print(f"All plates have been processed and saved to: {save_to}")


# In[ ]:


# Check the feature-selected output file for each plate
for plate_id in plate_ids:
    test_df = pd.read_parquet(saved_dir / f"{plate_id}_sc_feature_selected.parquet")
    print(plate_id, test_df.shape, "| Plate:", test_df.Metadata_Plate.unique())

print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[ ]:


# Check the annotated output file for each plate
for plate_id in plate_ids:
    test_df = pd.read_parquet(saved_dir / f"{plate_id}_sc_annotated.parquet")
    print(plate_id, test_df.shape, "| Plate:", test_df.Metadata_Plate.unique())

print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

