#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the plates for all samples as we only have three wells associated with the healthy controls, which is insufficient for normalization.

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
from pycytominer import annotate, feature_select, normalize
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[2]:


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
output_dir.mkdir(parents=True, exist_ok=True)

# Directory on bandicoot to save a copy of the single-cell profiles
bandicoot_output_dir = pathlib.Path(
    "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/single_cell_profiles/"
).expanduser()

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

# suffixes of the files saved for each plate
output_suffixes = ["sc_annotated", "sc_normalized", "sc_feature_selected"]


# ## Process data with pycytominer

# In[3]:


# Load the platemap once since all plates share the same layout
platemap_df = pd.read_csv(platemap_path)

for plate_id in plate_ids:
    print("Performing preprocessing on", plate_id)

    output_annotated_file = output_dir / f"{plate_id}_sc_annotated.parquet"
    output_normalized_file = output_dir / f"{plate_id}_sc_normalized.parquet"
    output_feature_select_file = output_dir / f"{plate_id}_sc_feature_selected.parquet"

    profile_df = pd.read_parquet(qc_labeled_dir / f"{plate_id}_qc_labeled.parquet")

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Annotation
    # (annotate also adds the Metadata_ prefix to the Image_Metadata_ columns)
    print("Performing annotation for", plate_id, "...")
    annotate(
        profiles=profile_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_file=output_annotated_file,
        output_type="parquet",
    )

    # Step 2: Normalization
    print("Performing normalization for", plate_id, "...")
    normalize(
        profiles=output_annotated_file,
        method="standardize",
        samples="all",
        output_file=output_normalized_file,
        output_type="parquet",
    )

    # Step 3: Feature selection
    print("Performing feature selection for", plate_id, "...")
    feature_select(
        profiles=output_normalized_file,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
        output_file=output_feature_select_file,
        output_type="parquet",
    )

    print(f"Annotation, normalization, and feature selection complete for {plate_id}")

print(f"All plates have been processed and saved to: {output_dir.resolve()}")


# ## Save a copy of the profiles to bandicoot

# In[4]:


for plate_id in plate_ids:
    for suffix in output_suffixes:
        file_name = f"{plate_id}_{suffix}.parquet"
        output(
            df=pd.read_parquet(output_dir / file_name),
            output_filename=bandicoot_output_dir / file_name,
            output_type="parquet",
        )

print(f"Saved a copy of all profiles to: {bandicoot_output_dir}")

