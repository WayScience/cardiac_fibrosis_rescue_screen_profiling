#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles profiles
# 
# NOTE: We are normalizing the bulk-profile plates to the negative controls as the "standard".

# ## Import libraries

# In[1]:


import pathlib
import sys
from functools import partial

import pandas as pd

from pycytominer import aggregate, annotate, feature_select, normalize

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

# Local directory to save bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles/")

# Directory on bandicoot to save bulk profiles
bandicoot_output_dir = pathlib.Path(
    "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/bulk_profiles/"
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
# 

# In[3]:


# Load the platemap once since all plates share the same layout
platemap_df = pd.read_csv(platemap_path)

# Negative controls used as the reference population
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"

for plate_id in plate_ids:
    print("Performing preprocessing on", plate_id)

    # loading profiles
    profile_df = pd.read_parquet(qc_labeled_dir / f"{plate_id}_qc_labeled.parquet")

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Aggregate single-cell data to the well-level using the median
    print("Performing aggregation for", plate_id, "...")
    aggregated_df = aggregate(
        population_df=profile_df,
        operation="median",
        strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
    )

    # Step 2: Annotation
    print("Performing annotation for", plate_id, "...")
    annotated_df = annotate(
        profiles=aggregated_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
    )

    # Rename site column to have the Metadata prefix
    annotated_df = annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"})
    save_profile(annotated_df, f"{plate_id}_bulk_annotated.parquet")

    # Step 3: Normalization (mad robustize)
    # Normalize using the negative controls as the reference population
    print("Performing normalization for", plate_id, "...")
    normalized_df = normalize(
        profiles=annotated_df,
        method="mad_robustize",
        samples=neg_control_query,
    )
    save_profile(normalized_df, f"{plate_id}_bulk_normalized.parquet")

    # Step 4: Feature selection
    print("Performing feature selection for", plate_id, "...")
    feature_select_df = feature_select(
        profiles=normalized_df,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
        corr_threshold=0.95,
        freq_cut=0.05,
    )
    save_profile(feature_select_df, f"{plate_id}_bulk_feature_selected.parquet")

    # Step 4b: Remove features with too little variation inside the exact control
    # population used to fit spherization.
    print(
        f"Feature selecting {plate_id} with variance threshold "
        "within negative controls only..."
    )
    zero_negcon_var_fs_df = feature_select(
        profiles=feature_select_df,
        operation="variance_threshold",
        freq_cut=0.05,
        unique_cut=0.01,
        samples=neg_control_query,
    )

    # Step 5: Spherize/whiten all profiles using the negative controls as the
    # reference population.
    print(f"Sphering {plate_id} using negative controls...")
    spherized_df = normalize(
        profiles=zero_negcon_var_fs_df,
        method="spherize",
        samples=neg_control_query,
        spherize_center=True,
        spherize_method="ZCA-cor",
        spherize_epsilon=1e-6,
    )
    save_profile(spherized_df, f"{plate_id}_bulk_spherized.parquet")

    print(
        "Saved annotated, normalized, feature-selected, and spherized profiles for",
        plate_id,
    )

print(f"All plates have been processed and saved to: {save_to}")


# In[4]:


# Check the spherized output file for each plate
for plate_id in plate_ids:
    test_df = pd.read_parquet(saved_dir / f"{plate_id}_bulk_spherized.parquet")
    print(plate_id, test_df.shape, "| Plate:", test_df.Metadata_Plate.unique())

print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[5]:


# Check the feature-selected output file for each plate
for plate_id in plate_ids:
    test_df = pd.read_parquet(saved_dir / f"{plate_id}_bulk_feature_selected.parquet")
    print(plate_id, test_df.shape, "| Plate:", test_df.Metadata_Plate.unique())

print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

