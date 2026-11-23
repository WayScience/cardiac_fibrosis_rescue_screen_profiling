#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles
# 
# For each plate, we aggregate the QC-passing single cells to the well level (median per well), annotate the wells with the plate map, normalize to the negative controls, and perform feature selection.
# We then pool the normalized wells of **all validation plates** and perform one feature selection on the pooled profiles, as in the screen ([`4.image_based_profiling/4.bulk_processing.ipynb`](../../4.image_based_profiling/4.bulk_processing.ipynb)).
# 
# NOTE: We normalize the bulk profiles to the negative controls (failing cells treated with DMSO) as the "standard".
# 
# **Input** (per plate, from step 1): `<plate>_qc_labeled.parquet`
# 
# **Output** (in `data/bulk_profiles/`):
# 
# - `<plate>_bulk_annotated.parquet`: well-level profiles of one plate, annotated with the plate map
# - `<plate>_bulk_normalized.parquet`: annotated profiles normalized to the plate's negative controls (`mad_robustize`)
# - `<plate>_bulk_feature_selected.parquet`: normalized profiles after feature selection on that plate
# - `bulk_feature_selected.parquet`: pooled normalized profiles of all plates after one feature selection

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
from pycytominer import aggregate, annotate, feature_select, normalize
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

# Local directory to save bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles/")
output_dir.mkdir(parents=True, exist_ok=True)
output_pooled_feature_select_file = output_dir / "bulk_feature_selected.parquet"

# Directory on bandicoot to save a copy of the bulk profiles
bandicoot_output_dir = pathlib.Path(
    "~/mnt/bandicoot/CFReT_screening_data/validation_profiles/bulk_profiles/"
).expanduser()

# Directory with one platemap per plate (named by plate barcode)
platemap_dir = pathlib.Path("../metadata/platemaps/").resolve(strict=True)

# Negative controls used as the reference population for normalization
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"

# settings for feature selection, used for each plate and for the pooled profiles
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

# suffixes of the files saved for each plate
output_suffixes = ["bulk_annotated", "bulk_normalized", "bulk_feature_selected"]


# ## Process each plate with pycytominer
# 
# For each plate, we drop the QC-failed cells, aggregate to the well level, annotate, normalize to the negative controls, and feature select.

# In[3]:


for plate_id in plate_ids:
    # Load the platemap for this plate (each plate has its own compounds)
    platemap_df = pd.read_csv(
        (platemap_dir / f"{plate_id}_platemap.csv").resolve(strict=True)
    )

    print("Performing preprocessing on", plate_id)

    output_annotated_file = output_dir / f"{plate_id}_bulk_annotated.parquet"
    output_normalized_file = output_dir / f"{plate_id}_bulk_normalized.parquet"
    output_feature_select_file = output_dir / f"{plate_id}_bulk_feature_selected.parquet"

    profile_df = pd.read_parquet(qc_labeled_dir / f"{plate_id}_qc_labeled.parquet")

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Aggregate single-cell data to the well level using the median
    print("Performing aggregation for", plate_id, "...")
    aggregated_df = aggregate(
        population_df=profile_df,
        operation="median",
        strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
    )

    # Step 2: Annotation
    # (annotate also adds the Metadata_ prefix to the Image_Metadata_ columns)
    print("Performing annotation for", plate_id, "...")
    annotate(
        profiles=aggregated_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_file=output_annotated_file,
        output_type="parquet",
    )

    # Step 3: Normalization (mad robustize) to the negative controls
    print("Performing normalization for", plate_id, "...")
    normalize(
        profiles=output_annotated_file,
        method="mad_robustize",
        samples=neg_control_query,
        output_file=output_normalized_file,
        output_type="parquet",
    )

    # Step 4: Feature selection on this plate
    print("Performing feature selection for", plate_id, "...")
    feature_select(
        profiles=output_normalized_file,
        output_file=output_feature_select_file,
        output_type="parquet",
        **feature_select_kwargs,
    )

print(f"All plates have been processed and saved to: {output_dir.resolve()}")


# ## Pool the wells of all plates and feature select
# 
# We concatenate the normalized well-level profiles of all plates and perform one feature selection on the pooled profiles, so that all wells share one set of features.

# In[4]:


pooled_df = pd.concat(
    [
        pd.read_parquet(output_dir / f"{plate_id}_bulk_normalized.parquet")
        for plate_id in plate_ids
    ],
    ignore_index=True,
)
print("Pooled well-level profiles:", pooled_df.shape)
print("Plates:", pooled_df["Metadata_Plate"].nunique())
print("Negative control wells:", len(pooled_df.query(neg_control_query)))

print("Feature selecting the pooled profiles...")
feature_select(
    profiles=pooled_df,
    output_file=output_pooled_feature_select_file,
    output_type="parquet",
    **feature_select_kwargs,
)

print(f"Saved pooled feature-selected profiles to {output_pooled_feature_select_file}")


# ## Save a copy of the profiles to bandicoot

# In[5]:


# copy every per-plate and pooled bulk profile saved by this notebook
output_files = [
    output_dir / f"{plate_id}_{suffix}.parquet"
    for plate_id in plate_ids
    for suffix in output_suffixes
] + [output_pooled_feature_select_file]

for output_file in output_files:
    output(
        df=pd.read_parquet(output_file),
        output_filename=bandicoot_output_dir / output_file.name,
        output_type="parquet",
    )

print(f"Saved a copy of {len(output_files)} bulk profiles to: {bandicoot_output_dir}")

