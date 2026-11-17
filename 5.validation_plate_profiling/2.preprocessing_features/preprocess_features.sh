#!/bin/bash

# create or update the environment from uv.lock (needs uv: https://docs.astral.sh/uv/)
# --locked stops with an error instead of changing uv.lock when pyproject.toml and the lockfile disagree
uv sync --locked

# convert notebooks to scripts
uv run --locked jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# convert the CellProfiler SQLite outputs to parquet files (all plates)
uv run --locked python nbconverted/0.convert_cytotable.py

# run single-cell QC for each plate as a script (no image rendering, so it is fast)
# (thresholds per plate are set in sc_qc_thresholds.json)
for plate_dir in ../1.cellprofiler_processing/cp_output/CARD*/; do
    plate_id=$(basename "$plate_dir")
    echo "Running QC for plate ${plate_id}"
    PLATE_ID="${plate_id}" uv run --locked python nbconverted/1.sc_quality_control.py
done
echo "QC completed for all plates"

# process single-cell and bulk profiles for every QC-labeled plate
# (profiles are saved locally and a copy is saved to bandicoot)
uv run --locked python nbconverted/2.single_cell_processing.py
uv run --locked python nbconverted/3.bulk_processing.py

echo "Preprocessing complete."
