#!/bin/bash

# create or update the environment from uv.lock (needs uv: https://docs.astral.sh/uv/)
# --locked stops with an error instead of changing uv.lock when pyproject.toml and the lockfile disagree
uv sync --locked

# convert notebooks to scripts
uv run --locked jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run the notebooks in order
uv run --locked python nbconverted/0.convert_cytotable.py
uv run --locked python nbconverted/1.sc_quality_control.py
uv run --locked python nbconverted/2.single_cell_processing.py
uv run --locked python nbconverted/3.bulk_processing.py

echo "Preprocessing complete."
