"""
This module provides functions for loading, writing, and modifying files in this data
repository. It centralizes input/output operations, ensuring efficient and consistent
handling of data files throughout the project.
"""

import pathlib
from collections import defaultdict

import pandas as pd


def load_barcodes(barcode_path: str | pathlib.Path) -> dict:
    """Load barcode data from a given file and organize it into batches.

    This function reads barcode data, groups it by the 'platemap_file' column, and
    organizes the grouped data into batches. Each batch is assigned a unique batch ID
    (e.g., 'batch_1', 'batch_2', etc.).

    Parameters
    ----------
    barcode_path : str | pathlib.Path
        Path to the file containing barcode data. This path can be a string or a pathlib.Path
        object.

    Returns
    -------
    dict
        A dictionary where keys are batch IDs (e.g., "batch_1", "batch_2") and values are
        dictionaries.  Each inner dictionary maps a 'platemap_file' name to the
        corresponding 'plate_barcode' data.
    """

    # type checking
    if not isinstance(barcode_path, (str | pathlib.Path)):
        raise TypeError("'barcode_path' must a str or pathlib.Path")
    if isinstance(barcode_path, str):
        barcode_path = pathlib.Path(barcode_path)

    # load barcode as csv
    barcodes = pd.read_csv(barcode_path)

    # Initialize an empty dictionary to store barcode contents
    barcode_contents = defaultdict(lambda: None)

    # Group barcodes by 'platemap_file' and iterate through each group
    for idx, (platemap_name, df) in enumerate(
        barcodes.groupby("platemap_file"), start=1
    ):
        # Generate a batch ID and plate ID
        batch_id = f"batch_{idx}"

        # Collect plate barcodes for the current platemap
        batch_content = {platemap_name: df["plate_barcode"].values.tolist()}

        # Store the batch content in the main dictionary
        barcode_contents[batch_id] = batch_content

    return barcode_contents


def save_parquet(
    df: pd.DataFrame,
    file_name: str,
    save_to: str,
    local_dir: str | pathlib.Path,
    bandicoot_dir: str | pathlib.Path,
) -> list[pathlib.Path]:
    """Save a DataFrame as a parquet file locally, to bandicoot, or to both.

    Parameters
    ----------
    df : pd.DataFrame
        The data to save.
    file_name : str
        Name of the parquet file (e.g., "PLATE_sc_annotated.parquet").
    save_to : str
        Where to save the file: "local", "bandicoot", or "both".
    local_dir : str | pathlib.Path
        Local directory to save into (created if it does not exist).
    bandicoot_dir : str | pathlib.Path
        Directory on the mounted bandicoot drive to save into. The parent directory
        must already exist, so a missing mount raises an error instead of silently
        writing to the local disk.

    Returns
    -------
    list[pathlib.Path]
        The paths the file was saved to.
    """

    if save_to not in {"local", "bandicoot", "both"}:
        raise ValueError("'save_to' must be one of 'local', 'bandicoot', or 'both'")

    output_dirs = []
    if save_to in {"local", "both"}:
        local_dir = pathlib.Path(local_dir)
        local_dir.mkdir(parents=True, exist_ok=True)
        output_dirs.append(local_dir)
    if save_to in {"bandicoot", "both"}:
        bandicoot_dir = pathlib.Path(bandicoot_dir).expanduser()
        if not bandicoot_dir.parent.is_dir():
            raise FileNotFoundError(
                f"Bandicoot directory not found: {bandicoot_dir.parent}. "
                "Is bandicoot mounted?"
            )
        bandicoot_dir.mkdir(exist_ok=True)
        output_dirs.append(bandicoot_dir)

    saved_paths = []
    for output_dir in output_dirs:
        output_path = output_dir / file_name
        df.to_parquet(output_path, index=False)
        saved_paths.append(output_path)

    return saved_paths
