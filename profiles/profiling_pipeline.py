"""
Perform the profiling pipeline (defined in profile.py) given all plates.
"""

import os
import pathlib
import subprocess
import pandas as pd

# Load constants
project = "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad"
batch = "2016_04_01_a549_48hr_batch1"
profile_dir = pathlib.PurePath(
    "/home/ubuntu/bucket/projects/", project, "workspace/backend"
)
barcode_platemap_dir = pathlib.PurePath("../metadata/platemaps")

# Load barcode platemap information
barcode_platemap_file = pathlib.PurePath(barcode_platemap_dir, "barcode_platemap.csv")
barcode_platemap_df = pd.read_csv(barcode_platemap_file)

# Load platemap information
platemap_dir = pathlib.PurePath(barcode_platemap_dir, "platemap")

# Load plate information
plate_dir = pathlib.PurePath(profile_dir, batch)
plates = [x for x in os.listdir(plate_dir) if x.startswith("SQ")]

# Load and check MOA information
moa_file = pathlib.PurePath(
    "../metadata/moa/repurposing_info_external_moa_map_resolved.tsv"
)
moa_df = pd.read_csv(moa_file, sep="\t")
assert isinstance(
    moa_df, pd.DataFrame
), "Error, MOA file does not exist. Is the path updated?"

# Process every plate
for plate in plates:
    print(f"Now processing... Plate: {plate}")
    output_dir = pathlib.PurePath("backend", plate)
    cell_count_dir = pathlib.PurePath("analysis", plate)

    platemap_id = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]

    platemap_file = pathlib.PurePath(platemap_dir, f"{platemap_id}.txt")
    sql_base = pathlib.PurePath(profile_dir, batch, plate, f"{plate}.sqlite")
    sql_file = f"sqlite:////{sql_base}"

    cmd = [
        "python",
        "profile.py",
        "--sql_file",
        sql_file,
        "--batch",
        batch,
        "--plate_name",
        plate,
        "--platemap_file",
        platemap_file,
        "--barcode_platemap_file",
        barcode_platemap_file,
        "--moa_file",
        moa_file,
        "--output_dir",
        output_dir,
        "--cell_count_dir",
        cell_count_dir,
    ]
    subprocess.call(cmd)
