"""
Perform the image-based profiling pipeline to process data
"""

import os
import sys
import pathlib
import pandas as pd
from pycytominer import aggregate, annotate, normalize, feature_select, cyto_utils
from pycytominer.cyto_utils.cells import SingleCells

from profile_utils import get_args

sys.path.append("../utils")
from dose import recode_dose


# Load Command Line Arguments
args = get_args()

sql_file = args.sql_file
batch = args.batch
plate_name = args.plate_name
platemap_file = args.platemap_file
barcode_platemap_file = args.barcode_platemap_file
moa_file = args.moa_file
cell_count_dir = args.cell_count_dir
cell_id = args.cell_id
output_dir = args.output_dir
plate_col = args.plate_col  # Default is "Image_Metadata_Plate"
well_col = args.well_col  # Default is "Image_Metadata_Well"

# Initialize profile processing
os.makedirs(output_dir, exist_ok=True)
os.makedirs(cell_count_dir, exist_ok=True)

aggregate_method = "median"
norm_method = "mad_robustize"
compression = {"method": "gzip", "mtime": 1}
float_format = "%.5g"
strata = [plate_col, well_col]
feature_select_ops = [
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
]
primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]

# Define external metadata to add to annotation
moa_df = pd.read_csv(moa_file, sep="\t")
barcode_platemap_df = pd.read_csv(barcode_platemap_file).query(
    "Assay_Plate_Barcode == @plate_name"
)

# Aggregate profiles
out_file = pathlib.PurePath(output_dir, f"{plate_name}.csv.gz")
sc = SingleCells(
    file_or_conn=sql_file, strata=strata, aggregation_operation=aggregate_method
)
sc.aggregate_profiles(
    output_file=out_file, float_format=float_format, compression_options=compression
)

# Count cells
count_file = pathlib.PurePath(cell_count_dir, f"{plate_name}_cell_count.csv")
cell_count_df = sc.count_cells()
cell_count_df.to_csv(count_file, sep=",", index=False)

del sc

# Annotate profiles - Level 3 Data
anno_file = pathlib.PurePath(output_dir, f"{plate_name}_augmented.csv.gz")
anno_df = annotate(
    profiles=out_file,
    platemap=platemap_file,
    join_on=["Metadata_well_position", well_col],
    format_broad_cmap=True,
    external_metadata=moa_df,
    external_join_left=["Metadata_broad_sample"],
    external_join_right=["Metadata_broad_sample"],
    cmap_args={"cell_id": cell_id, "perturbation_mode": "chemical"},
)

# Rename columns
anno_df = anno_df.rename(
    {"Image_Metadata_Plate": "Metadata_Plate", "Image_Metadata_Well": "Metadata_Well"},
    axis="columns",
)

# Add barcode platemap info
try:
    anno_df = anno_df.assign(
        Metadata_Assay_Plate_Barcode=barcode_platemap_df.Assay_Plate_Barcode.values[0],
        Metadata_Plate_Map_Name=barcode_platemap_df.Plate_Map_Name.values[0],
        Metadata_Batch_Number=barcode_platemap_df.Batch_Number.values[0],
        Metadata_Batch_Date=barcode_platemap_df.Batch_Date.values[0],
    )
except AttributeError:
    anno_df = anno_df.assign(
        Metadata_Assay_Plate_Barcode=barcode_platemap_df.Assay_Plate_Barcode.values[0],
        Metadata_Plate_Map_Name=barcode_platemap_df.Plate_Map_Name.values[0],
    )

# Add dose recoding information
anno_df = anno_df.assign(
    Metadata_dose_recode=(
        anno_df.Metadata_mmoles_per_liter.apply(
            lambda x: recode_dose(x, primary_dose_mapping, return_level=True)
        )
    )
)

# Reoroder columns
metadata_cols = cyto_utils.infer_cp_features(anno_df, metadata=True)
cp_cols = cyto_utils.infer_cp_features(anno_df)
reindex_cols = metadata_cols + cp_cols
anno_df = anno_df.reindex(reindex_cols, axis="columns")

# Output annotated file
cyto_utils.output(
    df=anno_df,
    output_filename=anno_file,
    float_format=float_format,
    compression_options=compression,
)

# Normalize Profiles (DMSO Control) - Level 4A Data
norm_dmso_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized_dmso.csv.gz")
normalize(
    profiles=anno_df,
    samples="Metadata_broad_sample == 'DMSO'",
    method=norm_method,
    output_file=norm_dmso_file,
    float_format=float_format,
    compression_options=compression,
)

# Normalize Profiles (Whole Plate) - Level 4A Data
norm_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized.csv.gz")
normalize(
    profiles=anno_df,
    samples="all",
    method=norm_method,
    output_file=norm_file,
    float_format=float_format,
    compression_options=compression,
)

# Feature Selection (DMSO Control) - Level 4B Data
feat_dmso_file = pathlib.PurePath(
    output_dir, f"{plate_name}_normalized_feature_select_dmso.csv.gz"
)
feature_select(
    profiles=norm_dmso_file,
    features="infer",
    operation=feature_select_ops,
    output_file=feat_dmso_file,
    float_format=float_format,
    compression_options=compression,
)

# Feature Selection (Whole Plate) - Level 4B Data
feat_file = pathlib.PurePath(
    output_dir, f"{plate_name}_normalized_feature_select.csv.gz"
)
feature_select(
    profiles=norm_file,
    features="infer",
    operation=feature_select_ops,
    output_file=feat_file,
    float_format=float_format,
    compression_options=compression,
)
