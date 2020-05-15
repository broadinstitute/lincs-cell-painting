"""
Perform the image-based profiling pipeline to process data
"""

import os
import pathlib
import pandas as pd
from pycytominer.aggregate import AggregateProfiles
from pycytominer import (
    aggregate,
    annotate,
    normalize,
    feature_select,
    audit,
    cyto_utils,
)

from profile_utils import get_args

# Load Command Line Arguments
args = get_args()

sql_file = args.sql_file
batch = args.batch
plate_name = args.plate_name
platemap_file = args.platemap_file
barcode_platemap_file = args.barcode_platemap_file
moa_file = args.moa_file
cell_count_dir = args.cell_count_dir
output_dir = args.output_dir

# Initialize profile processing
os.makedirs(output_dir, exist_ok=True)
os.makedirs(cell_count_dir, exist_ok=True)
cell_id = "A549"
aggregate_method = "median"
norm_method = "mad_robustize"
compression = "gzip"
float_format = "%.5g"
strata = ["Image_Metadata_Plate", "Image_Metadata_Well"]
feature_select_ops = [
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
    "blacklist",
]

# Define external metadata to add to annotation
moa_df = pd.read_csv(moa_file, sep="\t")
barcode_platemap_df = pd.read_csv(barcode_platemap_file).query(
    "Assay_Plate_Barcode == @plate_name"
)

# Aggregate profiles
out_file = pathlib.PurePath(output_dir, f"{plate_name}.csv.gz")
ap = AggregateProfiles(sql_file=sql_file, strata=strata, operation=aggregate_method)
ap.aggregate_profiles(
    output_file=out_file, float_format=float_format, compression="gzip"
)

# Count cells
count_file = pathlib.PurePath(cell_count_dir, f"{plate_name}_cell_count.csv")
cell_count_df = ap.count_cells()
cell_count_df.to_csv(count_file, sep=",", index=False)

del ap

# Annotate profiles - Level 3 Data
anno_file = pathlib.PurePath(output_dir, f"{plate_name}_augmented.csv.gz")
anno_df = annotate(
    profiles=out_file,
    platemap=platemap_file,
    join_on=["Metadata_well_position", "Image_Metadata_Well"],
    cell_id=cell_id,
    format_broad_cmap=True,
    perturbation_mode="chemical",
    external_metadata=moa_df,
    external_join_left=["Metadata_broad_sample"],
    external_join_right=["Metadata_broad_sample"],
)

# Rename columns
anno_df = anno_df.rename(
    {"Image_Metadata_Plate": "Metadata_Plate", "Image_Metadata_Well": "Metadata_Well"},
    axis="columns",
)

# Add barcode platemap info
anno_df = anno_df.assign(
    Metadata_Assay_Plate_Barcode=barcode_platemap_df.Assay_Plate_Barcode.values[0],
    Metadata_Plate_Map_Name=barcode_platemap_df.Plate_Map_Name.values[0],
    Metadata_Batch_Number=barcode_platemap_df.Batch_Number.values[0],
    Metadata_Batch_Date=barcode_platemap_df.Batch_Date.values[0],
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
    compression=compression,
)

# Normalize Profiles (DMSO Control) - Level 4A Data
norm_dmso_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized_dmso.csv.gz")
normalize(
    profiles=anno_df,
    samples="Metadata_broad_sample == 'DMSO'",
    method=norm_method,
    output_file=norm_dmso_file,
    float_format=float_format,
    compression=compression,
)

# Normalize Profiles (Whole Plate) - Level 4A Data
norm_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized.csv.gz")
normalize(
    profiles=anno_df,
    samples="all",
    method=norm_method,
    output_file=norm_file,
    float_format=float_format,
    compression=compression,
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
    compression=compression,
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
    compression=compression,
)
