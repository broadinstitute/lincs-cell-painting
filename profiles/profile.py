"""
Perform the image-based profiling pipeline to process data
"""

import os
import pathlib
import pandas as pd
from pycytominer.aggregate import AggregateProfiles
from pycytominer.annotate import annotate
from pycytominer.normalize import normalize
from pycytominer.feature_select import feature_select
from pycytominer.audit import audit
from pycytominer.cyto_utils.output import output

from utils import get_args

# Load Command Line Arguments
args = get_args()

sql_file = args.sql_file
batch = args.batch
plate_name = args.plate_name
platemap_file = args.platemap_file
barcode_platemap_file = args.barcode_platemap_file
moa_file = args.moa_file
output_dir = args.output_dir
cell_count_dir = args.cell_count_dir

# Initialize profile processing
os.makedirs(output_dir, exist_ok=True)
os.makedirs(cell_count_dir, exist_ok=True)
compression = "gzip"
float_format = "%.5g"
strata = ["Image_Metadata_Plate", "Image_Metadata_Well"]
feature_select_ops = [
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
    "blacklist",
]
moa_df = pd.read_csv(moa_file, sep='\t')
barcode_platemap_df = pd.read_csv(barcode_platemap_file)
external_metadata_df = moa_df.merge(barcode_platemap_df)

# Step 1 - Aggregate profiles
out_file = pathlib.PurePath(output_dir, f"{plate_name}.csv.gz")
ap = AggregateProfiles(sql_file=sql_file, strata=strata)
ap.aggregate_profiles(output_file=out_file, float_format=float_format, compression="gzip")

# Step 2 - Count cells
count_file = pathlib.PurePath(cell_count_dir, f"{plate_name}_cell_count.tsv")
cell_count_df = ap.count_cells()
cell_count_df.assign(platemap=plate_name).to_csv(count_file, sep="\t")

del ap

# Step 3 - Annotate profiles
anno_file = pathlib.PurePath(output_dir, f"{plate_name}_augmented.csv.gz")
anno_df = annotate(
    profiles=out_file,
    platemap=platemap_file,
    cell_id="A549",
    format_broad_cmap=True,
    perturbation_mode="chemical",
    external_metadata=external_metadata_df,
    external_join_left=
    external_join_right=
    join_on=["Metadata_well_position", "Image_Metadata_Well"],
)

# Realign annotation metadata columns




# Output annotated file
output(df=anno_df, output_filename=anno_file, float_format=float_format, compression="gzip")

# Step 3A - Normalize Profiles (DMSO Control)
norm_dmso_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized_dmso.csv.gz")
normalize(
    profiles=anno_df,
    samples="Metadata_broad_sample == 'DMSO'",
    output_file=norm_dmso_file,
    float_format=float_format,
    compression="gzip"
)

# Step 3B - Normalize Profiles (Whole Plate)
norm_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized.csv.gz")
normalize(
    profiles=anno_df,
    samples="all",
    output_file=norm_file,
    float_format=float_format,
    compression="gzip"
)

# Step 4A - Feature Selection (DMSO Control)
feat_dmso_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized_feature_select_dmso.csv.gz")
feature_select(
    profiles=norm_dmso_file,
    features="infer",
    operation=feature_select_ops,
    output_file=feat_dmso_file,
)

# Step 4B - Feature Selection (Whole Plate)
feat_dmso_file = pathlib.PurePath(output_dir, f"{plate_name}_normalized_feature_select.csv.gz")
feature_select(
    profiles=norm_file,
    features="infer",
    operation=feature_select_ops,
    output_file=feat_dmso_file,
)
