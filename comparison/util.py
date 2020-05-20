import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


def build_file_dictionary(base_dir, tool="pycytominer"):
    file_match = {}
    if tool == "pycytominer":
        file_match["level_3"] = "augmented.csv.gz"
        file_match["level_4a"] = "normalized_dmso.csv.gz"
        file_match["level_4b"] = "normalized_feature_select_dmso.csv.gz"
        file_match["pycytominer_select"] = "normalized_feature_select_dmso.csv.gz"

    elif tool == "cytominer":
        file_match["level_3"] = "augmented.csv"
        file_match["level_4a"] = "normalized.csv"
        file_match["level_4b"] = "normalized_variable_selected.csv"
        file_match["pycytominer_select"] = "normalized.csv"

    file_dict = {}
    for plate in base_dir.iterdir():
        plate_name = plate.name
        if plate_name == ".DS_Store":
            continue
        file_dict[plate_name] = {}
        for plate_file in plate.iterdir():
            plate_file_name = plate_file.name
            if file_match["level_3"] in plate_file_name:
                file_dict[plate_name]["level_3"] = plate_file
            if file_match["level_4a"] in plate_file_name:
                file_dict[plate_name]["level_4a"] = plate_file
            if file_match["level_4b"] in plate_file_name:
                file_dict[plate_name]["level_4b"] = plate_file
            if file_match["pycytominer_select"] in plate_file_name:
                file_dict[plate_name]["pycytominer_select"] = plate_file

    return file_dict


def load_data(
    plate,
    pycyto_dict,
    cyto_dict,
    level,
    round_decimals,
    well_col="Metadata_Well",
    plate_col="Metadata_Plate",
    sample_col="Metadata_broad_sample",
):
    # Extract file from file dictionary
    pycyto_file = pycyto_dict[plate][level]
    cyto_file = cyto_dict[plate][level]

    # Load data
    pycyto_df = pd.read_csv(pycyto_file)
    cyto_df = pd.read_csv(cyto_file)

    # Confirm metadata are aligned
    pd.testing.assert_series_equal(pycyto_df.loc[:, well_col], cyto_df.loc[:, well_col])
    pd.testing.assert_series_equal(
        pycyto_df.loc[:, plate_col], cyto_df.loc[:, plate_col]
    )
    pd.testing.assert_series_equal(
        pycyto_df.loc[:, sample_col], cyto_df.loc[:, sample_col]
    )

    # Align to CP Features only
    pycyto_features = infer_cp_features(pycyto_df)
    cyto_features = infer_cp_features(cyto_df)

    # Features must be the same before feature selection
    if level in ["level_3", "level_4a"]:
        assert set(pycyto_features) == set(cyto_features), "features should be aligned!"

    # Reindex and round data
    pycyto_df = pycyto_df.reindex(set(pycyto_features), axis="columns").round(
        round_decimals
    )
    cyto_df = cyto_df.reindex(set(cyto_features), axis="columns").round(round_decimals)

    # If we're testing pycytominer feature selection procedure,
    # align cyto data with pycyto features
    if level == "pycytominer_select":
        cyto_df = cyto_df.reindex(set(pycyto_features), axis="columns")

    # Return a tuple of (pycyto data, cyto data) with aligned feature indices
    return (pycyto_df, cyto_df)


def build_filenames(output_dir, level, metrics=["median", "mean", "sum"]):
    output_files = {}
    for metric in metrics:
        output_files[metric] = pathlib.Path(
            f"{output_dir}/comparison_result_{level}_{metric}.tsv.gz"
        )

    return output_files
