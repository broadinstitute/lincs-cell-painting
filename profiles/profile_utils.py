import os
import pathlib
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sql_file", help="the filename of the sqlite file")
    parser.add_argument("-b", "--batch", help="string indicating the batch name")
    parser.add_argument("-p", "--plate_name", help="string indicating the platename")
    parser.add_argument("-f", "--platemap_file", help="path of platmap identifiers")
    parser.add_argument("-a", "--barcode_platemap_file", help="path of plate info")
    parser.add_argument("-m", "--moa_file", help="path of moa/target data storing")
    parser.add_argument("-o", "--output_dir", help="the directory to output the files")
    parser.add_argument("-c", "--cell_count_dir", help="directory to save cell counts")
    args = parser.parse_args()

    return args


def get_pipeline_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="reprocess and overwrite all data",
    )
    args = parser.parse_args()

    return args


def find_incomplete_plates(
    plates, output_dir="backend", file_match="normalized_feature_select.csv.gz"
):
    incomplete_plates = []
    for plate in plates:
        plate_data_dir = pathlib.PurePath(f"{output_dir}/{plate}/")
        dir_contents = os.listdir(plate_data_dir)
        if not any([file_match in x for x in dir_contents]):
            incomplete_plates.append(plate)
    return incomplete_plates
