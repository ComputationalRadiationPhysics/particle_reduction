import sys
import csv
import numpy
import argparse


def write_values_into_csv_file(metric_values, csv_file_name):

    row = []

    with open(csv_file_name, 'a') as csvFile:
        writer = csv.writer(csvFile)
        for value in metric_values:

            writer.writerow(value)

    csvFile.close()


def norm_array01(vector, max_value, min_value):
    array_norm = []
    for point in vector:
        if abs(max_value - min_value) == 0:
            array_norm.append(0.)
            continue

        value = (point - min_value) / (max_value - min_value)

        array_norm.append(value)
    return array_norm








if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="file_comparation")

    parser.add_argument("-first_hdf", metavar='first_hdf', type=str,
                        help="first hdf file to compare")

    parser.add_argument("-second_hdf", metavar='second_hdf', type=str,
                        help="second hdf file to compare")

    parser.add_argument("-csv_file", metavar='csv_file', type=str,
                        help="csv file")

    parser.add_argument("-file_indx", metavar='file_indx', type=str,
                        help="information about name iteration")

    args = parser.parse_args()

    base_corparation(args.first_hdf, args.second_hdf, args.csv_file, args.file_indx)
