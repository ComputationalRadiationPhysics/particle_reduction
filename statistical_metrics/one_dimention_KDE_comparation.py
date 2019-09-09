import sys
sys.path.append("../")
import read_hdf_file
import csv
import numpy
import math
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


def compute_kernel_1d(array_x, weights):

    copy_weights = copy.deepcopy(weights)

    x_min = numpy.min(array_x)
    x_max = numpy.max(array_x)
    array_kernell = []
    if abs(x_max - x_min) == 0:
        return array_kernell

    array_x_norm = norm_array01(array_x, x_min, x_max)

    x_min = 0.
    x_max = 1.

    X = numpy.mgrid[x_min:x_max:250j]

    positions = numpy.vstack([X.ravel()])
    values = numpy.vstack([array_x_norm])
    kernel = stats.gaussian_kde(values, weights=copy_weights)
    array_kernell = kernel(positions)

    return array_kernell


def read_group_values(particles_groups, idx):

    group_first = particles_groups.particles_groups[idx]

    data_first, weights_first, dimensions_first, unit_si_position_first, unit_si_momentum_first \
        = read_hdf_file.read_points_group(group_first)

    if len(data_first) == 0:
        return [], [], []

    position_offset_first, unit_si_offset_first = read_hdf_file.read_position_offset(group_first)

    absolute_coordinates_first = read_hdf_file.get_absolute_coordinates(data_first, position_offset_first,
                                                                        unit_si_offset_first,
                                                                        unit_si_position_first, dimensions_first,
                                                                        unit_si_momentum_first)
    absolute_coordinates_first = numpy.array(absolute_coordinates_first)

    return absolute_coordinates_first, dimensions_first, weights_first


def compute_max_differences(values_kde_first, values_kde_second):

    difference_values = []
    for i in range(0, len(values_kde_first)):
        difference_values.append(abs(values_kde_first[i] - values_kde_second[i]))
    max_difference = numpy.max(difference_values)

    return max_difference


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
