import sys
sys.path.append("../")
import read_hdf_file
import scipy.spatial as sp
import h5py
import csv
import numpy
from scipy import stats
import copy
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


def compute_norm_eu_values(values_kde_first, values_kde_second):

    eucludian_values = []
    for i in range(0, len(values_kde_first)):
        eucludian_values.append(
            (values_kde_first[i] - values_kde_second[i]) * (values_kde_first[i] - values_kde_second[i]))

    sum_of_values = numpy.sum(eucludian_values) / len(eucludian_values)

    sqrt_value = math.sqrt(sum_of_values)

    return sqrt_value


def compute_stats_metrics(values_kde_first, values_kde_second):

    print("values_kde_first "+str(values_kde_first))
    print("values_kde_second " + str(values_kde_second))
    if len(values_kde_first) == 0 or len(values_kde_second) == 0:
        return 0, 0, 0, 0, 0, 0, 0, 0

    eu_metric = compute_norm_eu_values(values_kde_first, values_kde_second)
    print('my euclidian metric ' + str(eu_metric))

    max_difference = compute_max_differences(values_kde_first, values_kde_second)
    print('max diff' + str(max_difference))
    print('-----------------------------------')
    minkowski = sp.distance.minkowski(values_kde_first, values_kde_second)
    print('minkovski')
    print(minkowski)

    print('cosine')
    cosine = sp.distance.cosine(values_kde_first, values_kde_second)
    print(cosine)

    print('chebyshev')
    chebyshev = sp.distance.chebyshev(values_kde_first, values_kde_second)
    print(chebyshev)

    print('correlation')
    correlation = sp.distance.correlation(values_kde_first, values_kde_second)
    print(correlation)

    print('braycurtis')
    braycurtis = sp.distance.braycurtis(values_kde_first, values_kde_second)
    print(braycurtis)

    print('cityblock')
    cityblock = sp.distance.cityblock(values_kde_first, values_kde_second)
    print(cityblock)

    return eu_metric, max_difference, minkowski, cosine, chebyshev, correlation, braycurtis, cityblock


def compute_1d_kernels(idx_start, idx_end, values, weights):


    kernel_x = None
    kernel_y = None
    kernel_z = None

    dimension = int(idx_end - idx_start)


    if dimension == 3:
        kernel_x = compute_kernel_1d(values[:, idx_start], weights)
        kernel_y = compute_kernel_1d(values[:, idx_start + 1], weights)
        kernel_z = compute_kernel_1d(values[:, idx_start + 2], weights)
    elif dimension == 2:
        kernel_x = compute_kernel_1d(values[:, idx_start], weights)
        kernel_y = compute_kernel_1d(values[:, idx_start + 1], weights)

    return kernel_x, kernel_y, kernel_z


def compute_1_d_stats_metrics(kernel_coords_first_x, kernel_coords_first_y, kernel_coords_first_z,
                              kernel_coords_second_x, kernel_coords_second_y, kernel_coords_second_z, file_name, name_metrics):

    result_rows = []

    eu_metric_x, max_difference_x, minkowski_x, cosine_x, chebyshev_x, correlation_x, braycurtis_x, cityblock_x = \
        compute_stats_metrics(kernel_coords_first_x, kernel_coords_second_x)

    result_name = file_name + name_metrics + "X"
    vector_values = [result_name, str(eu_metric_x), str(max_difference_x), str(minkowski_x), str(cosine_x),
                     str(chebyshev_x), str(correlation_x), str(braycurtis_x), str(cityblock_x)]

    result_rows.append(vector_values)

    result_name = file_name + name_metrics + "Y"
    eu_metric_y, max_difference_y, minkowski_y, cosine_y, chebyshev_y, correlation_y, braycurtis_y, cityblock_y = \
        compute_stats_metrics(kernel_coords_first_y, kernel_coords_second_y)

    vector_values = [result_name, str(eu_metric_y), str(max_difference_y), str(minkowski_y), str(cosine_y),
                     str(chebyshev_y), str(correlation_y), str(braycurtis_y), str(cityblock_y)]

    result_rows.append(vector_values)

    if kernel_coords_first_z != None and kernel_coords_second_z != None:

        result_name = file_name + name_metrics + "Z"
        eu_metric_z, max_difference_z, minkowski_z, cosine_z, chebyshev_z, correlation_z, braycurtis_z, cityblock_z = \
            compute_stats_metrics(kernel_coords_first_z, kernel_coords_second_z)

        vector_values = [result_name, str(eu_metric_z), str(max_difference_z),
                         str(minkowski_z), str(cosine_z), str(chebyshev_z),
                         str(correlation_z), str(braycurtis_z), str(cityblock_z)]
        result_rows.append(vector_values)

    return result_rows


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
