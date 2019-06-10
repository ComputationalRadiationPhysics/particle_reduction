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


def norm_array01(vector, max_value, min_value):
    array_norm = []
    for point in vector:
        value = (point - min_value) /(max_value - min_value)
        array_norm.append(value)
    return array_norm



def compute_kernel_2d(array_x, array_y, weights_values):

    copy_weights = copy.deepcopy(weights_values)

    x_min = numpy.min(array_x)
    x_max = numpy.max(array_x)
    array_x_norm = norm_array01(array_x, x_min, x_max)

    y_min = numpy.min(array_y)
    y_max = numpy.max(array_y)
    array_y_norm = norm_array01(array_y, y_min, y_max)

    x_min = 0.
    x_max = 1.

    y_min = 0.
    y_max = 1.

    X, Y = numpy.mgrid[x_min:x_max:50j, y_min:y_max:50j]
    lover_values = [x_min, y_min]
    upper_values = [x_max, y_max]

    positions = numpy.vstack([X.ravel(), Y.ravel()])
    values = numpy.vstack([array_x_norm, array_y_norm])
    kernel = stats.gaussian_kde(values, weights=copy_weights)
    value_intergrate = kernel.integrate_box(lover_values, upper_values)
    array_kernell = kernel(positions)

    return array_kernell



def compute_kernel_3d(array_x, array_y, array_z, weights_values):
    copy_weights = copy.deepcopy(weights_values)

    x_min = numpy.min(array_x)
    x_max = numpy.max(array_x)

    array_x_norm = norm_array01(array_x, x_min, x_max)
    y_min = numpy.min(array_y)
    y_max = numpy.max(array_y)

    array_y_norm = norm_array01(array_y, y_min, y_max)

    z_min = numpy.min(array_z)
    z_max = numpy.max(array_z)

    array_z_norm = norm_array01(array_z, z_min, z_max)
    x_min = 0.
    x_max = 1.

    y_min = 0.
    y_max = 1.

    z_min = 0.
    z_max = 1.

    X, Y, Z = numpy.mgrid[x_min:x_max:50j, y_min:y_max:50j, z_min:z_max:50j ]
    lover_values = [x_min, y_min, z_min]
    upper_values = [x_max, y_max, z_max]

    positions = numpy.vstack([X.ravel(), Y.ravel(), Z.ravel()])

    values = numpy.vstack([array_x_norm, array_y_norm, array_z_norm])
    kernel = stats.gaussian_kde(values, weights=copy_weights)
#    value_intergrate = kernel.integrate_box(lover_values, upper_values)
    array_kernell = kernel(positions)

    return array_kernell


def compute_position_kernel(absolute_coordinates, dimensions_first, weights):


    if dimensions_first == 3:
        return compute_kernel_3d(absolute_coordinates[:, 0], absolute_coordinates[:, 1], absolute_coordinates[:, 2], weights)

    else:
        return compute_kernel_2d(absolute_coordinates[:, 0], absolute_coordinates[:, 1], weights)


def write_values_into_csv_file(metric_values, csv_file_name):

    row = []

    for value in metric_values:
        row.append(str(value))

    with open(csv_file_name, 'a') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(row)

    csvFile.close()

def base_corparation(first_hdf_file_name, second_hdf_file_name, csv_file_name):


    #first_hdf_file_name = "/home/kseniia/Documents/measuring_qality/checkpoint_0.h5"
    #second_hdf_file_name = "/home/kseniia/Documents/measuring_qality/random_reduction_0.1.h5"

    first_hdf_file = h5py.File(first_hdf_file_name, 'a')
    second_hdf_file = h5py.File(second_hdf_file_name, 'a')

    particles_name_first = read_hdf_file.get_particles_name(first_hdf_file)

    particles_groups_first = read_hdf_file.ParticlesGroups(particles_name_first)
    first_hdf_file.visititems(particles_groups_first)

    particles_name_second = read_hdf_file.get_particles_name(second_hdf_file)

    particles_groups_second = read_hdf_file.ParticlesGroups(particles_name_second)
    second_hdf_file.visititems(particles_groups_second)

    idx = first_hdf_file_name.rfind('/')
    name_of_file = first_hdf_file_name[idx + 1:len(first_hdf_file_name) - 3]


    for i in range(0, len(particles_groups_first.particles_groups)):

        name_of_group = particles_groups_first.particles_groups[i].name
        group_first = particles_groups_first.particles_groups[i]
        group_second = particles_groups_second.particles_groups[i]
        print(name_of_group)

        idx = name_of_group.rfind('/')
        substr = name_of_group[idx + 1:len(name_of_group)]
        name_of_iteration = substr + '_' + name_of_file
        print(name_of_iteration)

        absolute_coordinates_first, dimensions_first, weights_first = read_hdf_file.get_relative_coordinates(group_first)
        absolute_coordinates_secound, dimensions_secound, weights_secound = read_hdf_file.get_relative_coordinates(group_second)

        if len(absolute_coordinates_first) == 0:
            continue

        values_kde_first = compute_position_kernel(absolute_coordinates_first, dimensions_first, weights_first)

        values_kde_second = compute_position_kernel(absolute_coordinates_secound, dimensions_secound, weights_secound)

        print('-----------------------------------')
        difference_values = []
        for i in range(0, len(values_kde_first)):
            difference_values.append(abs(values_kde_first[i] - values_kde_second[i]))
        max_difference = numpy.max(difference_values)

        eucludian_values = []
        for i in range(0, len(values_kde_first)):
            eucludian_values.append(
                (values_kde_first[i] - values_kde_second[i]) * (values_kde_first[i] - values_kde_second[i]))

        sum_of_values = numpy.sum(eucludian_values) / len(eucludian_values)
        # print('sum_of_values  ' + str(sum_of_values))
        sqrt_value = math.sqrt(sum_of_values)
        print('my euclidian metric ' + str(sqrt_value))

        dist = sp.distance.euclidean(values_kde_first, values_kde_second)

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



        vector_values = [name_of_iteration, str(dist), str(minkowski), str(correlation), str(braycurtis),
                               str(cityblock), str(cosine)]

        write_values_into_csv_file(vector_values, csv_file_name)
