import read_hdf_file
import h5py
import numpy
import math


class WeightReader():

    def __init__(self):
        self.weight = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):
            if node.name.endswith('weighting'):
                self.weight = node.value


def compute_weight_sum(weight, values):

    dimension = values.get_dimension()

    sum_x = 0.
    sum_y = 0.
    sum_z = 0.

    if dimension == 3:
        size_of_position = len(values.vector_x)
        for i in range(0, size_of_position):
            sum_x += values.vector_x[i] * weight[i]
            sum_y += values.vector_y[i] * weight[i]
            sum_z += values.vector_z[i] * weight[i]

    elif dimension == 2:
        size_of_position = len(values.vector_x)
        for i in range(0, size_of_position):
            sum_x += values.vector_x[i] * weight[i]
            sum_y += values.vector_y[i] * weight[i]

    return sum_x, sum_y, sum_z


def get_dataset_values(group, name_dataset):

    values = read_hdf_file.Dataset_Reader(name_dataset)
    group.visititems(values)
    weight_reader = WeightReader()
    group.visititems(weight_reader)

    return weight_reader.weight, values


def count_weight_difference(weight_first, values_first, weight_second, values_second):

    sum_first = compute_weight_sum(weight_first, values_first)

    sum_second = compute_weight_sum(weight_second, values_second)

    for i in range(0, values_first.get_dimension()):
        relative_error = (sum_first[i] - sum_second[i]) / sum_first[i]
        print(relative_error)
        assert math.fabs(relative_error) < 1e-6, "Big relative error, reduction is wrong"

    return sum_first, sum_second


def compute_standard_deviation(weights, coords):

    sum_weights = numpy.sum(weights)
    sum_coords = 0.
    for i in range(0, len(weights)):
        sum_coords += weights[i] * coords[i]

    average_value = sum_coords / sum_weights

    sum_sq = 0.

    for i in range(0, len(weights)):
        sum_sq += (weights[i] * coords[i] - average_value) * (weights[i] * coords[i] - average_value)

    norm_sq = sum_sq / sum_weights

    norm_sq = math.sqrt(norm_sq)

    return norm_sq


def compute_momentum_standart_deviation(weights, momentum_values):

    deviation_values = []
    if momentum_values.get_dimension() == 3:
        deviation_values.append(compute_standard_deviation(weights, momentum_values.vector_x))
        deviation_values.append(compute_standard_deviation(weights, momentum_values.vector_y))
        deviation_values.append(compute_standard_deviation(weights, momentum_values.vector_z))

    if momentum_values.get_dimension() == 2:
        deviation_values.append(compute_standard_deviation(weights, momentum_values.vector_x))
        deviation_values.append(compute_standard_deviation(weights, momentum_values.vector_y))

    return deviation_values


def compare_weight_coordinates(first_hdf_file_name, second_hdf_file_name):

    first_hdf_file = h5py.File(first_hdf_file_name, 'a')
    second_hdf_file = h5py.File(second_hdf_file_name, 'a')

    particles_name_first = read_hdf_file.get_particles_name(first_hdf_file)
    particles_groups_first = read_hdf_file.ParticlesGroups(particles_name_first)
    first_hdf_file.visititems(particles_groups_first)

    particles_name_second = read_hdf_file.get_particles_name(second_hdf_file)
    particles_groups_second = read_hdf_file.ParticlesGroups(particles_name_second)
    second_hdf_file.visititems(particles_groups_second)

    size_groups = len(particles_groups_second.particles_groups)

    for i in range(0, size_groups):
        weight_first, positions_first = get_dataset_values(particles_groups_first.particles_groups[i], 'position')
        weight_first, momentum_first = get_dataset_values(particles_groups_first.particles_groups[i], 'momentum')

        weight_second, positions_second = get_dataset_values(particles_groups_second.particles_groups[i], 'position')
        weight_second, momentum_second = get_dataset_values(particles_groups_second.particles_groups[i], 'momentum')

        count_weight_difference(weight_first, positions_first, weight_second, positions_second)
        count_weight_difference(weight_first, momentum_first, weight_second, momentum_second)



def test_function():

    compare_weight_coordinates()


