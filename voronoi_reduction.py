""" Voronoi reduction algorithm """

import argparse
import os
import h5py
import re
import numpy
import math

dict_parametrs_names = {'position/x', 'position/y', 'position/z', 'momentum/x', 'momentum/y', 'momentum/z'}


class Dimentions_data():

    def __init__(self, name, vector):
        self.vector = numpy.array(vector)
        self.name = name
        self.standard_deviation = -1.
        self.lenght = -1.
        self.weights = []
        self.sum_weights = -1.


    def set_weights(self, weights):
        self.weights = weights
        self.sum_weights = numpy.sum(weights)

    def get_statistical_average(self):
        vector_mult = self.vector * self.weights
        sum_mult = numpy.sum(vector_mult)
        normalised_average = sum_mult / self.sum_weights
        return normalised_average

    def get_standard_deviation(self):
        average_value = self.get_statistical_average()

        sum_differences = 0.

        for i in range(0, len(self.weights)):
            sum_differences += self.weights[i] * math.pow((self.vector[i] - average_value), 2)

        normalised_values = math.sqrt(1 / self.sum_weights * sum_differences)

        return normalised_values

    def get_lenght(self):
        max_value = max(self.vector)
        min_value = min(self.vector)

        return max_value - min_value

    def get_coefficient_variation(self):
        deviation = self.get_standard_deviation()
        lenght = self.get_lenght()

        return deviation/lenght


class Particles_data():

    def __init__(self, name, data):
        self.name_particles = name
        self.data = data


class Parametrs_reader():

    def __init__(self):
        self.data = {}
        self.weights = []


    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):
            for name in dict_parametrs_names:
                if node.name.endswith(name):

                    current_data = Dimentions_data(name, node.value)
                    self.data[name] = current_data
                    break

            if node.name.endswith('weighting'):
                self.weights = node.value

    def get_parametrs(self):
        for demention in self.data:
            self.data[demention].set_weights(self.weights)


class Particles_groups():
    """ Collect values from datasets in hdf file """

    def __init__(self, particles_name):
        self.particles_groups = []
        self.positions = []
        self.name_particles = particles_name

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            print('is instance group ')
            name_idx = node.name.find(self.name_particles)
            if name_idx != -1:
                group_particles_name = node.name[name_idx + len(self.name_particles) + 1:]
                if group_particles_name.endswith('position'):
                    self.positions.append(node)
                if group_particles_name.find('/') == -1 and len(group_particles_name) != 0:
                    self.particles_groups.append(node)
        return None


def get_particles_name(hdf_file):
    """ Get name of particles group """

    particles_name = ''
    if hdf_file.attrs.get('particlesPath') != None:
        particles_name = hdf_file.attrs.get('particlesPath')
        particles_name = decode_name(particles_name)
    else:
        particles_name = 'particles'
    return particles_name


def decode_name(attribute_name):
    """ Decode name from binary """

    decoding_name = attribute_name.decode('ascii', errors='ignore')
    decoding_name = re.sub(r'\W+', '', decoding_name)
    return decoding_name


def collect_cell_parametrs(hdf_file_name):

    hdf_file = h5py.File(hdf_file_name)
    particles_name = get_particles_name(hdf_file)
    hdf_datasets = Particles_groups(particles_name)
    hdf_file.visititems(hdf_datasets)

    links_array = []

    for key in hdf_datasets.particles_groups:
        particles_parametrs = Parametrs_reader()
        key.visititems(particles_parametrs)
        data = particles_parametrs.get_parametrs()
        current_particles_data = Particles_data(key.name, data)
        links_array.append(current_particles_data)

    return links_array


def calculate_averages(weights, cell_values, sum_weights):

    vector_mult = cell_values * weights
    sum_mult = numpy.sum(vector_mult)
    normalised_coordinate = sum_mult / sum_weights
    return normalised_coordinate


def calculate_standard_deviation(sum_weights, average_value, weights, cell_values):

    sum_differences = 0.

    for i in range(0, len(weights)):
        sum_differences += weights[i] * math.pow((cell_values[i] - average_value), 2.)

    normalised_values = math.sqrt(1/sum_weights * sum_differences)

    return normalised_values



def devide_cells(parameter_max_dimention, particles):

    first_hyperline = particles
    secound_hyperline = particles
    first_hyperline[parameter_max_dimention].get_first_hyperiline()
    secound_hyperline[parameter_max_dimention].get_second_hyperiline()
    return first_hyperline, secound_hyperline


def voronoi_algorithm(hdf_file, hdf_file_reduction, tolerance_momentum, tolerance_position):

    links_array = collect_cell_parametrs(hdf_file)


    print(' hdf file: ' + str(hdf_file))
    print(' hdf file reduction ' + str(hdf_file_reduction))
    print(' tolerance momentum  ' + str(tolerance_momentum))
    print(' tolerance position ' + str(tolerance_position))


def voronoi_reduction(hdf_file, hdf_file_reduction, tolerance_momentum, tolerance_position):

    name_hdf_file_reduction = ''

    if hdf_file != '':
        if os.path.exists(hdf_file):
            name = hdf_file[:-4]
            idx_of_name = name.rfind('/')
            if idx_of_name != -1:
                name_hdf_file_reduction = hdf_file_reduction + hdf_file[idx_of_name + 1: -6] + 'reduction.h5'
            else:
                name_hdf_file_reduction = hdf_file_reduction + hdf_file[:-3] + '.h5'

            voronoi_algorithm(hdf_file, hdf_file_reduction, tolerance_momentum, tolerance_position)
        else:
            print('The .hdf file does not exist')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="voronoi reduction")

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="reducted hdf file")

    parser.add_argument("-momentum_tol", metavar='tolerance_momentum', type=float,
                        help="tolerance of momentum")

    parser.add_argument("-momentum_pos", metavar='tolerance_position', type=float,
                        help="tolerance of position")

    args = parser.parse_args()
    voronoi_reduction(args.hdf, args.hdf_re, args.momentum_tol, args.momentum_pos)

