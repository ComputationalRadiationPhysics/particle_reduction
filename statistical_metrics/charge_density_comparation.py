
import h5py
import numpy
import math
import argparse
import csv
class get_fields():
    """

    Collect particles groups from hdf file
    particles_name -- name of main partcles group

    """

    def __init__(self):
        self.fields = []

    def __call__(self, name, node):
        if node.name.endswith('fields'):
            self.fields = node
        return None


class Read_density():
    def __init__(self):
        self.charge_density = {}
        self.density = {}
        self.energy_density = {}
        self.name_particles = 'fields'

    def get_name_group(self, node_name):
        name_idx = node_name.find(self.name_particles)
        group_particles_name = ''
        if name_idx != -1:
            group_particles_name = node_name[name_idx + len(self.name_particles) + 1:]

        return group_particles_name

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):
            name_group = self.get_name_group(node.name)
            if node.name.endswith('_chargeDensity'):
                self.charge_density[name_group] = node.value
            if node.name.endswith('_density'):
                self.density[name_group] = node.value
            if node.name.endswith('_energyDensity'):
                self.energy_density[name_group] = node.value


def compute_difference(values_first, values_secound):

    difference_values = []
    for i in range(0, len(values_first)):
        for j in range(0, len(values_first[0])):
            difference_values.append(abs(values_first[i][j] - values_secound[i][j]))
    max_difference = numpy.max(difference_values)

    return max_difference


def compute_eucledian(values_first, values_secound):

    eucludian_values = []
    for i in range(0, len(values_first)):
        for j in range(0, len(values_first[0])):
            eucludian_values.append(
                (values_first[i][j] - values_secound[i][j]) * (values_first[i][j] - values_secound[i][j]))

    sum_of_values = numpy.sum(eucludian_values) / len(eucludian_values)

    sqrt_value = math.sqrt(sum_of_values)
    return sqrt_value


def write_values_into_csv_file(metric_values, csv_file_name):

    row = []

    for value in metric_values:
        row.append(str(value))

    with open(csv_file_name, 'a') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(row)

    csvFile.close()


def compute_difference_charge_density(type_electrons, charge_density_first, charge_density_secound):

    name_values = type_electrons + '_chargeDensity'
    values_first = charge_density_first[name_values]
    values_secound = charge_density_secound[name_values]
    max_difference = compute_difference(values_first, values_secound)
    eucledian_values = compute_eucledian(values_first, values_secound)
    return max_difference, eucledian_values

