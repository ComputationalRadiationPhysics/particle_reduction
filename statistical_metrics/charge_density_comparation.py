
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
