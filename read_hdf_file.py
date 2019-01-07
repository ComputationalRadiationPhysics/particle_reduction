import h5py
import re

class Particles_functor():
    """ Collect values from datasets in hdf file """

    def __init__(self):
        self.particles_groups = []
        self.positions = []
        self.momentum = []
        self.mass = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            if node.name.endswith('position'):
                self.positions.append(node)

            if node.name.endswith('momentum'):
                self.momentum.append(node)

            if node.name.endswith('mass'):
                self.mass = node.attrs['value']
        return None


class Particles_groups():
    def __init__(self, particles_name):

        self.name_particles = particles_name
        self.particles_groups = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            name_idx = node.name.find(self.name_particles)
            if name_idx != -1:
                group_particles_name = node.name[name_idx + len(self.name_particles) + 1:]
                if group_particles_name.find('/') == -1 and group_particles_name != '':
                    self.particles_groups.append(node)
        return None

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):

            dataset_x = self.name_dataset + '/x'
            dataset_y = self.name_dataset + '/y'
            dataset_z = self.name_dataset + '/z'

            if node.name.endswith('position/z'):
                self.z_coord = node.value

        return None


class Dataset_reader():
    """ Collect values from datasets in hdf file """

    def __init__(self, name_dataset):
        self.vector_x = []
        self.vector_y = []
        self.vector_z = []
        self.name_dataset = name_dataset

    def __call__(self, name, node):

        dataset_x = self.name_dataset + '/x'
        dataset_y = self.name_dataset + '/y'
        dataset_z = self.name_dataset + '/z'
        if isinstance(node, h5py.Dataset):
            if node.name.endswith(dataset_x):
                self.vector_x = node.value

            if node.name.endswith(dataset_y):
                self.vector_y = node.value

            if node.name.endswith(dataset_z):
                self.vector_z = node.value

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

