import h5py
import re

class Particles_groups():
    """ Collect values from datasets in hdf file """

    def __init__(self, particles_name):
        self.particles_groups = []
        self.positions = []
        self.name_particles = particles_name
        self.momentum = []
        self.mass = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            name_idx = node.name.find(self.name_particles)
            if name_idx != -1:
                group_particles_name = node.name[name_idx + len(self.name_particles) + 1:]
                if group_particles_name.endswith('position'):
                    self.positions.append(node)

                if group_particles_name.endswith('momentum'):
                    self.momentum.append(node)

                if group_particles_name.endswith('mass'):
                    self.mass = node.attrs['value']
        return None


class Position_reader():
    

    def __init__(self):
        self.x_coord = []
        self.y_coord = []
        self.z_coord = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):
            print('name == ' + str(node.name))
            if node.name.endswith('position/x'):
                self.x_coord = node.value

            if node.name.endswith('position/y'):
                self.y_coord = node.value

            if node.name.endswith('position/z'):
                self.z_coord = node.value

        return None


class Momentum_reader():
    """ Collect values from datasets in hdf file """

    def __init__(self):
        self.x_momentum = []
        self.y_momentum = []
        self.z_momentum = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):
            print('name == ' + str(node.name))
            if node.name.endswith('momentum/x'):
                self.x_momentum = node.value

            if node.name.endswith('momentum/y'):
                self.y_momentum = node.value

            if node.name.endswith('momentum/z'):
                self.z_momentum = node.value

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

