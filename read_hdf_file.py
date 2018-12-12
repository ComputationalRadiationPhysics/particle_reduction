import h5py
import re

class Particles_groups():
    """ Collect values from datasets in hdf file """

    def __init__(self, particles_name):
        self.particles_groups = []
        self.positions = []
        self.name_particles = particles_name
        self.momentum = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            name_idx = node.name.find(self.name_particles)
            if name_idx != -1:
                group_particles_name = node.name[name_idx + len(self.name_particles) + 1:]
                if group_particles_name.endswith('position'):
                    print('position   ' + group_particles_name)
                    self.positions.append(node)

                if group_particles_name.endswith('momentum'):
                    print('momentum  ' + group_particles_name)
                    self.momentum.append(node)
        return None

