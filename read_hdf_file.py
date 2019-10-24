import h5py
import re
import numpy

class Dimensions:
    def __init__(self, dimension_position, dimension_momentum):
        self.dimension_position = dimension_position
        self.dimension_momentum = dimension_momentum


class ParticlesFunctor():
    """

    Collect values(weighting, position, momentum) from paticle dataset in hdf file.
    positions -- group of position coords
    momentum -- group of momentum coords
    weightins -- values of weights for particles

    """

    def __init__(self):
        self.positions = []
        self.momentum = []
        self.weighting = []
        self.bound_electrons = []
        self.position_offset = []

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):
            if node.name.endswith('weighting'):
                self.weighting = node.value

            if node.name.endswith('boundElectrons'):
                self.bound_electrons = node.value

        if isinstance(node, h5py.Group):
            if node.name.endswith('position'):

                self.positions.append(node)

            if node.name.endswith('momentum'):
                self.momentum.append(node)

            if node.name.endswith('positionOffset'):
                self.position_offset.append(node)

        return None


class MassFunctor():
    """

    """
    def __init__(self):
        self.mass = []

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):
            if node.name.endswith('mass'):
                self.weighting = node.value

        if isinstance(node, h5py.Group):
            if node.name.endswith('mass'):
                self.mass = node.attrs['value']

        return None

class ParticlesGroups():
    """

    Collect particles groups from hdf file
    particles_name -- name of main partcles group

    """

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


class ReadPatchGroup():
    def __init__(self):
        self.patch_group = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Group):
            if node.name.endswith('particlePatches'):
                self.patch_group.append(node)


class ReadPatchValues():
    def __init__(self):
        self.numParticles = []
        self.numParticlesOffset = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):# numParticles
            if node.name.endswith('numParticles'):
                self.numParticles = node.value

            if node.name.endswith('numParticlesOffset'):
                self.numParticlesOffset = node.value


class points_writer():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, library_datasets, name_dataset):
        self.dataset_x = name_dataset + '/x'
        self.dataset_y = name_dataset + '/y'
        self.dataset_z = name_dataset + '/z'

        self.vector_x = library_datasets[self.dataset_x]
        self.vector_y = library_datasets[self.dataset_y]
        self.vector_z = library_datasets[self.dataset_z]
        self.hdf_file = hdf_file

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):

            if node.name.endswith(self.dataset_x):
                attributes = {}

                for attr in node.attrs.keys():
                    attributes[attr] = node.attrs[attr]
                node_name = node.name

                current_dtype = self.hdf_file[node.name].dtype
                del self.hdf_file[node.name]

                dset = self.hdf_file.create_dataset(node_name, data=self.vector_x, dtype=current_dtype)
                for attr in attributes:
                    dset.attrs[attr] = attributes[attr]

            elif node.name.endswith(self.dataset_y):
                attributes = {}
                for attr in node.attrs.keys():
                    attributes[attr] = node.attrs[attr]
                node_name = node.name
                current_dtype = self.hdf_file[node.name].dtype
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_y, dtype=current_dtype)
                for attr in attributes:
                    dset.attrs[attr] = attributes[attr]

            elif node.name.endswith(self.dataset_z):
                attributes = {}
                for attr in node.attrs.keys():
                    attributes[attr] = node.attrs[attr]
                node_name = node.name
                current_dtype = self.hdf_file[node.name].dtype
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_z, dtype=current_dtype)
                for attr in attributes:
                    dset.attrs[attr] = attributes[attr]

        return None


class vector_writer():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, vector_values, name_dataset, start_dimension):
        self.dataset_x = name_dataset + '/x'
        self.dataset_y = name_dataset + '/y'
        self.dataset_z = name_dataset + '/z'

        self.vector_values = vector_values
        self.hdf_file = hdf_file
        self.is_first_part = True
        self.start_dimension = start_dimension

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):
            if node.name.endswith(self.dataset_x):
                self.choose_writing_type(node, self.start_dimension)

            elif node.name.endswith(self.dataset_y):
                self.choose_writing_type(node, self.start_dimension + 1)

            elif node.name.endswith(self.dataset_z):
                self.choose_writing_type(node, self.start_dimension + 2)

    def choose_writing_type(self, node, type_idx):
        if self.is_first_part:
            self.first_part_writing(node, self.vector_values[:, type_idx])
        else:
            self.resize_writing(node, self.vector_values[:, type_idx])

    def first_part_writing(self, node, values):

        attributes = {}
        for attr in node.attrs.keys():
            attributes[attr] = node.attrs[attr]

        node_name = node.name
        current_dtype = self.hdf_file[node.name].dtype
        del self.hdf_file[node.name]
        dset = self.hdf_file.create_dataset(node_name, maxshape=(None,), data=values, dtype=current_dtype, chunks=True)
        for attr in attributes:
            dset.attrs[attr] = attributes[attr]

        return None

    def resize_writing(self, node, values):

        node_name = node.name
        self.hdf_file[node_name].resize((self.hdf_file[node_name].shape[0] + values.shape[0]), axis=0)
        self.hdf_file[node_name][-values.shape[0]:] = values



class dataset_writer():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, values, name_dataset):

        self.values = values
        self.hdf_file = hdf_file
        self.name_dataset = name_dataset

    def __call__(self, name, node):



        if isinstance(node, h5py.Dataset):

            if node.name.endswith(self.name_dataset):
                attributes = {}
                for attr in node.attrs.keys():
                    attributes[attr] = node.attrs[attr]
                node_name = node.name
                current_dtype = self.hdf_file[node.name].dtype
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.values, dtype=current_dtype)
                for attr in attributes:
                    dset.attrs[attr] = attributes[attr]
        return None



class PatchValuesWriter():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file,  numParticlesOffset, numParticles):

        self.numParticles = numParticles
        self.numParticlesOffset = numParticlesOffset
        self.hdf_file = hdf_file

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):

            if node.name.endswith('numParticles'):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.numParticles)

            if node.name.endswith('numParticlesOffset'):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.numParticlesOffset)
        return None


class DatasetReader():
    """

     Read datasets values from hdf file
     name_dataset -- name of base group

    """

    def __init__(self, name_dataset):
        self.vector_x = []
        self.vector_y = []
        self.vector_z = []
        self.unit_SI_x = 1.
        self.unit_SI_y = 1.
        self.unit_SI_z = 1.
        self.name_dataset = name_dataset

    def __call__(self, name, node):

        dataset_x = self.name_dataset + '/x'
        dataset_y = self.name_dataset + '/y'
        dataset_z = self.name_dataset + '/z'
        if isinstance(node, h5py.Dataset):
            if node.name.endswith(dataset_x):
                self.vector_x = node.value
                self.unit_SI_x = node.attrs["unitSI"]

            if node.name.endswith(dataset_y):
                self.vector_y = node.value
                self.unit_SI_y = node.attrs["unitSI"]

            if node.name.endswith(dataset_z):
                self.vector_z = node.value
                self.unit_SI_z = node.attrs["unitSI"]

        return None

    def get_unit_si_array(self):

        array_unit_SI = []
        if self.get_dimension() == 3:
            array_unit_SI = [self.unit_SI_x, self.unit_SI_y, self.unit_SI_z]
        elif self.get_dimension() == 2:
            array_unit_SI = [self.unit_SI_x, self.unit_SI_y]

        return array_unit_SI

    def get_dimension(self):
        """

         get dimension of particles datasets

        """

        size = 0
        if len(self.vector_x) > 0:
            if len(self.vector_y) > 0:
                if len(self.vector_z) > 0:
                    size = 3
                else:
                    size = 2
            else:
                size = 1

        return size


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


def create_points_array(coord_collection, momentum_collection):
    """

    create array of 2-d, 3-d points from datasets
    coord_collection -- datasets from hdf file

    """

    vector_coords = []

    dimension_coord = coord_collection.get_dimension()
    dimension_momentum = momentum_collection.get_dimension()



    if dimension_coord == 3 and dimension_momentum == 3:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y, coord_collection.vector_z,
                             momentum_collection.vector_x, momentum_collection.vector_y, momentum_collection.vector_z)]

    elif dimension_coord == 3 and dimension_momentum == 2:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y, coord_collection.vector_z,
                             momentum_collection.vector_x, momentum_collection.vector_y)]

    elif dimension_coord == 2 and dimension_momentum == 3:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y,
                             momentum_collection.vector_x, momentum_collection.vector_y, momentum_collection.vector_z)]

    elif dimension_coord == 2 and dimension_momentum == 2:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y,
                             momentum_collection.vector_x, momentum_collection.vector_y)]

    return vector_coords


def create_points_array_bound_electrons(coord_collection, momentum_collection, bound_electrons):
    """

    create array of 2-d, 3-d points from datasets
    coord_collection -- datasets from hdf file

    """

    vector_coords = []

    dimension_coord = coord_collection.get_dimension()
    dimension_momentum = momentum_collection.get_dimension()



    if dimension_coord == 3 and dimension_momentum == 3:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y, coord_collection.vector_z,
                             momentum_collection.vector_x, momentum_collection.vector_y, momentum_collection.vector_z,
                             bound_electrons)]

    elif dimension_coord == 3 and dimension_momentum == 2:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y, coord_collection.vector_z,
                             momentum_collection.vector_x, momentum_collection.vector_y, bound_electrons)]

    elif dimension_coord == 2 and dimension_momentum == 3:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y,
                             momentum_collection.vector_x, momentum_collection.vector_y, momentum_collection.vector_z,
                             bound_electrons)]

    elif dimension_coord == 2 and dimension_momentum == 2:
        vector_coords = [list(x) for x in
                         zip(coord_collection.vector_x, coord_collection.vector_y,
                             momentum_collection.vector_x, momentum_collection.vector_y, bound_electrons)]

    return vector_coords


def create_dataset_from_point_array(points, name_dataset):

    vector_x = []
    vector_y = []
    vector_z = []
    weighting = []

    dimension = len(points[name_dataset][0].coords)
    coordinates_dataset = points[name_dataset]
    for i in range(0, len(coordinates_dataset)):
        if dimension == 3:
            vector_x.append(coordinates_dataset[i].coords[0])
            vector_y.append(coordinates_dataset[i].coords[1])
            vector_z.append(coordinates_dataset[i].coords[2])
            weighting.append(coordinates_dataset[i].weight)

        if dimension == 2:
            vector_x.append(coordinates_dataset[i].coords[0])
            vector_y.append(coordinates_dataset[i].coords[1])

            weighting.append(coordinates_dataset[i].weight)

    return vector_x, vector_y, vector_z, weighting


def create_library_of_datasets(points):


    datasets = {}

    position_x, position_y, position_z, weighting = create_dataset_from_point_array(points, 'position')
    momentum_x, momentum_y, momentum_z, weighting = create_dataset_from_point_array(points, 'momentum')

    datasets['position/x'] = position_x
    datasets['position/y'] = position_y
    datasets['position/z'] = position_z

    datasets['momentum/x'] = momentum_x
    datasets['momentum/y'] = momentum_y
    datasets['momentum/z'] = momentum_z

    datasets['weighting'] = weighting

    return datasets


def create_points_library(coord_collect, momentum_collect, weighting):
    """

    create set of postion points and momentum points

    """

    points_coords = create_point_array(coord_collect, weighting)
    points_momentum = create_point_array(momentum_collect, weighting)

    points = {}
    points['position'] = points_coords
    points['momentum'] = points_momentum

    return points


def read_group_values(group):
    """

    convert values from position and momentum datasets into points
    group -- base group of points from hdf file

    """

    hdf_datasets = ParticlesFunctor()
    group.visititems(hdf_datasets)
    weighting = hdf_datasets.weighting
    position_values = DatasetReader('position')
    momentum_values = DatasetReader('momentum')

    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)
    points = create_points_library(position_values, momentum_values, weighting)
    return points


def read_mass(group):
    hdf_mass = MassFunctor()
    group.visititems(hdf_mass)
    return hdf_mass.mass


def read_points_group(group):
    """

    convert values from position and momentum datasets into points
    group -- base group of points from hdf file

    """

    hdf_datasets = ParticlesFunctor()
    group.visititems(hdf_datasets)
    weighting = hdf_datasets.weighting
    bound_electrons = hdf_datasets.bound_electrons
    position_values = DatasetReader('position')
    momentum_values = DatasetReader('momentum')

    if len(hdf_datasets.positions) == 0 or len(hdf_datasets.momentum) == 0:
        return [], [], [], [], []
    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)

    points = []

    if len(bound_electrons) == 0:
        points = create_points_array(position_values, momentum_values)
    else:
        points = create_points_array_bound_electrons(position_values, momentum_values, bound_electrons)

    dimention_position = position_values.get_dimension()
    unit_SI_position = position_values.get_unit_si_array()
    unit_SI_momentum = momentum_values.get_unit_si_array()

    dimention_momentum = momentum_values.get_dimension()
    dimensions = Dimensions(dimention_position, dimention_momentum)

    return points, weighting, dimensions, unit_SI_position, unit_SI_momentum


def read_position_offset(group):

    offset_values = []

    hdf_datasets = ParticlesFunctor()
    group.visititems(hdf_datasets)
    position_offset_values = DatasetReader('positionOffset')
    position_offset_group = hdf_datasets.position_offset[0]
    position_offset_group.visititems(position_offset_values)

    if position_offset_values.get_dimension() == 2:
        offset_values = [list(x) for x in
                         zip(position_offset_values.vector_x, position_offset_values.vector_y)]

    elif position_offset_values.get_dimension() == 3:
        offset_values = [list(x) for x in
                         zip(position_offset_values.vector_x, position_offset_values.vector_y,
                             position_offset_values.vector_z)]

    offset_unit_si = position_offset_values.get_unit_si_array()

    return offset_values, offset_unit_si


def write_group_values(hdf_file_reduction, group, library_datasets, offset):
    """

    write values from point library to hdf file
    hdf_file_reduction -- result file
    group -- base group of partilces from original file
    result -- library points

    """

    hdf_datasets = ParticlesFunctor()
    group.visititems(hdf_datasets)
    position_values = DatasetReader('position')
    momentum_values = DatasetReader('momentum')
    position_offset_values = DatasetReader('positionOffset')
    position_offset_group = hdf_datasets.position_offset[0]
    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)
    position_offset_group.visititems(position_offset_values)

    writen_position = points_writer(hdf_file_reduction, library_datasets, 'position')
    writen_momentum = points_writer(hdf_file_reduction, library_datasets, 'momentum')

    writen_weighting = dataset_writer(hdf_file_reduction, library_datasets)
    position_group.visititems(writen_position)
    momentum_group.visititems(writen_momentum)
    group.visititems(writen_weighting)


def create_datasets_from_vector(reduced_data, dimensions, position_offset):

    library_datasets = {}

    size_values = len(reduced_data[0]) - 1

    position_x = reduced_data[:, 0]
    position_y = reduced_data[:, 1]
    position_z = []
    if dimensions.dimension_position == 3:
        position_z = reduced_data[:, 2]

    momentum_x = reduced_data[:, dimensions.dimension_position]
    momentum_y = reduced_data[:, dimensions.dimension_position + 1]
    momentum_z = []
    if dimensions.dimension_momentum == 3:
        momentum_z = reduced_data[:, dimensions.dimension_position + 2]

    bound_electrons = reduced_data[:, size_values]

    position_offset_x = position_offset[:, 0]
    position_offset_y = position_offset[:, 1]
    position_offset_z = []
    if dimensions.dimension_position == 3:
        position_offset_z = position_offset[:, 2]

    library_datasets['position/x'] = position_x
    library_datasets['position/y'] = position_y
    library_datasets['position/z'] = position_z

    library_datasets['positionOffset/x'] = position_offset_x
    library_datasets['positionOffset/y'] = position_offset_y
    library_datasets['positionOffset/z'] = position_offset_z

    library_datasets['momentum/x'] = momentum_x
    library_datasets['momentum/y'] = momentum_y
    library_datasets['momentum/z'] = momentum_z

    library_datasets['boundElectrons'] = bound_electrons

    return library_datasets


def write_group_values(hdf_file_reduction, group, reduced_data, weights):
    """

    write values from point library to hdf file
    hdf_file_reduction -- result file
    group -- base group of partilces from original file
    result -- library points

    """

    hdf_datasets = ParticlesFunctor()
    group.visititems(hdf_datasets)
    position_values = DatasetReader('position')
    momentum_values = DatasetReader('momentum')
    position_offset_values = DatasetReader('positionOffset')
    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_offset_group = hdf_datasets.position_offset[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)
    position_offset_group.visititems(position_offset_values)

    write_position = points_writer(hdf_file_reduction, reduced_data, 'position')
    write_momentum = points_writer(hdf_file_reduction, reduced_data, 'momentum')
    write_position_offset = points_writer(hdf_file_reduction, reduced_data, 'positionOffset')
    write_weighting = dataset_writer(hdf_file_reduction, weights, 'weighting')
    write_bound_electrons = dataset_writer(hdf_file_reduction, reduced_data['boundElectrons'], 'boundElectrons')

    position_group.visititems(write_position)
    momentum_group.visititems(write_momentum)
    position_offset_group.visititems(write_position_offset)
    group.visititems(write_weighting)
    group.visititems(write_bound_electrons)


def write_patch_group(group, hdf_file_reduction, num_particles_offset, num_particles):
    patch_group = ReadPatchGroup()
    group.visititems(patch_group)
    patch_writter = PatchValuesWriter(hdf_file_reduction, num_particles_offset, num_particles)
    patch_group.patch_group[0].visititems(patch_writter)


def read_patches_values(group):

    patch_group = ReadPatchGroup()
    group.visititems(patch_group)
    patch_values = ReadPatchValues()
    patch_group.patch_group[0].visititems(patch_values)
    return patch_values.numParticles, patch_values.numParticlesOffset

def get_absolute_coordinates(data, position_offset, unit_si_offset, unit_si_position, dimensions, unit_si_momentum):

    absolute_result = []

    unit_si_position = numpy.array(unit_si_position)
    unit_si_offset = numpy.array(unit_si_offset)
    unit_si_momentum = numpy.array(unit_si_momentum)

    i = 0
    for point in data:
        offset = position_offset[i]
        coordinates = numpy.array(point[0:dimensions.dimension_position])

        absolute_coordinates = coordinates * unit_si_position + offset * unit_si_offset
        momentum = point[dimensions.dimension_position:dimensions.dimension_momentum + dimensions.dimension_position]
        absolute_momentum = momentum * unit_si_momentum
        other_values = point[dimensions.dimension_momentum + dimensions.dimension_position: len(point)]

        absolute_point = numpy.append(absolute_coordinates, absolute_momentum)
        absolute_point = numpy.append(absolute_point, other_values)

        absolute_result.append(absolute_point.tolist())
        i+=1

    return absolute_result


def get_relative_coordinates(absolute_coordinates, unit_si_offset,
                             unit_si_position, dimensions, unit_si_momentum):

    relative_result = []
    offset = []

    unit_si_position = numpy.array(unit_si_position)
    unit_si_offset = numpy.array(unit_si_offset)
    unit_si_momentum = numpy.array(unit_si_momentum)

    for point in absolute_coordinates:
        coordinates = numpy.array(point[0:dimensions.dimension_position])
        position_offset = numpy.divide(coordinates, unit_si_position)
        position_offset = position_offset.astype(int)

        offset.append(position_offset.tolist())

        relative_coordinates = numpy.divide((coordinates - position_offset * unit_si_offset), unit_si_position)

        momentum = point[dimensions.dimension_position:dimensions.dimension_momentum + dimensions.dimension_position]

        relative_momentum = numpy.divide(momentum, unit_si_momentum)

        other_values = point[dimensions.dimension_momentum + dimensions.dimension_position:
                             len(point)]

        relative_point = numpy.append(relative_coordinates, relative_momentum)
        relative_point = numpy.append(relative_point, other_values)

        relative_result.append(relative_point.tolist())

    relative_result = numpy.array(relative_result)
    offset = numpy.array(offset)

    return relative_result, offset