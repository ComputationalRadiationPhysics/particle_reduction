import h5py
import re
import Voronoi_algorithm


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

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):
            if node.name.endswith('weighting'):
                self.weighting = node.value

        if isinstance(node, h5py.Group):
            if node.name.endswith('position'):
                self.positions.append(node)

            if node.name.endswith('momentum'):
                self.momentum.append(node)

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


class CoverterVoronoiToPoints():
    def __init__(self, result_points):

        self.points = {}
        self.coords =[]
        self.momentum = []
        for point in result_points:
            if point.points['position'] != None:
                vector_coords = point.points['position'][0].coords
                weight = point.points['position'][0].weight
                self.coords.append(Voronoi_algorithm.Point(vector_coords, weight))

            if point.points['momentum'] != None:
                vector_coords = point.points['momentum'][0].coords
                weight = point.points['momentum'][0].weight
                self.momentum.append(Voronoi_algorithm.Point(vector_coords, weight))


        self.points['position'] = self.coords
        self.points['momentum'] = self.momentum

    def get_points(self):
        return self.points



class DatasetWriter_Voronoi_cells():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, result_points, name_dataset):

        self.vector_x = []
        self.vector_y = []
        self.vector_z = []
        self.weighting = []
        self.hdf_file = hdf_file
        self.name_dataset = name_dataset
        self.demention = len(result_points[0].points[self.name_dataset][0].coords)

        for point in result_points:

            if point.points[self.name_dataset] != None:

                vector_coords = point.points[self.name_dataset][0].coords
                if self.demention == 2:
                    self.vector_x.append(vector_coords[0])
                    self.vector_y.append(vector_coords[1])
                if self.demention == 3:
                    self.vector_x.append(vector_coords[0])
                    self.vector_y.append(vector_coords[1])
                    self.vector_y.append(vector_coords[2])

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):

            dataset_x = self.name_dataset + '/x'

            dataset_y = self.name_dataset + '/y'
            dataset_z = self.name_dataset + '/z'

            if node.name.endswith(dataset_x):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_x)
            elif node.name.endswith(dataset_y):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_y)

            elif node.name.endswith(dataset_z):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_z)

            elif node.name.endswith('weighting'):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.weighting)

        return None


class DatasetWriter():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, result_points, name_dataset):
        self.dataset_x = name_dataset + '/x'
        self.dataset_y = name_dataset + '/y'
        self.dataset_z = name_dataset + '/z'

        self.vector_x = result_points[self.dataset_x]
        self.vector_y = result_points[self.dataset_y]
        self.vector_z = result_points[self.dataset_z]
        self.hdf_file = hdf_file

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):

            if node.name.endswith(self.dataset_x):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_x)
            elif node.name.endswith(self.dataset_y):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_y)

            elif node.name.endswith(self.dataset_z):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.vector_z)

        return None


class WeightWriter():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, weight):

        self.weighting = weight
        self.hdf_file = hdf_file

    def __call__(self, name, node):

        if isinstance(node, h5py.Dataset):

            if node.name.endswith('weighting'):
                node_name = node.name
                del self.hdf_file[node.name]
                dset = self.hdf_file.create_dataset(node_name, data=self.weighting, dtype=float)
        return None


class PatchValuesWriter():
    """

    Write dataset into result hdf file
    name_dataset -- name recorded dataset
    hdf_file -- result hdf file
    result_points -- points to write to hdf file

    """

    def __init__(self, hdf_file, numParticles, numParticlesOffset):

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


def create_point_array(coord_collection, weighting):
    """

    create array of 2-d, 3-d points from datasets
    coord_collection -- datasets from hdf file

    """

    point_array = []

    dimension = coord_collection.get_dimension()

    if dimension == 3:
        for i in range(0, len(coord_collection.vector_x)):
            point_array.append(Voronoi_algorithm.Point([coord_collection.vector_x[i], coord_collection.vector_y[i],
                                      coord_collection.vector_z[i]], weighting[i]))

    elif dimension == 2:
        for i in range(0, len(coord_collection.vector_x)):
            point_array.append(Voronoi_algorithm.Point([coord_collection.vector_x[i], coord_collection.vector_y[i]], weighting[i]))

    return point_array


def create_points_array_ver2(coord_collection, momentum_collection):
    """

    create array of 2-d, 3-d points from datasets
    coord_collection -- datasets from hdf file

    """

    vector_coords = []

    dimension_coord = coord_collection.get_dimension()

    dimension_momentum = momentum_collection.get_dimension()

    if dimension_coord == 3 and dimension_momentum == 3:
        for i in range(0, len(coord_collection.vector_x)):
            vector_coords.append([coord_collection.vector_x[i], coord_collection.vector_y[i],
                                  coord_collection.vector_z[i], momentum_collection.vector_x[i],
                                  momentum_collection.vector_y[i], momentum_collection.vector_z[i]])

    elif dimension_coord == 3 and dimension_momentum == 2:
        for i in range(0, len(coord_collection.vector_x)):
            vector_coords.append([coord_collection.vector_x[i], coord_collection.vector_y[i],
                                  coord_collection.vector_z[i], momentum_collection.vector_x[i],
                                  momentum_collection.vector_y[i]])

    elif dimension_coord == 2 and dimension_momentum == 3:
        for i in range(0, len(coord_collection.vector_x)):
            vector_coords.append([coord_collection.vector_x[i], coord_collection.vector_y[i],
                                  momentum_collection.vector_x[i], momentum_collection.vector_y[i],
                                  momentum_collection.vector_z[i]])

    elif dimension_coord == 2 and dimension_momentum == 2:
        for i in range(0, len(coord_collection.vector_x)):
            vector_coords.append([coord_collection.vector_x[i], coord_collection.vector_y[i],
                                  momentum_collection.vector_x[i], momentum_collection.vector_y[i]])

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


def write_group_values(hdf_file_reduction, group, library_datasets, num_particles_offset=None, num_particles=None):
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
    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)
    writen_position = DatasetWriter(hdf_file_reduction, library_datasets, 'position')
    writen_momentum = DatasetWriter(hdf_file_reduction, library_datasets, 'momentum')
    writen_weighting = WeightWriter(hdf_file_reduction, library_datasets)
    position_group.visititems(writen_position)
    momentum_group.visititems(writen_momentum)
    group.visititems(writen_weighting)


def create_datasets_from_vector(reduced_data, dimensions):

    library_datasets = {}

    position_x = []
    position_y = []
    position_z = []

    momentum_x = []
    momentum_y = []
    momentum_z = []

    for point in reduced_data:
        if dimensions.dimension_position == 3:
            position_x.append(point[0])
            position_y.append(point[1])
            position_z.append(point[2])

        if dimensions.dimension_position == 2:
            position_x.append(point[0])
            position_z.append(point[1])

        if dimensions.dimension_momentum == 3:
            momentum_x.append(point[dimensions.dimension_position])
            momentum_y.append(point[dimensions.dimension_position + 1])
            momentum_z.append(point[dimensions.dimension_position + 2])

        if dimensions.dimension_momentum == 2:
            momentum_x.append(point[0])
            momentum_y.append(point[dimensions.dimension_position + 1])

        library_datasets['position/x'] = position_x
        library_datasets['position/y'] = position_y
        library_datasets['position/z'] = position_z

        library_datasets['momentum/x'] = momentum_x
        library_datasets['momentum/y'] = momentum_y
        library_datasets['momentum/z'] = momentum_z

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
    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_group.visititems(position_values)
    momentum_group.visititems(momentum_values)

    writen_position = DatasetWriter(hdf_file_reduction, reduced_data,  'position')
    writen_momentum = DatasetWriter(hdf_file_reduction, reduced_data, 'momentum')
    writen_weighting = WeightWriter(hdf_file_reduction, weights)
    position_group.visititems(writen_position)
    momentum_group.visititems(writen_momentum)
    group.visititems(writen_weighting)


        group.visititems(patch_group)
        patch_writter = PatchValuesWriter(hdf_file_reduction, num_particles_offset, num_particles)
        patch_group.patch_group[0].visititems(patch_writter)


def read_patches_values(group):

    patch_group = ReadPatchGroup()
    group.visititems(patch_group)
    patch_values = ReadPatchValues()
    patch_group.patch_group[0].visititems(patch_values)
    return patch_values.numParticles, patch_values.numParticlesOffset