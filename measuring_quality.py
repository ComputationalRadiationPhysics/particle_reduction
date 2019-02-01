import read_hdf_file
import h5py


class WeightReader():

    def __init__(self):
        self.weight = []

    def __call__(self, name, node):
        if isinstance(node, h5py.Dataset):
            if node.name.endswith('weighting'):
                self.weight = node.value


def convert_mass_to_array(mass, size_of_positions):

    dimension = values.get_dimension()

def compute_weight_positions_sum(mass, position_values):

    demention = position_values.get_demention()

    sum_x = 0.
    sum_y = 0.
    sum_z = 0.

    if demention == 3:
        size_of_position = len(position_values.vector_x)
        for i in range(0, size_of_position):
            sum_x += position_values.vector_x[i] * mass[i]
            sum_y += position_values.vector_y[i] * mass[i]
            sum_z += position_values.vector_z[i] * mass[i]

    elif demention == 2:
        size_of_position = len(position_values.vector_x)
        for i in range(0, size_of_position):
            sum_x += position_values.vector_x[i] * mass[i]
            sum_y += position_values.vector_y[i] * mass[i]

    return sum_x, sum_y, sum_z


def count_weight_coordinates_difference(first_group, second_group):

    hdf_datasets = read_hdf_file.Particles_functor()
    first_group.visititems(hdf_datasets)
    position_group_first = hdf_datasets.positions[0]
    position_values_first = read_hdf_file.Dataset_reader('position')
    position_group_first.visititems(position_values_first)
    mass_reader_first = Mass_reader(first_group)
    first_group.visititems(mass_reader_first)

    size_of_positions_first = len(position_values_first.vector_x)
    position_values_first.get_demention()
    first_mass = convert_mass_to_array(mass_reader_first.mass, size_of_positions_first)
    sum_first = compute_weight_positions_sum(first_mass, position_values_first)

    hdf_datasets = read_hdf_file.Particles_functor()
    second_group.visititems(hdf_datasets)
    position_group_second = hdf_datasets.positions[0]
    position_values_second = read_hdf_file.Dataset_reader('position')
    position_group_second.visititems(position_values_second)
    mass_reader_second = Mass_reader(second_group)
    second_group.visititems(mass_reader_second)
    size_of_positions_second = len(position_values_second.vector_x)
    second_mass = convert_mass_to_array(mass_reader_second.mass, size_of_positions_second)

    sum_second = compute_weight_positions_sum(second_mass, position_values_second)

def count_weight_momentum_difference(first_group, second_group):

    print('weight momentum')

    hdf_datasets = read_hdf_file.Particles_functor()
    first_group.visititems(hdf_datasets)
    momentum_group_first = hdf_datasets.momentum[0]
    momentum_values_first = read_hdf_file.Dataset_reader('momentum')
    momentum_group_first.visititems(momentum_values_first)
    mass_reader_first = Mass_reader(first_group)
    first_group.visititems(mass_reader_first)

    size_of_positions_first = len(momentum_values_first.vector_x)
    momentum_values_first.get_demention()
    first_mass = convert_mass_to_array(mass_reader_first.mass, size_of_positions_first)
    sum_first = compute_weight_positions_sum(first_mass, momentum_values_first)

    print('first sum == ' + str(sum_first))

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
    particles_groups_first = read_hdf_file.Particles_groups(particles_name_first)
    first_hdf_file.visititems(particles_groups_first)

    particles_name_second = read_hdf_file.get_particles_name(second_hdf_file)
    particles_groups_second = read_hdf_file.Particles_groups(particles_name_second)
    second_hdf_file.visititems(particles_groups_second)

    size_groups = len(particles_groups_second.particles_groups)


    for i in range(0, size_groups):
        count_weight_coordinates_difference(particles_groups_first.particles_groups[i],
                                            particles_groups_second.particles_groups[i])
        count_weight_momentum_difference(particles_groups_first.particles_groups[i],
                                            particles_groups_second.particles_groups[i])


