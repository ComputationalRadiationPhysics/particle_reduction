import read_hdf_file
import h5py


def count_weight_coordinates_difference(first_group, second_group):

    hdf_datasets = read_hdf_file.Particles_functor()
    first_group.visititems(hdf_datasets)
    position_group_first = hdf_datasets.positions[0]
    position_values_first = read_hdf_file.Dataset_reader('position')
    position_group_first.visititems(position_values_first)
    size_of_positions_first = len(position_values_first.vector_x)
    position_values_first.get_demention()
    hdf_datasets = read_hdf_file.Particles_functor()
    second_group.visititems(hdf_datasets)
    position_group_second = hdf_datasets.positions[0]
    position_values_second = read_hdf_file.Dataset_reader('position')
    position_group_second.visititems(position_values_second)


    size_of_positions_second = len(position_values_second.vector_x)



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


