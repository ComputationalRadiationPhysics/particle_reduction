import read_hdf_file
import h5py


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

