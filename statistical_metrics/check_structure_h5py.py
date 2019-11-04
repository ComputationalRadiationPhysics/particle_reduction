import sys
sys.path.append("../")
import numpy
import read_hdf_file
import argparse
import h5py


def get_patches_ranges(hdf_file, group, position_collection):

    num_particles, num_particles_offset = read_hdf_file.read_patches_values(group)
    ranges_patches = num_particles_offset
    size = hdf_file[position_collection.vector[0]][()].size
    ranges_patches = numpy.append(ranges_patches, int(size))
    ranges_patches.astype(int)

    return ranges_patches

def get_particles_sizes(hdf_file, group):

    hdf_datasets = read_hdf_file.ParticlesFunctor()
    group.visititems(hdf_datasets)

    position_collection, momentum_collection, weighting, bound_electrons = read_hdf_file.read_points_group_v2(hdf_datasets)

    if position_collection.get_dimension() == 0 or momentum_collection.get_dimension() == 0:
        return

    size = hdf_file[position_collection.vector[0]][()].size
    print("current size: "+ str(size))


    position_offset, unit_si_offset = read_hdf_file.read_position_offset(hdf_datasets)

    ranges_patches = get_patches_ranges(hdf_file, group, position_collection)

    print("ranges_patches  ")
    print(ranges_patches)
    dimension_position = position_collection.get_dimension()
    dimension_momentum = momentum_collection.get_dimension()

    print("dimension_position "+ str(dimension_position))
    print("dimension_momentum "+ str(dimension_momentum))


def base_reading_function(hdf_file_name):
    hdf_file = h5py.File(hdf_file_name, 'a')

    particles_name = read_hdf_file.get_particles_name(hdf_file)
    particles_collect = read_hdf_file.ParticlesGroups(particles_name)
    print(particles_collect)
    hdf_file.visititems(particles_collect)

    for group in particles_collect.particles_groups:
        print('name group ' + str(group.name))

        get_particles_sizes(hdf_file, group)


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="main reduction")

    parser.add_argument("-hdf_file", metavar='hdf_file', type=str,
                        help="hdf file")

    args = parser.parse_args()

    base_reading_function(args.hdf_file)

