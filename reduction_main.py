from shutil import copyfile
import read_hdf_file
import argparse
import os
import h5py

def voronoi_reduction(hdf_file, hdf_file_reduction, tolerance_momentum, tolerance_position):

    name_hdf_file_reduction = ''

    if hdf_file != '':
        if os.path.exists(hdf_file):
            name = hdf_file[:-4]
            idx_of_name = name.rfind('/')
            if idx_of_name != -1:
                name_hdf_file_reduction = hdf_file_reduction + hdf_file[idx_of_name + 1: -6] + 'reduction.h5'
            else:
                name_hdf_file_reduction = hdf_file_reduction + hdf_file[:-3] + '.h5'

            tolerances = [tolerance_momentum, tolerance_position]
            voronoi_algorithm(hdf_file, hdf_file_reduction, tolerances)
        else:
            print('The .hdf file does not exist')


def voronoi_algorithm(hdf_file_name, hdf_file_reduction_name, tolerances):

    copyfile(hdf_file_name, hdf_file_reduction_name)

    hdf_file = h5py.File(hdf_file_name, 'a')
    hdf_file_reduction = h5py.File(hdf_file_reduction_name, 'a')
    particles_name = read_hdf_file.get_particles_name(hdf_file)

    particles_collect = read_hdf_file.Particles_groups(particles_name)
    hdf_file.visititems(particles_collect)
    for group in particles_collect.particles_groups:
        read_hdf_file.particle_group_iteration(group, hdf_file_reduction, tolerances)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="voronoi reduction")

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="reducted hdf file")

    parser.add_argument("-momentum_tol", metavar='tolerance_momentum', type=float,
                        help="tolerance of momentum")

    parser.add_argument("-momentum_pos", metavar='tolerance_position', type=float,
                        help="tolerance of position")

    args = parser.parse_args()
    voronoi_reduction(args.hdf, args.hdf_re, args.momentum_tol, args.momentum_pos)


