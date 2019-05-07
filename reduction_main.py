from shutil import copyfile
import numpy
import read_hdf_file
import os
import copy
import h5py
import argparse
import Algorithms.Random_thinning_algorithm as Random_thinning_algorithm
import Algorithms.Number_conservative_thinning_algorithm as Number_conservative_thinning_algorithm
import Algorithms.Energy_conservative_thinning_algorithm as Energy_conservative_thinning_algorithm
import Algorithms.k_means_clustering_algorithm as k_means_clustering_algorithm
import Algorithms.k_means_merge_average_algorithm as k_means_merge_average_algorithm
import Algorithms.Voronoi_algorithm as Voronoi_algorithm
import Algorithms.Leveling_thinning_algorithm as Leveling_thinning_algorithm

class Algorithm:
    # Create based on class name:
    def factory( type, parameters, mass):

        if type == "random": return Random_thinning_algorithm.Random_thinning_algorithm(parameters.reduction_percent)
        if type == "number_conservative": return Number_conservative_thinning_algorithm.Number_conservative_thinning_algorithm(parameters.reduction_percent)
        if type == "energy_conservative":
            return Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm(parameters.reduction_percent, mass)
        if type == "kmeans":
            divisions = [10, 10]
            return k_means_clustering_algorithm.K_means_clustering_algorithm(parameters.reduction_percent, parameters.max_iterations, parameters.tolerance,
                                                                             divisions)
        if type == "kmeans-avg":
            divisions = [10, 10]
            return k_means_merge_average_algorithm.K_means_merge_average_algorithm(parameters.reduction_percent, parameters.max_iterations, parameters.tolerance,
                                                                                   divisions)
        if type == "voronoi":
            return Voronoi_algorithm.VoronoiMergingAlgorithm(parameters.tolerance)
        if type == "leveling":
            return Leveling_thinning_algorithm.Leveling_thinning_algorithm(parameters.leveling_coefficient)

        if type == "vranic":pass
        assert 0, "Bad type_algoritm: " + type
    factory = staticmethod(factory)


def base_reduction_function(hdf_file_name, hdf_file_reduction_name, type, parameters):

    particles_collect, hdf_file_reduction = get_particles_groups(hdf_file_name, hdf_file_reduction_name)

    for group in particles_collect.particles_groups:
        print('name group '+ str(group.name))
        process_reduction_group(type, group, hdf_file_reduction, parameters)


def process_reduction_group(type, group, hdf_file_reduction, parameters):
    mass = read_hdf_file.read_mass(group)
    algorithm = Algorithm.factory(type, parameters, mass)
    process_patches_in_group(hdf_file_reduction, group, algorithm)


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

def get_relative_coordinates(absolute_coordinates, unit_SI_offset,
                             unit_SI_position, dimensions, unit_SI_momentum):

    relative_coordinates = []
    offset = []

    for i in range(0, len(absolute_coordinates)):

        relative_point = []
        offset_point = []
        for j in range(0, dimensions.dimension_position):
            current_position_offset = int(absolute_coordinates[i][j]/unit_SI_offset[j])

            offset_point.append(current_position_offset)
            current_relative_coordinate = (absolute_coordinates[i][j]-current_position_offset * unit_SI_offset[j])\
                                          /unit_SI_position[j]

            relative_point.append(current_relative_coordinate)

        for j in range(dimensions.dimension_position, len(absolute_coordinates[0])):
            current_relative_momentum = absolute_coordinates[i][j]/unit_SI_momentum[int(j - dimensions.dimension_position)]
            relative_point.append(current_relative_momentum)

        relative_coordinates.append(relative_point)
        offset.append(offset_point)

    return relative_coordinates, offset


def process_patches_in_group(hdf_file_reduction, group, algorithm):

    data, weights, dimensions, unit_si_position, unit_si_momentum \
        = read_hdf_file.read_points_group(group)

    if len(data) == 0:
        return

    position_offset, unit_si_offset = read_hdf_file.read_position_offset(group)
    num_particles_offset, num_particles_offset = read_hdf_file.read_patches_values(group)
    absolute_coordinates = get_absolute_coordinates(data, position_offset, unit_si_offset, unit_si_position, dimensions,
                                                    unit_si_momentum)

    algorithm.dimensions = dimensions

    reduced_data, reduced_weights, result_num_particles = \
        iterate_patches(absolute_coordinates, weights, num_particles_offset, algorithm)

    relative_coordinates, offset = get_relative_coordinates(reduced_data, unit_si_offset,
                                                            unit_si_position, dimensions, unit_si_momentum)

    result_num_particles_offset = numpy.cumsum(result_num_particles[0:len(result_num_particles) - 1], dtype=int)
    result_num_particles_offset = numpy.insert(result_num_particles_offset, 0, 0)

    read_hdf_file.write_patch_group(group, hdf_file_reduction, result_num_particles_offset, result_num_particles)

    library_datasets = read_hdf_file.create_datasets_from_vector(relative_coordinates, dimensions, offset)

    read_hdf_file.write_group_values(hdf_file_reduction, group, library_datasets, reduced_weights)


def get_particles_groups(hdf_file_name, hdf_file_reduction_name):
    copyfile(hdf_file_name, hdf_file_reduction_name)
    hdf_file = h5py.File(hdf_file_name, 'a')
    hdf_file_reduction = h5py.File(hdf_file_reduction_name, 'a')
    particles_name = read_hdf_file.get_particles_name(hdf_file_reduction)
    particles_collect = read_hdf_file.ParticlesGroups(particles_name)
    hdf_file.visititems(particles_collect)

    return particles_collect, hdf_file_reduction


def random_thinning_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):

    parameters = Random_thinning_algorithm.Random_thinning_algorithm_parameters(reduction_percent)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "random", parameters)


def number_conservative_thinning_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):
    parameters = Number_conservative_thinning_algorithm.\
        Number_conservative_thinning_algorithm_parameters(reduction_percent)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "number_conservative", parameters)


def leveling_thinning_algorithm(hdf_file_name, hdf_file_reduction_name, leveling_coefficient):
    algorithm = Leveling_thinning_algorithm.Leveling_thinning_algorithm(leveling_coefficient)
    thinning_base_procedure(hdf_file_name, hdf_file_reduction_name, algorithm)


def energy_conservative_thinning_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):
    parameters = Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm_parameters(
        reduction_percent)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "energy_conservative", parameters)


def k_means_cluster_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):
    parameters = k_means_clustering_algorithm.K_means_clustering_algorithm_parameters(reduction_percent)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "kmeans", parameters)


def k_means_avg_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):
    parameters = k_means_merge_average_algorithm.K_means_merge_average_algorithm_parameters(reduction_percent)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "kmeans_avg", parameters)


def iterate_patches(data, weights, num_particles_offset, algorithm, bound_electrons):
    print('iterate_patches ')

    ranges_patches = num_particles_offset

    ranges_patches = numpy.append(ranges_patches, len(data))
    ranges_patches.astype(int)


    reduced_data = []
    reduced_weights = []
    result_num_particles = []


    for i in range(0, len(ranges_patches) - 1):
        start = int(ranges_patches[i])
        end = int(ranges_patches[i + 1])
        print('start '+ str(start))
        print('end ' + str(end))
        copy_data = copy.deepcopy(data[start:end])

        copy_weights = copy.deepcopy(weights[start:end])
        reduced_data_patch, reduced_weight_patch = algorithm._run(copy_data, copy_weights)


        for point in reduced_data_patch:
            reduced_data.append(point)

        for weight in reduced_weight_patch:
            reduced_weights.append(weight)
        result_num_particles.append(len(reduced_data_patch))

    return reduced_data, reduced_weights, result_num_particles


if __name__ == "__main__":
    """ Parse arguments from command line """

    parser = argparse.ArgumentParser(description="voronoi reduction")

    parser.add_argument("-algorithm", metavar='algorithm', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="reducted hdf file")

    parser.add_argument("-reduction_percent", metavar='reduction_percent', type=float,
                        help="part of the particles to reduce")

    parser.add_argument("-sample_amount", metavar='sample_amount', type=int,
                        help="amount of sample ")

    parser.add_argument("-momentum_tol", metavar='tolerance_momentum', type=float,
                        help="tolerance of momentum")

    parser.add_argument("-momentum_pos", metavar='tolerance_position', type=float,
                        help="tolerance of position")

    parser.add_argument("-particles_type", metavar='particles_type', type=str,
                        help="types of particles")

    parser.add_argument("-leveling_coefficient", metavar='leveling_coefficient', type=float,
                        help="leveling_coefficient")

    parser.add_argument("-k_means_subdivision", metavar='leveling_coefficient', type=list,
                        help="leveling_coefficient")

    args = parser.parse_args()

    if args.algorithm == 'voronoi':
        tolerance = [args.momentum_tol, args.momentum_pos]
        parameters = Voronoi_algorithm.VoronoiMergingAlgorithmParameters(tolerance)
        base_reduction_function(args.hdf, args.hdf_re, "voronoi", parameters)

    elif args.algorithm == 'random':
        parameters = Random_thinning_algorithm.Random_thinning_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "random", parameters)

    elif args.algorithm == 'number_conservative':
        parameters = Number_conservative_thinning_algorithm.Number_conservative_thinning_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "number_conservative", parameters)

    elif args.algorithm == 'energy_conservative':
        parameters = Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "energy_conservative", parameters)

    elif args.algorithm == 'kmeans':
        parameters = k_means_clustering_algorithm.K_means_clustering_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "kmeans", parameters)

    elif args.algorithm == 'kmeans_avg':
        parameters = k_means_merge_average_algorithm.K_means_merge_average_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "kmeans_avg", parameters)

    elif args.algorithm == 'vranic_algorithm':
        Vranic_algorithm_algorithm(args.hdf, args.hdf_re, args.momentum_tol, args.particles_type)

    elif args.algorithm == 'leveling':
        parameters = Leveling_thinning_algorithm.Leveling_thinning_algorithm_parameters(args.leveling_coefficient)
        base_reduction_function(args.hdf, args.hdf_re, "leveling", parameters)


