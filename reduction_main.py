from shutil import copyfile
import numpy
import read_hdf_file
import os
import copy
import h5py
import argparse
import time
import math
import Algorithms.Random_thinning_algorithm as Random_thinning_algorithm
import Algorithms.Number_conservative_thinning_algorithm as Number_conservative_thinning_algorithm
import Algorithms.Energy_conservative_thinning_algorithm as Energy_conservative_thinning_algorithm
import Algorithms.k_means_clustering_algorithm as k_means_clustering_algorithm
import Algorithms.k_means_merge_average_algorithm as k_means_merge_average_algorithm
import Algorithms.Voronoi_algorithm as Voronoi_algorithm
import Algorithms.Leveling_thinning_algorithm as Leveling_thinning_algorithm
import Algorithms.Voronoi_probabilistic_algorithm as Voronoi_probabilistic_algorithm

import openpmd_api

from openpmd_api import Series, Access_Type, Dataset, Mesh_Record_Component, \
    Unit_Dimension

class Dimensions:
    def __init__(self, dimension_position, dimension_momentum):
        self.dimension_position = dimension_position
        self.dimension_momentum = dimension_momentum


class Algorithm:
    # Create based on class name:
    def factory( type, parameters, mass):

        if type == "random": return Random_thinning_algorithm.Random_thinning_algorithm(parameters.reduction_percent)
        if type == "number_conservative": return Number_conservative_thinning_algorithm.Number_conservative_thinning_algorithm(parameters.reduction_percent)
        if type == "energy_conservative":
            return Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm(parameters.reduction_percent, mass)
        if type == "kmeans":
            divisions = [16, 40]
            return k_means_clustering_algorithm.K_means_clustering_algorithm(parameters.reduction_percent, parameters.max_iterations, parameters.tolerance,
                                                                             divisions)
        if type == "kmeans_avg":
            divisions = [32, 60]
            return k_means_merge_average_algorithm.K_means_merge_average_algorithm(parameters.reduction_percent, parameters.max_iterations, parameters.tolerance,
                                                                                   divisions)
        if type == "voronoi":
            return Voronoi_algorithm.VoronoiMergingAlgorithm(parameters.tolerance)
        if type == "voronoi_prob":
            voronoi_parameters = Voronoi_probabilistic_algorithm.Voronoi_probabilistic_algorithm_parameters\
                (parameters.reduction_percent, parameters.ratio_left_particles)
            return Voronoi_probabilistic_algorithm.Voronoi_probabilistic_algorithm(voronoi_parameters)
        if type == "leveling":
            return Leveling_thinning_algorithm.Leveling_thinning_algorithm(parameters.leveling_coefficient)

        if type == "vranic":pass
        assert 0, "Bad type_algoritm: " + type
    factory = staticmethod(factory)


def copy_all_root_attributes(series_hdf, series_hdf_reduction):

    series_hdf_reduction.set_author(series_hdf.author)

    series_hdf_reduction.set_date(series_hdf.date)
    series_hdf_reduction.set_iteration_encoding(series_hdf.iteration_encoding)
    series_hdf_reduction.set_iteration_format(series_hdf.iteration_format)
    series_hdf_reduction.set_meshes_path(series_hdf.meshes_path)

    series_hdf_reduction.set_openPMD(series_hdf.openPMD)
    series_hdf_reduction.set_openPMD_extension(series_hdf.openPMD_extension)
    series_hdf_reduction.set_particles_path(series_hdf.particles_path)
    series_hdf_reduction.set_software(series_hdf.software)
    series_hdf_reduction.set_software_version(series_hdf.software_version)


def copy_attributes(start_obj, end_obj):

    for attr in start_obj.attributes:
        end_obj.set_attribute(attr, start_obj.get_attribute(attr))


def copy_iteration_parameters(current_iteration, reduction_iteration):

    time_unit_SI = current_iteration.time_unit_SI()
    time = current_iteration.time()
    dt = current_iteration.dt()

    copy_attributes(current_iteration, reduction_iteration)
    reduction_iteration.set_time(time) \
        .set_dt(dt) \
        .set_time_unit_SI(time_unit_SI)


def copy_meshes(series_hdf, reduction_series, current_iteration, reduction_iteration):

    mesh_record_start = current_iteration.meshes
    mesh_record_end = reduction_iteration.meshes

    copy_attributes(mesh_record_start, mesh_record_end)

    for mesh in current_iteration.meshes:
        current_mesh = current_iteration.meshes[mesh]
        reduction_mesh = reduction_iteration.meshes[mesh]
        copy_attributes(current_mesh, reduction_mesh)
        for component in current_mesh:
            mesh_component = current_mesh[component]
            component_values = current_mesh[component][()]
            reduction_values = reduction_mesh[component]
            series_hdf.flush()
            dset = Dataset(component_values.dtype, extent=[len(component_values), len(component_values[0])])

            copy_attributes(mesh_component, reduction_values)
            reduction_values.reset_dataset(dset)

            reduction_mesh[component][()] = component_values
            reduction_series.flush()


def base_reduction_function(hdf_file_name, hdf_file_reduction_name, type_algorithm, parameters):

    series_hdf = openpmd_api.Series(hdf_file_name, openpmd_api.Access_Type.read_only)
    series_hdf_reduction = openpmd_api.Series(hdf_file_reduction_name, openpmd_api.Access_Type.create)

    copy_all_root_attributes(series_hdf, series_hdf_reduction)
    algorithm = Algorithm.factory(type_algorithm, parameters, 1.)
    for iteration in series_hdf.iterations:
        current_iteration = series_hdf.iterations[iteration]
        reduction_iteration = series_hdf_reduction.iterations[iteration]
        copy_iteration_parameters(current_iteration, reduction_iteration)
        copy_meshes(series_hdf, series_hdf_reduction, current_iteration, reduction_iteration)
        process_iteration_group(algorithm, current_iteration, series_hdf, series_hdf_reduction, reduction_iteration)


def base_reduction_voronoi(hdf_file_name, hdf_file_reduction_name, type, parameters):
    particles_collect, hdf_file_reduction = get_particles_groups(hdf_file_name, hdf_file_reduction_name)

    for group in particles_collect.particles_groups:
        process_iteration_group(type, group, hdf_file_reduction, parameters)


def check_item_exist(particle_species, name_item):

    item_exist = False
    for value in particle_species.items():
        if value[0] == name_item:
            item_exist = True

    return item_exist


def process_iteration_group(algorithm, iteration, series_hdf, series_hdf_reduction, reduction_iteration):

    for name_group in iteration.particles:

        if not (check_item_exist(iteration.particles[name_group], "momentum") and
                check_item_exist(iteration.particles[name_group], "position")):
            continue
        process_patches_in_group_v2(iteration.particles[name_group], series_hdf,
                                    series_hdf_reduction, reduction_iteration.particles[name_group], algorithm)


def process_patches_in_group(hdf_file_reduction, group, algorithm):

    data, weights, dimensions, unit_si_position, unit_si_momentum \
        = read_hdf_file.read_points_group(group)

    if len(data) == 0:
        return

    position_offset, unit_si_offset = read_hdf_file.read_position_offset(group)
    num_particles_offset, num_particles_offset = read_hdf_file.read_patches_values(group)

    absolute_coordinates = read_hdf_file.get_absolute_coordinates(data, position_offset, unit_si_offset, unit_si_position, dimensions,
                                                                  unit_si_momentum)

    algorithm.dimensions = dimensions

    reduced_data, reduced_weights, result_num_particles = \
        iterate_patches(absolute_coordinates, weights, num_particles_offset, algorithm)

    relative_coordinates, offset = read_hdf_file.get_relative_coordinates(reduced_data, unit_si_offset,
                                                            unit_si_position, dimensions, unit_si_momentum)

    result_num_particles_offset = numpy.cumsum(result_num_particles[0:len(result_num_particles) - 1], dtype=int)
    result_num_particles_offset = numpy.insert(result_num_particles_offset, 0, 0)

    read_hdf_file.write_patch_group(group, hdf_file_reduction, result_num_particles_offset, result_num_particles)

    library_datasets = read_hdf_file.create_datasets_from_vector(relative_coordinates, dimensions, offset)
    read_hdf_file.write_group_values(hdf_file_reduction, group, library_datasets, reduced_weights)


def get_dimensions(position_values, momentum_values):

    dimension_position = position_values.get_dimension()
    dimension_momentum = momentum_values.get_dimension()
    dimensions = Dimensions(dimension_position, dimension_momentum)

    return dimensions


def writing_reduced_information(rewrite_start_node, hdf_file_reduction, relative_coordinates,
                                hdf_datasets, group, relative_position_offset, dimensions, bound_electrons, weights):

    position_group = hdf_datasets.positions[0]
    momentum_group = hdf_datasets.momentum[0]
    position_offset_group = hdf_datasets.position_offset[0]

    write_position = read_hdf_file.vector_writer(hdf_file_reduction, relative_coordinates, 'position', 0)

    write_position.is_first_part = rewrite_start_node
    position_group.visititems(write_position)

    write_momentum = read_hdf_file.vector_writer(hdf_file_reduction, relative_coordinates, 'momentum',
                                                 dimensions.dimension_position)

    write_momentum.is_first_part = rewrite_start_node
    momentum_group.visititems(write_momentum)

    write_position_offset = read_hdf_file.vector_writer(hdf_file_reduction, relative_position_offset, 'positionOffset',
                                                        0)

    write_position_offset.is_first_part = rewrite_start_node
    position_offset_group.visititems(write_position_offset)

    write_weighting = read_hdf_file.dataset_writer(hdf_file_reduction, weights, 'weighting')

    write_weighting.is_first_part = rewrite_start_node
    group.visititems(write_weighting)

    if bound_electrons != "":
        idx_bound_electrons = len(relative_coordinates[0]) - 1
        write_bound_electrons = read_hdf_file.dataset_writer(hdf_file_reduction, relative_coordinates[:, idx_bound_electrons], 'boundElectrons')

        write_bound_electrons.is_first_part = rewrite_start_node
        group.visititems(write_bound_electrons)


def absolute_item(values, offset, unit_si_offset, unit_si_position):

    absolute_result = []

    i = 0

    for point in values:

        absolute_coord = point * unit_si_position + offset[i] * unit_si_offset
        absolute_result.append(absolute_coord)
        i = +1

    return absolute_result


def get_absolute_coordinates(series, position, position_offset, idx_start, idx_end):

    absolute_coordinates = []

    for value in position.items():
        name_value = value[0]
        position_axis = position[name_value]
        position_offset_axis = position_offset[str(name_value)]
        position_dataset = position_axis[idx_start:idx_end]
        position_offset_dataset = position_offset_axis[idx_start:idx_end]
        series.flush()
        item_values = absolute_item(position_dataset, position_offset_dataset, position_offset_axis.unit_SI, position_axis.unit_SI)
        absolute_coordinates.append(item_values)

    absolute_coordinates = numpy.transpose(absolute_coordinates)

    return absolute_coordinates


def absolute_momentum_array(array_dataset, unit_si_momentum):

    absolute_momentum = []
    for point in array_dataset:
        absolute_momentum.append(point * unit_si_momentum)

    return absolute_momentum


def get_absolute_momentum(series, momentum_values, idx_start, idx_end):

    absolute_momentum = []
    for value in momentum_values.items():
        name_value = value[0]
        momentum_axis = momentum_values[name_value]

        array_dataset = momentum_axis[idx_start:idx_end]
        series.flush()
        unit_si_momentum = momentum_axis.unit_SI
        item_absolute_values = absolute_momentum_array(array_dataset, unit_si_momentum)

        absolute_momentum.append(item_absolute_values)

    absolute_momentum = numpy.transpose(absolute_momentum)

    return absolute_momentum


def create_input_data(absolute_coordinates, absolute_momentum, bound_elctrons):

    result_data = []

    for i in range(0, len(absolute_coordinates)):
        value = []
        if len(bound_elctrons) == 0:
            value = numpy.concatenate((absolute_coordinates[0], absolute_momentum[0]), axis=None)
        else:
            value = numpy.concatenate((absolute_coordinates[i], absolute_momentum[i], bound_elctrons[i]), axis=None)

        result_data.append(value)

    return result_data


def get_relative_coordinates(absolute_values, unit_si_offset,
                             unit_si_position, unit_si_momentum):

    relative_result = []
    offset = []

    unit_si_position = numpy.array(unit_si_position)
    dimension_position = len(unit_si_position)
    unit_si_offset = numpy.array(unit_si_offset)
    unit_si_momentum = numpy.array(unit_si_momentum)
    dimension_momentum = len(unit_si_momentum)

    for point in absolute_values:
        coordinates = numpy.array(point[0:dimension_position])
        position_offset = numpy.divide(coordinates, unit_si_position)

        position_offset = position_offset.astype(int)

        offset.append(position_offset.tolist())

        relative_coordinates = numpy.divide((coordinates - position_offset * unit_si_offset), unit_si_position)

        momentum = point[dimension_position:dimension_momentum + dimension_position]

        relative_momentum = numpy.divide(momentum, unit_si_momentum)

        other_values = point[dimension_momentum + dimension_position:
                             len(point)]

        relative_point = numpy.append(relative_coordinates, relative_momentum)
        relative_point = numpy.append(relative_point, other_values)

        relative_result.append(relative_point.tolist())

    relative_result = numpy.array(relative_result)
    offset = numpy.array(offset)

    return relative_result, offset


def get_units_si(position, position_offset, momentum):

    position_si = []
    position_offset_si = []
    for value in position.items():
        name_value = value[0]
        position_axis = position[name_value]
        position_offset_axis = position_offset[str(name_value)]
        position_si.append(position_axis.unit_SI)
        position_offset_si.append(position_offset_axis.unit_SI)

    momentum_si = []
    for value in momentum.items():
        name_value = value[0]
        momentum_axis = momentum[name_value]
        momentum_si.append(momentum_axis.unit_SI)

    return position_si, position_offset_si, momentum_si


def is_vector_exist(vector_name, particle_species):

    for key in particle_species:
        if key == vector_name:
            return True

    return False


def make_vector_structures(particle_species, particle_species_reduction, name_of_structure, struction_size):

    position = particle_species[name_of_structure]
    position_redutcion = particle_species_reduction[name_of_structure]

    for vector in position:

        if struction_size == -1:
            struction_size = position[vector].shape[0]
        dtype = position[vector].dtype
        d = Dataset(dtype, [struction_size])
        unit_si = position[vector].unit_SI

        current_vector = position_redutcion[vector]
        current_vector.reset_dataset(d)
        current_vector.set_unit_SI(unit_si)

    for i in range(0, len(ranges_patches) - 1):

        idx_start = int(ranges_patches[i])
        idx_end = int(ranges_patches[i + 1])

        unit_si_position = position_collection.get_unit_si_array()
        unit_si_momentum = momentum_collection.get_unit_si_array()

        data = read_hdf_file.create_points_array(hdf_file, idx_start, idx_end, position_collection, momentum_collection, bound_electrons)

        if len(data) == 0:
            return

        weights = read_hdf_file.get_weights(hdf_file, idx_start, idx_end, hdf_datasets.weighting)

        position_offset_vector = read_hdf_file.get_position_offset(hdf_file, idx_start, idx_end, position_offset)
        absolute_coordinates = read_hdf_file.get_absolute_coordinates(data, position_offset_vector, unit_si_offset,
                                                                      unit_si_position, dimensions, unit_si_momentum)
        copy_data = copy.deepcopy(absolute_coordinates)

        copy_weights = copy.deepcopy(weights)
        reduced_data, reduced_weight = algorithm._run(copy_data, copy_weights)

        relative_coordinates, relative_position_offset = read_hdf_file.get_relative_coordinates(reduced_data, unit_si_offset,
                                                                              unit_si_position, dimensions,
                                                                              unit_si_momentum)

        writing_reduced_information(rewrite_start_node, hdf_file_reduction, relative_coordinates,
                                    hdf_datasets, group, relative_position_offset, dimensions, bound_electrons, reduced_weight)
        rewrite_start_node = False



def get_particles_groups(hdf_file, hdf_file_reduction_name):

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


def voronoi_prob_algorithm(hdf_file_name, hdf_file_reduction_name, reduction_percent):
    ratio_left_particles = 20
    voronoi_parameters = Voronoi_probabilistic_algorithm.Voronoi_probabilistic_algorithm_parameters \
        (reduction_percent, ratio_left_particles)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "voronoi_prob", voronoi_parameters)

    
def voronoi_algorithm(hdf_file_name, hdf_file_reduction_name, momentum_tolerance, position_tolerance):

    tolerance = [momentum_tolerance, position_tolerance]
    parameters = Voronoi_algorithm.VoronoiMergingAlgorithmParameters(tolerance)
    base_reduction_function(hdf_file_name, hdf_file_reduction_name, "voronoi", parameters)


def iterate_patches(data, weights, num_particles_offset, algorithm):

    ranges_patches = num_particles_offset

    ranges_patches = numpy.append(ranges_patches, len(data))
    ranges_patches.astype(int)


    reduced_data = []
    reduced_weights = []
    result_num_particles = []


    for i in range(0, len(ranges_patches) - 1):
        start = int(ranges_patches[i])
        end = int(ranges_patches[i + 1])

        start_time = time.time()

        copy_data = copy.deepcopy(data[start:end])
        copy_weights = copy.deepcopy(weights[start:end])
        reduced_data_patch, reduced_weight_patch = algorithm._run(copy_data, copy_weights)


        for point in reduced_data_patch:
            reduced_data.append(point)

        for weight in reduced_weight_patch:
            reduced_weights.append(weight)

        end_time = time.time()
        print('TIME ' + str(end_time - start_time))
        result_num_particles.append(len(reduced_data_patch))

    return reduced_data, reduced_weights, result_num_particles


if __name__ == "__main__":
    """ Parse arguments from command line """

    parser = argparse.ArgumentParser(description="main reduction")

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

    parser.add_argument("-divide_particles", metavar='size_of_divide_particles', type=float,
                        help="size_of_divide_particles")

    args = parser.parse_args()

    if args.algorithm == 'voronoi':
        tolerance = [args.momentum_tol, args.momentum_pos]
        parameters = Voronoi_algorithm.VoronoiMergingAlgorithmParameters(tolerance)
        base_reduction_function(args.hdf, args.hdf_re, "voronoi", parameters)

    elif args.algorithm == 'voronoi_prob':
        parameters = Voronoi_probabilistic_algorithm.Voronoi_probabilistic_algorithm_parameters(args.reduction_percent, args.divide_particles)
        base_reduction_function(args.hdf, args.hdf_re, "voronoi_prob", parameters)

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


