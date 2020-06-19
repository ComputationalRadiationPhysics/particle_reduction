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
    def factory(type, parameters, mass):

        if type == "random": return Random_thinning_algorithm.Random_thinning_algorithm(parameters.reduction_percent)
        if type == "number_conservative": return Number_conservative_thinning_algorithm.Number_conservative_thinning_algorithm(parameters.reduction_percent)
        if type == "energy_conservative":
            return Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm(parameters.reduction_percent, mass)
        if type == "kmeans":
            divisions = [16, 16, 4]
            return k_means_clustering_algorithm.K_means_clustering_algorithm(parameters.reduction_percent, parameters.max_iterations, parameters.tolerance,
                                                                             divisions)
        if type == "kmeans_avg":
            divisions = [32, 60, 16]
            return k_means_merge_average_algorithm.K_means_merge_average_algorithm(parameters.reduction_percent, divisions,
                                                                                   parameters.max_iterations, parameters.tolerance)
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
            dset = Dataset(component_values.dtype, extent=component_values.shape)
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
            value = numpy.concatenate((absolute_coordinates[i], absolute_momentum[i]), axis=None)
        else:
            value = numpy.concatenate((absolute_coordinates[i], absolute_momentum[i], bound_elctrons[i]), axis=None)

        result_data.append(value)

    return result_data


def get_relative_coordinates(data, dict_data_indexes, unit_si_offset,
                             unit_si_position):

    particle_index_range = dict_data_indexes["position"]

    offset = []
    relative_coordinates = []

    for i in range(0, len(data)):
        coordinates = data[i][particle_index_range[0]:particle_index_range[1]]
        position_offset = numpy.divide(coordinates, unit_si_offset)
        position_offset = position_offset.astype(int)
        offset.append(position_offset)

        relative_point = numpy.divide((coordinates - position_offset * unit_si_offset), unit_si_position)
        relative_coordinates.append(relative_point)

    offset = numpy.array(offset)
    relative_coordinates = numpy.array(relative_coordinates)

    return relative_coordinates, offset



def get_relative_data(data, partcles_spices, dict_data_indexes, weights):


    relative_data = []
    for record_name in dict_data_indexes:
        if record_name == "position":
            continue
        indexes_data = dict_data_indexes[record_name]
        current_record = partcles_spices[record_name]
        values = data[:, indexes_data[0]:indexes_data[1]]
        unmacroweighted = get_unmacroweighted(values, current_record, weights)
        unit_si = get_unit_SI(current_record)
        revative_values = get_relative_values(unmacroweighted, unit_si)
        relative_data.append(revative_values)

    return relative_data



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


def make_particle_patches_structure(particle_species, particle_species_reduction):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    dtype_patches = particle_species.particle_patches["numParticles"][SCALAR].dtype
    extent = particle_species.particle_patches["numParticles"][SCALAR].shape[0]


    dset = Dataset(dtype_patches, extent=[extent])

    particle_species_reduction.particle_patches["numParticles"][SCALAR].reset_dataset(dset)
    particle_species_reduction.particle_patches["numParticlesOffset"][SCALAR]. \
        reset_dataset(dset)

    patches_structure = particle_species.particle_patches
    patches_structure_reduction = particle_species_reduction.particle_patches
    make_vector_structures(patches_structure, patches_structure_reduction, "offset", -1)
    make_vector_structures(patches_structure, patches_structure_reduction, "extent", -1)


def make_copy_vector_structures(particle_species, particle_species_reduction, name_of_dataset, name_of_copy_dataset):

    base_record = particle_species[name_of_dataset]
    record_redutcion = particle_species_reduction[name_of_copy_dataset]

    for vector in base_record:

        struction_size = base_record[vector].shape[0]
        dtype = base_record[vector].dtype
        d = Dataset(dtype, [struction_size])
        current_vector = record_redutcion[vector]
        current_vector.reset_dataset(d)


def create_copy_dataset_structures(particle_species, particle_species_reduction):

    for record_component_name, record_component in particle_species.items():
        copy_record_name = record_component_name + "_copy"
        make_copy_vector_structures(particle_species, particle_species_reduction, record_component_name, copy_record_name)



def copy_unit_dimension(obj, reduction_obj):

    unit_dimension = obj.unit_dimension

    dict_units = [Unit_Dimension.L, Unit_Dimension.M, Unit_Dimension.T, Unit_Dimension.I,
                  Unit_Dimension.theta, Unit_Dimension.N, Unit_Dimension.J]

    result_unit_dimesion = {}
    i = 0
    for dict_value in dict_units:
        result_unit_dimesion[dict_value] = unit_dimension[i]
        i = i + 1

    reduction_obj.set_unit_dimension(result_unit_dimesion)


def copy_momentum_parameters(current_momentum, reduction_momentum):

    for attr in current_momentum.attributes:
        if attr != "unitDimension" and attr != "timeOffset":
            reduction_momentum.set_attribute(attr, current_momentum.get_attribute(attr))

    macro_weighted = current_momentum.time_offset
    reduction_momentum.set_time_offset(macro_weighted)

    copy_unit_dimension(current_momentum, reduction_momentum)


def create_dataset_structures(particle_species, particle_species_reduction, reduction_size):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    weights = particle_species["weighting"][SCALAR]

    d = Dataset(weights.dtype, [reduction_size])
    weights_reduction = particle_species_reduction["weighting"][SCALAR]
    weights_reduction.reset_dataset(d)
    make_vector_structures(particle_species, particle_species_reduction, "position", reduction_size)
    make_vector_structures(particle_species, particle_species_reduction, "momentum", reduction_size)
    make_vector_structures(particle_species, particle_species_reduction, "positionOffset", reduction_size)

    if is_vector_exist("numParticles", particle_species):
        make_particle_patches_structure(particle_species, particle_species_reduction)

    copy_momentum_parameters(particle_species["position"], particle_species_reduction["position"])
    copy_momentum_parameters(particle_species["momentum"], particle_species_reduction["momentum"])
    copy_momentum_parameters(particle_species["positionOffset"], particle_species_reduction["positionOffset"])

    if is_vector_exist("boundElectrons", particle_species):

        bound_electrons = particle_species["boundElectrons"][SCALAR]
        d = Dataset(bound_electrons.dtype, [reduction_size])
        weights_reduction = particle_species_reduction["boundElectrons"][SCALAR]
        weights_reduction.reset_dataset(d)


def count_reduction_size(reduction_percent, num_particles):

    result_size = 0
    num_particles_reduction = [0]

    for value in num_particles:
        patch_size = int(round(value * (1 -reduction_percent)))
        result_size += patch_size
        num_particles_reduction.append(patch_size)
    num_particles_reduction = numpy.cumsum(num_particles_reduction)

    return result_size, num_particles_reduction


def write_reduction_vector(particle_species, name_vector, series_hdf_reduction, reduction_vector, pos_in_vector,
                           previos_idx, current_idx):

    position_offset = particle_species[name_vector]
    for vector in position_offset:
        current_offset = reduction_vector[:, pos_in_vector].astype(numpy.int32)
        position_offset[vector][previos_idx:current_idx] = current_offset
        series_hdf_reduction.flush()
        idx_position_offset_vector = idx_position_offset_vector + 1


def write_record_copy(particle_record, values_to_write, series_hdf_reduction, previos_idx, current_idx):

    pos_vector_in_reduction_data = 0


    for vector in particle_record:
        current_type = particle_record[vector].dtype
        current_reduced_data = values_to_write[:, pos_vector_in_reduction_data].astype(current_type)
        particle_record[vector][previos_idx:current_idx] = current_reduced_data
        series_hdf_reduction.flush()
        pos_vector_in_reduction_data = pos_vector_in_reduction_data + 1

def write_draft_position(particle_species, series_hdf_reduction, position_values,
                     offset, previos_idx, current_idx):

    series_hdf_reduction.flush()
    position = particle_species["position_copy"]
    position_offset = particle_species["positionOffset_copy"]

    write_record_copy(position, position_values, series_hdf_reduction, previos_idx, current_idx)
    write_record_copy(position_offset, offset, series_hdf_reduction, previos_idx, current_idx)


def write_draft_copy(data, dict_data_indexes, series_hdf_reduction, base_record,
                     previos_idx, current_idx):
    print(len(data))

    for record_name in dict_data_indexes:
        if record_name == "position":
            continue

        particle_index_range = dict_data_indexes[record_name]
        values = data[particle_index_range[0]:particle_index_range[1]][0]
        draft_record_name = record_name +"_copy"
        current_record = base_record[draft_record_name]
        write_record_copy(current_record, values, series_hdf_reduction, previos_idx, current_idx)


def write_patches_information(particle_species, num_particles, num_particles_offset):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    for i in range(0, len(num_particles)):
        particle_species.particle_patches["numParticles"][SCALAR].store(i, numpy.array([num_particles[i]], dtype=numpy.ulonglong))

    for i in range(0, len(num_particles_offset)):
        particle_species.particle_patches["numParticlesOffset"][SCALAR].store(i, numpy.array([num_particles_offset[i]],
                                                                                       dtype=numpy.ulonglong))

def copy_main_version(series_hdf_reduction, particle_species_reduction, result_size):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    series_hdf_reduction.flush()
    copy_values = particle_species_reduction["weighting_copy"][SCALAR][0:result_size]
    series_hdf_reduction.flush()
    particle_species_reduction["weighting"][SCALAR][0:result_size] = copy_values
    series_hdf_reduction.flush()
    del particle_species_reduction["weighting_copy"]

    position = particle_species_reduction["position"]
    momentum = particle_species_reduction["momentum"]
    position_offset = particle_species_reduction["positionOffset"]
    position_draft = particle_species_reduction["position_copy"]
    momentum_draft = particle_species_reduction["momentum_copy"]
    position_offset_draft = particle_species_reduction["positionOffset_copy"]

    for vector in position_offset_draft:
        copy_values = position_offset_draft[vector][0:result_size]
        series_hdf_reduction.flush()
        position_offset[vector][0:result_size] = copy_values
        series_hdf_reduction.flush()

    del particle_species_reduction["positionOffset_copy"]

    for vector in position_draft:
        copy_values = position_draft[vector][0:result_size]
        series_hdf_reduction.flush()
        position[vector][0:result_size] = copy_values
        series_hdf_reduction.flush()

    del particle_species_reduction["position_copy"]

    for vector in momentum_draft:
        copy_values = momentum_draft[vector][0:result_size]
        series_hdf_reduction.flush()
        momentum[vector][0:result_size] = copy_values
        series_hdf_reduction.flush()

    del particle_species_reduction["momentum_copy"]

    if is_vector_exist("boundElectrons_copy", particle_species_reduction):
        bound_electrons_draft = particle_species_reduction["boundElectrons_copy"]
        bound_electrons = particle_species_reduction["boundElectrons"]
        test_values = bound_electrons_draft[SCALAR][0:result_size]
        series_hdf_reduction.flush()
        bound_electrons[SCALAR][0:result_size] = test_values
        series_hdf_reduction.flush()
        del particle_species_reduction["boundElectrons_copy"]


def set_scalar_group(particle_species, particle_species_reduction, name_group):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    group = particle_species[name_group]
    reduction_group = particle_species_reduction[name_group]
    mass_value = particle_species[name_group][SCALAR][0]
    mass_reduction = particle_species_reduction[name_group][SCALAR]
    mass_reduction.make_constant(mass_value)
    reduction_group.set_time_offset(group.time_offset)

    for attr in group.attributes:
        if attr != "unitDimension" and attr != "timeOffset":
            reduction_group.set_attribute(attr, group.get_attribute(attr))

    copy_unit_dimension(group, reduction_group)


def write_scalar_groups(particle_species, particle_species_reduction):

    if is_vector_exist("mass", particle_species):
        set_scalar_group(particle_species, particle_species_reduction, "mass")

    if is_vector_exist("charge", particle_species):
        set_scalar_group(particle_species, particle_species_reduction, "charge")


def get_chunk_sizes(particle_species, series_hdf):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    dataset_size_max = particle_species["position"]["x"].shape[0]
    max_chunk_size = 1e7

    ranges_patches = []

    if is_vector_exist("numParticlesOffset", particle_species):
        ranges_patches = particle_species.particle_patches["numParticlesOffset"][SCALAR].load()
        series_hdf.flush()
        ranges_patches = numpy.append(ranges_patches, dataset_size_max)
    else:
        num_full_chunks = dataset_size_max/max_chunk_size
        if num_full_chunks > 1:
            ranges_patches = numpy.full(dataset_size_max, num_full_chunks)
            ranges_patches = numpy.append(ranges_patches, dataset_size_max - max_chunk_size * num_full_chunks)
        else:
            ranges_patches = numpy.append(ranges_patches, 0)
            ranges_patches = numpy.append(ranges_patches, dataset_size_max)

    return ranges_patches

def get_dimentions(position, momentum):

    dimension_momentum = 0
    for obj in momentum.items():
        dimension_momentum = dimension_momentum + 1

    dimension_position = 0
    for obj in position.items():
        dimension_position = dimension_position + 1

    dimensions = Dimensions(dimension_position, dimension_momentum)

    return dimensions


def process_patches_in_group_v2(particle_species, series_hdf, series_hdf_reduction,
                                particle_species_reduction, algorithm):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    create_copy_dataset_structures(particle_species, particle_species_reduction)

    copy_attributes(particle_species, particle_species_reduction)


    position_offset = particle_species["positionOffset"]
    position = particle_species["position"]

    weights = particle_species["weighting"]

    weighted_position = get_macroweighted(series_hdf, position, weights, idx_start, idx_end)

    weighted_position_offset = get_macroweighted(series_hdf, position_offset, weights, idx_start, idx_end)
    unit_SI_position = get_unit_SI(position)
    unit_SI_position_offset = get_unit_SI(position_offset)
    absolute_coordinates = get_absolute_coordinates(weighted_position,
                                                    weighted_position_offset, unit_SI_position, unit_SI_position_offset)

    absolute_coordinates = numpy.array(absolute_coordinates)
    return absolute_coordinates, unit_SI_position, unit_SI_position_offset


def get_data(series_hdf, particle_species, weights, idx_start, idx_end):

    dict_data_indexes = {}
    idx_start_component = 0
    idx_end_component = 0
    data = []
    for record_component_name, record_component in particle_species.items():


        if record_component_name == "position" or record_component_name == "weighting" \
                or record_component_name == "positionOffset":
            continue

        weighted_values = get_macroweighted(series_hdf, record_component, weights, idx_start, idx_end)

        unit_SI = get_unit_SI(record_component)
        absolute_values = get_absolute_values(weighted_values, unit_SI)

        if len(absolute_values[0]) == 1:
            data.append(absolute_values[:, 0])
        else:
            for i in range(0, len(absolute_values[0])):
                data.append(absolute_values[:, i])

        for sub_record_component_name, sub_record_component in record_component.items():
            idx_end_component = idx_end_component + 1

        dict_data_indexes[record_component_name] = [idx_start_component, idx_end_component]
        idx_start_component = idx_end_component

    return data, dict_data_indexes, idx_end_component


def process_patches_in_group_v2(particle_species, series_hdf, series_hdf_reduction,
                                particle_species_reduction, algorithm):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    create_copy_dataset_structures(particle_species, particle_species_reduction)

    copy_attributes(particle_species, particle_species_reduction)


    ranges_patches = get_chunk_sizes(particle_species, series_hdf)

    previos_idx = 0
    current_idx = 0

    result_size = 0

    new_num_particles = []


    weights = particle_species["weighting"]

    for i in range(0, len(ranges_patches) - 1):

        idx_start = int(ranges_patches[i])
        idx_end = int(ranges_patches[i + 1])

        absolute_coordinates, unit_si_position, unit_si_offset = get_coordinates(series_hdf, particle_species, idx_start, idx_end)

        data, dict_data_indexes, last_idx = get_data(series_hdf, particle_species, weights, idx_start, idx_end)

        dict_data_indexes["position"] = [last_idx, last_idx + len(absolute_coordinates[0])]


        for i in range(0, len(absolute_coordinates[0])):
            data.append(absolute_coordinates[:, i])

        data = numpy.transpose(data)
        weights = weights[SCALAR][idx_start: idx_end]
        series_hdf.flush()


        reduced_data, reduced_weight = algorithm._run(data, weights, dict_data_indexes)

        relative_coordinates, offset = get_relative_coordinates(reduced_data, dict_data_indexes, unit_si_offset,
                                 unit_si_position)

        relative_data = get_relative_data(reduced_data, particle_species, dict_data_indexes, reduced_weight)


        new_num_particles.append(len(reduced_weight))
        current_idx += len(reduced_weight)



        write_draft_position(particle_species_reduction, series_hdf_reduction, relative_coordinates,
                             offset, previos_idx, current_idx)
        write_draft_copy(relative_data, dict_data_indexes, series_hdf_reduction, particle_species_reduction,
                         previos_idx, current_idx)

        previos_idx += len(reduced_weight)
        result_size = previos_idx

    new_num_particles_offset = numpy.cumsum(new_num_particles[0:len(new_num_particles) - 1], dtype=int)
    new_num_particles_offset = numpy.insert(new_num_particles_offset, 0, 0)

    create_dataset_structures(particle_species, particle_species_reduction, result_size)

    write_patches_information(particle_species_reduction, new_num_particles, new_num_particles_offset)
   # write_scalar_groups(particle_species, particle_species_reduction)

   # copy_main_version(series_hdf_reduction, particle_species_reduction, result_size)


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

        result_num_particles.append(len(reduced_data_patch))

    return reduced_data, reduced_weights, result_num_particles


if __name__ == "__main__":
    """ Parse arguments from command line     
    
    """

    parser = argparse.ArgumentParser(description="main reduction")

    parser.add_argument("-algorithm", metavar='algorithm', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file to be reducted")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="result reduction hdf file")

    parser.add_argument("-reduction_percent", metavar='reduction_percent', type=float,
                        help="part of the particles to reduce( used in Energy_conservative, random, number_conservative,"
                             "kmeans, kmeans_avg algorithms and in voronoi probalistic algorithm )")

    parser.add_argument("-momentum_tol", metavar='tolerance_momentum', type=float,
                        help="tolerance of momentum ( in SI ), used in voronoi algorithm")

    parser.add_argument("-position_lol", metavar='tolerance_position', type=float,
                        help="tolerance of position( in SI ), used in voronoi algorithm")


    parser.add_argument("-leveling_coefficient", metavar='leveling_coefficient', type=float,
                        help="leveling_coefficient, used in leveling algorithm")


    args = parser.parse_args()

    if args.algorithm == 'voronoi':
        tolerance = [args.momentum_tol, args.position_lol]
        parameters = Voronoi_algorithm.VoronoiMergingAlgorithmParameters(tolerance)
        voronoi_algorithm(args.hdf, args.hdf_re, args.momentum_tol, args.position_lol)

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

    elif args.algorithm == 'leveling':
        parameters = Leveling_thinning_algorithm.Leveling_thinning_algorithm_parameters(args.leveling_coefficient)
        base_reduction_function(args.hdf, args.hdf_re, "leveling", parameters)


