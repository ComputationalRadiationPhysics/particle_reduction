from shutil import copyfile
import numpy
import argparse
import math
from argparse import RawTextHelpFormatter
import textwrap
import Algorithms.Random_thinning_algorithm as Random_thinning_algorithm
import Algorithms.Number_conservative_thinning_algorithm as Number_conservative_thinning_algorithm
import Algorithms.Energy_conservative_thinning_algorithm as Energy_conservative_thinning_algorithm
import Algorithms.k_means_clustering_algorithm as k_means_clustering_algorithm
import Algorithms.k_means_merge_average_algorithm as k_means_merge_average_algorithm
import Algorithms.Voronoi_algorithm as Voronoi_algorithm
import Algorithms.Leveling_thinning_algorithm as Leveling_thinning_algorithm
import Algorithms.Voronoi_probabilistic_algorithm as Voronoi_probabilistic_algorithm
import Algorithms.Vranic_algorithm as Vranic_algorithm

import openpmd_api

from openpmd_api import Series, Access_Type, Dataset, Mesh_Record_Component,\
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
            return Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm(parameters.reduction_percent)
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

        if type == "vranic":
            return Vranic_algorithm.Vranic_merging_algorithm(parameters)
        assert 0, "Bad type_algoritm: " + type
    factory = staticmethod(factory)


def copy_all_root_attributes(series_hdf, series_hdf_reduction):
    series_hdf_reduction.set_iteration_encoding(series_hdf.iteration_encoding)
    series_hdf_reduction.set_date(series_hdf.date)

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


def transform_into_weighted_values(series, record_component, weights, idx_start, idx_end):

    weightingPower = record_component.get_attribute("weightingPower")
    weighted_values = []
    for name, values in record_component.items():
        current_values = values[idx_start: idx_end]
        series.flush()
        components = multiply_macroweighted(current_values, weights, weightingPower)
        weighted_values.append(components)

    weighted_values = numpy.transpose(weighted_values)

    return weighted_values

def transform_from_unweighted_values(values, record_component, weights):

    weightingPower = record_component.get_attribute("weightingPower")
    weighted_values = []
    for i in range(0, len(values)):
        weighted_value = values[i] /( weights[i] ** weightingPower)
        weighted_values.append(weighted_value)

   # weighted_values = numpy.transpose(weighted_values)

    return weighted_values

def get_non_transformed_values(series, record_component, idx_start, idx_end):

    weighted_values = []
    for name, values in record_component.items():
        current_values = values[idx_start: idx_end]
        series.flush()
        weighted_values.append(current_values)

    weighted_values = numpy.transpose(weighted_values)
    return weighted_values

def multiply_macroweighted(values, weights, weightingPower):

    multiply_values = []

    for i in range(0, len(values)):
        multiply_values.append(values[i] * weights[i]**weightingPower)

    return multiply_values


def get_macroweighted(series, record_component, weights, idx_start, idx_end):

    macroWeighted = record_component.get_attribute("macroWeighted")

    weighted_values = []

    if macroWeighted == 1:
        weighted_values = transform_into_weighted_values(series, record_component, weights, idx_start, idx_end)
    else:
        weighted_values = get_non_transformed_values(series, record_component, idx_start, idx_end)

    return weighted_values

def get_unmacroweighted(values, record_component, weights):

    macroWeighted = record_component.get_attribute("macroWeighted")

    weighted_values = []

    if macroWeighted == 1:
        weighted_values = transform_from_unweighted_values(values, record_component, weights)
    else:
        weighted_values = values

    return weighted_values


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


def get_dimensions(position_values, momentum_values):

    dimension_position = position_values.get_dimension()
    dimension_momentum = momentum_values.get_dimension()
    dimensions = Dimensions(dimension_position, dimension_momentum)

    return dimensions

def get_absolute_coordinates(position, position_offset, unit_si_offset, unit_si_position):

    absolute_coordinates = []

    size_position = len(position)

    for i in range(0, size_position):
        point = position[i]
        offset = position_offset[i]

        absolute_point = point * unit_si_position + offset * unit_si_offset
        absolute_coordinates.append(absolute_point)

    return absolute_coordinates



def get_absolute_values(record_values, unit_si):

    absolute_values = []
    size_values = len(record_values)

    for i in range(0, size_values):
        point = record_values[i]
        absolute_point = point * unit_si
        absolute_values.append(absolute_point)

    return numpy.array(absolute_values)

def get_relative_values(values, unit_si):

    relative_values = []
    size_values = len(values)

    for i in range(0, size_values):
        point = values[i]
        absolute_point = point / unit_si
        relative_values.append(absolute_point)

    return numpy.array(relative_values)

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

        relative_values = get_relative_values(unmacroweighted, unit_si)
        for i in range(0, len(relative_values[0])):
            relative_data.append(relative_values[:, i])

    relative_data = numpy.array(relative_data)
    relative_data = relative_data.T

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

def make_vector_structures(record, reducted_record, record_size):

    for vector in record:

        if record_size == -1:
            record_size = record[vector].shape[0]
        dtype = record[vector].dtype
        d = Dataset(dtype, [record_size])
        unit_si = record[vector].unit_SI

        current_vector = reducted_record[vector]
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
    make_vector_structures(patches_structure["offset"], patches_structure_reduction["offset"], -1)
    make_vector_structures(patches_structure["extent"], patches_structure_reduction["extent"], -1)


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


def copy_record_attributes(base_record, reduction_record):

    for attr in base_record.attributes:
        if attr != "unitDimension" and attr != "timeOffset" and attr != "shape":
            reduction_record.set_attribute(attr, base_record.get_attribute(attr))

    macro_weighted = base_record.time_offset
    reduction_record.set_time_offset(macro_weighted)

    copy_unit_dimension(base_record, reduction_record)

def create_dataset_structures(particle_species, particle_species_reduction, reduction_size):

    for record_name, record in particle_species.items():
        reducted_record = particle_species_reduction[record_name]
        make_vector_structures(record, reducted_record, reduction_size)
        copy_record_attributes(record, reducted_record)

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


def write_draft_copy(data, reduced_weight, dict_data_indexes, series_hdf_reduction, base_record,
                     previos_idx, current_idx):

    for record_name in dict_data_indexes:
        if record_name == "position":
            continue

        particle_index_range = dict_data_indexes[record_name]
        values = data[:, particle_index_range[0]:particle_index_range[1]]
        draft_record_name = record_name +"_copy"
        current_record = base_record[draft_record_name]
        write_record_copy(current_record, values, series_hdf_reduction, previos_idx, current_idx)


    weighting_record = base_record["weighting_copy"][Mesh_Record_Component.SCALAR]
    current_type = weighting_record.dtype
    reduced_weight = reduced_weight.astype(current_type)

    weighting_record[previos_idx:current_idx] = reduced_weight
    series_hdf_reduction.flush()

def write_patches_information(particle_species, num_particles, num_particles_offset):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    dset = Dataset(numpy.dtype("uint64"), extent=[len(num_particles)])
    particle_species.particle_patches["numParticles"][SCALAR].reset_dataset(dset)
    particle_species.particle_patches["numParticlesOffset"][SCALAR].reset_dataset(dset)

    for i in range(0, len(num_particles)):
        particle_species.particle_patches["numParticles"][SCALAR].store(i, numpy.array([num_particles[i]], dtype=numpy.ulonglong))

    for i in range(0, len(num_particles_offset)):
        particle_species.particle_patches["numParticlesOffset"][SCALAR].store(i, numpy.array([num_particles_offset[i]], dtype=numpy.ulonglong))


def copy_record_values_from_draft(record_draft, record_main, series_hdf_reduction, result_size,
                                  particle_species_reduction):

    for vector in record_draft:
        copy_values = record_draft[vector][0:result_size]
        series_hdf_reduction.flush()
        record_main[vector][0:result_size] = copy_values
        series_hdf_reduction.flush()

    del record_draft


def copy_main_version(series_hdf_reduction, particle_species, particle_species_reduction,result_size):

    for record_name, record in particle_species.items():
        draft_record_name = record_name + "_copy"
        record_main = particle_species_reduction[record_name]
        record_draft = particle_species_reduction[draft_record_name]

        copy_record_values_from_draft(record_draft, record_main, series_hdf_reduction, result_size,
                                      particle_species_reduction)

    for record_name, record in particle_species.items():
        draft_record_name = record_name + "_copy"
        del particle_species_reduction[draft_record_name]


def get_mass(particle_species):
    assert(is_vector_exist("mass", particle_species))
    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    mass_value = particle_species["mass"][SCALAR][0]
    return mass_value


def get_chunk_sizes(particle_species, series_hdf):

    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR

    dataset_size_max = particle_species["position"]["x"].shape[0]

    max_chunk_size = 1e6

    ranges_patches = []

    if particle_species.particle_patches.num_patches > 0:
        ranges_patches = particle_species.particle_patches["numParticlesOffset"][SCALAR].load()
        series_hdf.flush()
        ranges_patches = numpy.append(ranges_patches, dataset_size_max)
    else:
        num_full_chunks = dataset_size_max/max_chunk_size
        if num_full_chunks > 1:
            ranges_patches.append(0)
            for i in range(0, math.floor(num_full_chunks)):
                ranges_patches.append((i + 1) * max_chunk_size)
            ranges_patches = numpy.append(ranges_patches, dataset_size_max)
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


def get_unit_SI(record):

    unit_SI = []

    for idx, value in record.items():
        unit_SI.append(value.unit_SI)

    return unit_SI

def get_coordinates(series_hdf, particle_species, idx_start, idx_end):
    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    position_offset = particle_species["positionOffset"]
    position = particle_species["position"]


    weights = particle_species["weighting"][SCALAR][idx_start: idx_end]
    series_hdf.flush()

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
    SCALAR = openpmd_api.Mesh_Record_Component.SCALAR
    weights = particle_species["weighting"][SCALAR][idx_start: idx_end]
    series_hdf.flush()

    for record_component_name, record_component in particle_species.items():

        if record_component_name == "position" or record_component_name == "weighting" \
                or record_component_name == "positionOffset":
            continue

        if 'macroWeighted' in record_component.attributes:
            weighted_values = get_macroweighted(series_hdf, record_component, weights, idx_start, idx_end)
        else:
            weighted_values = get_non_transformed_values(series_hdf, record_component, idx_start, idx_end)

        if 'unitSI' in record_component.attributes:
            unit_SI = get_unit_SI(record_component)
            absolute_values = get_absolute_values(weighted_values, unit_SI)
        else:
            absolute_values = weighted_values

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

    if args.algorithm == 'vranic':
        mass_value = get_mass(particle_species)
        algorithm.parameters.mass = mass_value

    for i in range(0, len(ranges_patches) - 1):
        print("I == "+ str(i))
        idx_start = int(ranges_patches[i])
        idx_end = int(ranges_patches[i + 1])
        if idx_start == idx_end:
            continue

        absolute_coordinates, unit_si_position, unit_si_offset = get_coordinates(series_hdf, particle_species, idx_start, idx_end)
        data, dict_data_indexes, last_idx = get_data(series_hdf, particle_species, weights, idx_start, idx_end)

        dict_data_indexes["position"] = [last_idx, last_idx + len(absolute_coordinates[0])]

        for i in range(0, len(absolute_coordinates[0])):
            data.append(absolute_coordinates[:, i])

        data = numpy.transpose(data)

        weights_curent = weights[SCALAR][idx_start: idx_end]
        series_hdf.flush()
        reduced_data, reduced_weight = algorithm._run(data, weights_curent, dict_data_indexes)
        relative_coordinates, offset = get_relative_coordinates(reduced_data, dict_data_indexes, unit_si_offset,
                                 unit_si_position)
        relative_data = get_relative_data(reduced_data, particle_species, dict_data_indexes, reduced_weight)

        new_num_particles.append(len(reduced_weight))
        current_idx += len(reduced_weight)


        write_draft_position(particle_species_reduction, series_hdf_reduction, relative_coordinates,
                             offset, previos_idx, current_idx)

        write_draft_copy(relative_data, reduced_weight, dict_data_indexes, series_hdf_reduction, particle_species_reduction,
                         previos_idx, current_idx)
        previos_idx += len(reduced_weight)
        result_size = previos_idx


    new_num_particles_offset = numpy.cumsum(new_num_particles[0:len(new_num_particles) - 1], dtype=int)
    new_num_particles_offset = numpy.insert(new_num_particles_offset, 0, 0)

    create_dataset_structures(particle_species, particle_species_reduction, result_size)
    make_particle_patches_structure(particle_species, particle_species_reduction)
    write_patches_information(particle_species_reduction, new_num_particles, new_num_particles_offset)

    copy_main_version(series_hdf_reduction, particle_species, particle_species_reduction, result_size)

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

if __name__ == "__main__":
    """ Parse arguments from command line     
    
    """

    parser = argparse.ArgumentParser(description="main reduction", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-algorithm", metavar='algorithm', type=str,
                        help= textwrap.dedent('''\
                                      random : random thinning 
                                      number_conservative : number conservative thinning
                                      number_conservative : number conservative thinning
                                      kmeans : kmeans merging
                                      voronoi : voronoi merging
                                      voronoi_prob : voronoi algorithm, depending on ratio of deleted particles
                                      vranic : vranic-based merging
                                       '''))

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file to be reducted")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="result reduction hdf file")

    parser.add_argument("-ratio_deleted_particles", metavar='ratio_deleted_particles', type=float,
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
        divide_particles = 20
        parameters = Voronoi_probabilistic_algorithm.Voronoi_probabilistic_algorithm_parameters(1. - args.ratio_deleted_particles, divide_particles)
        base_reduction_function(args.hdf, args.hdf_re, "voronoi_prob", parameters)

    elif args.algorithm == 'random':
        parameters = Random_thinning_algorithm.Random_thinning_algorithm_parameters(args.ratio_deleted_particles)
        base_reduction_function(args.hdf, args.hdf_re, "random", parameters)

    elif args.algorithm == 'number_conservative':
        parameters = Number_conservative_thinning_algorithm.Number_conservative_thinning_algorithm_parameters(args.ratio_deleted_particles)
        base_reduction_function(args.hdf, args.hdf_re, "number_conservative", parameters)

    elif args.algorithm == 'energy_conservative':
        parameters = Energy_conservative_thinning_algorithm.Energy_conservative_thinning_algorithm_parameters(args.ratio_deleted_particles)
        base_reduction_function(args.hdf, args.hdf_re, "energy_conservative", parameters)

    elif args.algorithm == 'kmeans':
        parameters = k_means_clustering_algorithm.K_means_clustering_algorithm_parameters(args.ratio_deleted_particles)
        base_reduction_function(args.hdf, args.hdf_re, "kmeans", parameters)

    elif args.algorithm == 'kmeans_avg':
        parameters = k_means_merge_average_algorithm.K_means_merge_average_algorithm_parameters(args.reduction_percent)
        base_reduction_function(args.hdf, args.hdf_re, "kmeans_avg", parameters)

    elif args.algorithm == 'leveling':
        parameters = Leveling_thinning_algorithm.Leveling_thinning_algorithm_parameters(args.leveling_coefficient)
        base_reduction_function(args.hdf, args.hdf_re, "leveling", parameters)

    elif args.algorithm == 'vranic':
        parameters = Vranic_algorithm.Vranic_merging_algorithm_parameters(args.momentum_tol, args.position_lol)
        base_reduction_function(args.hdf, args.hdf_re, "vranic", parameters)
