import numpy as np
import h5py
from shutil import copyfile
import read_hdf_file
import random
import math

class Segment:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def is_in_segment(self, value):
        if value >=self.start and value <= self.end:
            return True
        else:
            return False


class Momentum_cell:

    def __init__(self):
        self.segment_x = 0.
        self.segment_y = 0.
        self.segment_z = 0.

        self.momentums = []

    def get_idixes(self):
        self.momentums = np.array(self.momentums)
        indexes = np.array(self.momentums[:,3], dtype=('int'))
        return indexes


def create_sorted_momentum_list(data, dimension):
    
    size = len(data)
    indexes = [i for i in range(len(data))]
    result = []
    for i in range(0, size):
        result.append([data[i][dimension.dimension_position], data[i][dimension.dimension_position+1],
                       data[i][dimension.dimension_position+2], indexes[i]])

    sort_list = sorted(result, key=lambda e: (e[0], e[1], e[2]))

    return sort_list


def get_segment(momentum_min, idx, step):
    start = momentum_min + (idx - 1) * step
    end = momentum_min + (idx) * step
    segment = Segment(start, end)
    return segment


def get_position_idx3d(x_patch, y_patch, z_patch, size_y, size_z):
    return (x_patch * size_y + y_patch) * size_z + z_patch


def get_position_idx2d(x_patch, y_patch, size_y):
    return x_patch * size_y + y_patch


def get_cell_idx(max_coord, min_coord, separator, x_current):
    """ Get name of particles group """
    lenght = max_coord - min_coord
    if lenght == 0:
        return 0
    else:
        return max(0, min(int((x_current - min_coord) * separator / lenght), separator - 1))


def create_momentum_cells(sorted_momentum_list, tolerance_momentum, segment_x, segment_y, segment_z):

    x_momentum_min = segment_x.start
    x_momentum_max = segment_x.end

    y_momentum_min = segment_y.start
    y_momentum_max = segment_y.end

    z_momentum_min = segment_z.start
    z_momentum_max = segment_z.end

    num_x_steps = int((x_momentum_max - x_momentum_min)/ tolerance_momentum) + 1
    num_y_steps = int((y_momentum_max - y_momentum_min) / tolerance_momentum) + 1
    num_z_steps = int((z_momentum_max - z_momentum_min) / tolerance_momentum) + 1

    size_cells = num_x_steps * num_y_steps * num_z_steps
    cells = [Momentum_cell() for i in range(size_cells)]

    for i in range(0, len(sorted_momentum_list)):
        x_idx = get_cell_idx(x_momentum_max, x_momentum_min, num_x_steps, sorted_momentum_list[i][0])

        y_idx = get_cell_idx(y_momentum_max, y_momentum_min, num_y_steps, sorted_momentum_list[i][1])
        z_idx = get_cell_idx(z_momentum_max, z_momentum_min, num_z_steps, sorted_momentum_list[i][2])

        current_idx = get_position_idx3d(x_idx, y_idx, z_idx, num_y_steps, num_z_steps)
        cells[current_idx].momentums.append(sorted_momentum_list[i])

    return cells


def get_start_ranges(data,dimensions):

    x_min = min(data[:,dimensions.dimension_position])
    x_max = max(data[:,dimensions.dimension_position])

    x_segment = Segment(x_min, x_max)

    y_min = min(data[:,dimensions.dimension_position+1])
    y_max = max(data[:,dimensions.dimension_position+1])

    y_segment = Segment(y_min, y_max)

    z_min = min(data[:,dimensions.dimension_position+2])
    z_max = max(data[:,dimensions.dimension_position+2])

    z_segment = Segment(z_min, z_max)

    return x_segment, y_segment, z_segment


def get_energy_value(momentum, mass):

    c = 299792458.
    if mass == 0:
        return c * np.sqrt(momentum[0]**2 + momentum[1]**2 + momentum[2]**2)
    momentum_norm = [momentum[0]/(mass * c), momentum[1]/(mass * c), momentum[2]/(mass * c) ]
    norm_vector = momentum_norm[0] * momentum_norm[0] + momentum_norm[1] * momentum_norm[1] + momentum_norm[2] * momentum_norm[2]
    energy = np.sqrt(1. + norm_vector) * (mass * c * c)
    return energy

def get_longitudinal_momentum(momentums, weights):

    values_x = 0.
    values_y = 0.
    values_z = 0.
    sum_weights = 0.

    for i in range(0, len(momentums)):
        values_x += momentums[i][0] * weights[i]
        values_y += momentums[i][1] * weights[i]
        values_z += momentums[i][2] * weights[i]
        sum_weights += weights[i]

    result = np.array([values_x, values_y, values_z])
    return result / sum_weights


def get_weighted_energy(momentums, weights, mass):

    sum_energy = 0.

    for i in range(0, len(momentums)):
        energy = get_energy_value(momentums[i], mass)
        sum_energy = sum_energy + energy * weights[i]

    return sum_energy


def get_momentum_vector_lenght(energy_t, mass):
    c = 299792458.
    if mass == 0:
        return energy_t / c
    return mass * c * math.sqrt(energy_t * energy_t / mass**2 / c**4  - 1.)


def get_transverse_momentum_length(momentums, mass, weights, longitudinal_momentum):
    sum_weights = sum(weights)

    energy_t = get_weighted_energy(momentums, weights, mass)
    energy_t = energy_t / sum_weights
    norm_vector = get_momentum_vector_lenght(energy_t, mass)
    long_norm_squared = longitudinal_momentum[0]**2 + longitudinal_momentum[1]**2 + longitudinal_momentum[2]**2
    return np.sqrt(norm_vector**2 - long_norm_squared)


def get_perp_vec(longitudinal_momentum, x_segment, y_segment, z_segment):
    # Get vector perpendicular to longitudinal momentum and in the plane defined by
    # longitudinal momentum and one of the cell diagonal

    long_momentum_length = np.sqrt(longitudinal_momentum[0]**2+longitudinal_momentum[1]**2+longitudinal_momentum[2]**2)
    long_momentum_norm = longitudinal_momentum / long_momentum_length

    cell_sizex = x_segment.end - x_segment.start
    cell_sizey = y_segment.end - y_segment.start
    cell_sizez = z_segment.end - z_segment.start
    diagonal_vec = np.array([cell_sizex,cell_sizey,cell_sizez])

    ## In principle we should cancel the merging when the cell has a zero size in two of the dimensions
    ## because in this case we cannot conserve both momentum and energy for a massive particle
    ## We'll ignore this rare case for now for simplicity

    sorted_directions = np.argsort(np.abs(long_momentum_norm))

    # This should ensure that the diagonal is not colinear with total_momentum
    diagonal_vec[sorted_directions[2]] = -diagonal_vec[sorted_directions[2]] 

    dot_product = long_momentum_norm[0]*diagonal_vec[0] + long_momentum_norm[1]*diagonal_vec[1] + long_momentum_norm[2]*diagonal_vec[2]

    # This calculates the double vector product total_momentum_norm x (total_momentum_norm x diagonal_vec)
    perp_vec = np.array([dot_product*long_momentum_norm[0]-diagonal_vec[0],
                         dot_product*long_momentum_norm[1]-diagonal_vec[1],
                         dot_product*long_momentum_norm[2]-diagonal_vec[2]])

    perp_vec_length = np.sqrt(perp_vec[0]**2+perp_vec[1]**2+perp_vec[2]**2)

    return perp_vec/perp_vec_length


def recalculate_momentum(momentums, weights, mass, x_segment, y_segment, z_segment):

    longitudinal_momentum = get_longitudinal_momentum(momentums, weights)
    transverse_momentum_length = get_transverse_momentum_length(momentums, mass, weights, longitudinal_momentum)
    perp_vec = get_perp_vec(longitudinal_momentum, x_segment, y_segment, z_segment)
    
    return longitudinal_momentum + transverse_momentum_length*perp_vec, longitudinal_momentum - transverse_momentum_length*perp_vec


def calculate_result_points(data, weights, idxes_array):
    first_point = []
    second_point = []
    sum_weights = sum(weights)
    result_weight = [sum_weights / 2., sum_weights / 2.]

    idxes_array.sort()

    idexes_coordinates = random.sample(idxes_array.data, 2)
    first_coordinates = data[idexes_coordinates[0]]
    second_coordinates = data[idexes_coordinates[1]]

    return first_coordinates, second_coordinates, result_weight


def merge_into_points(first_coordinates, second_coordinates, first_momentum, second_momentum, dimension):
    if dimension.dimension_position == 3:
        first_point = [first_coordinates[0], first_coordinates[1], first_coordinates[2],  first_momentum[0], first_momentum[1], first_momentum[2]]
        second_point = [second_coordinates[0], second_coordinates[1], second_coordinates[2], second_momentum[0],
                        second_momentum[1], second_momentum[2]]
    elif dimension.dimension_position == 2:
        merged_points = np.array([[first_coordinates[0], first_coordinates[1],  first_momentum[0], first_momentum[1], first_momentum[2]],
                      [second_coordinates[0], second_coordinates[1], second_momentum[0], second_momentum[1], second_momentum[2]]])
    else:
        assert(0)

    return merged_points


def recount_cells(data, weights, momentum_cells, mass, dimension, x_segment, y_segment, z_segment):
    result = np.zeros((0,data.shape[1]))
    weights_result = np.zeros((0))

    for i in range(0, len(momentum_cells)):
        if len(momentum_cells[i].momentums) != 0:
            idxes_array = momentum_cells[i].get_idixes()
            if len(idxes_array) > 2:
                first_coordinates, second_coordinates, result_weight = calculate_result_points(data, weights[idxes_array], idxes_array)
                print("momentum_cells[i].momentums",momentum_cells[i].momentums)
                first_momentum, second_momentum = recalculate_momentum(momentum_cells[i].momentums, weights[idxes_array], mass, x_segment, y_segment, z_segment)
                print("first_momentum",first_momentum)
                merged_points = merge_into_points(first_coordinates, second_coordinates, first_momentum, second_momentum, dimension)
                result = np.append(result,merged_points,axis=0)
                weights_result = np.append(weights_result,result_weight[0])
                weights_result = np.append(weights_result,result_weight[1])

            else:
                result = np.append(result,data[idxes_array],axis=0)
                weights_result = np.append(weights_result,weights[idxes_array])

    return result, weights_result

def sortParticlesByBins(data, weights, tolerance_position, dimension):
    """Returns a list of data and weight arrays sorted by bins of size tolerance_position in each direction"""

    if dimension.dimension_position == 2:
        return sortParticlesByBins2D(data, weights, tolerance_position)
    elif dimenion.dimension_position == 3:
        return sortParticlesByBins3D(data, weights, tolerance_position)
    else:
        assert(0)

def sortParticlesByBins2D(data, weights, tolerance_position):
    """Version of sortParticlesByBins for 2D simulations"""

    x_min = min(data[:,0])
    x_max = max(data[:,0])
    y_min = min(data[:,1])
    y_max = max(data[:,1])

    num_x_steps = int((x_max - x_min)/ tolerance_position) + 1
    num_y_steps = int((y_max - y_min)/ tolerance_position) + 1

    numbins = num_x_steps * num_y_steps

    result_data = [np.zeros((0,data.shape[1])) for i in range(numbins)]
    weight_data = [np.zeros((0)) for i in range(numbins)]

    for i in range(0, data.shape[0]):
        x_idx = get_cell_idx(x_max, x_min, num_x_steps, data[i][0])
        y_idx = get_cell_idx(y_max, y_min, num_y_steps, data[i][1])

        current_idx = get_position_idx2d(x_idx, y_idx, num_y_steps)
        result_data[current_idx] = np.append(result_data[current_idx],data[i,:][None,:],axis=0)
        weight_data[current_idx] = np.append(weight_data[current_idx],weights[i])

    return result_data, weight_data

def sortParticlesByBins3D(data, weights, tolerance_position):
    """Version of sortParticlesByBins for 3D simulations"""

    x_min = min(data[:,0])
    x_max = max(data[:,0])
    y_min = min(data[:,1])
    y_max = max(data[:,1])
    z_min = min(data[:,2])
    z_max = max(data[:,2])

    num_x_steps = int((x_max - x_min)/ tolerance_position) + 1
    num_y_steps = int((y_max - y_min)/ tolerance_position) + 1
    num_z_steps = int((z_max - z_min)/ tolerance_position) + 1

    numbins = num_x_steps * num_y_steps * num_z_steps

    result_data = [np.zeros((0,data.shape[1])) for i in range(numbins)]
    weight_data = [np.zeros((0)) for i in range(numbins)]

    for i in range(0, data.shape[0]):
        x_idx = get_cell_idx(x_max, x_min, num_x_steps, data[i][0])
        y_idx = get_cell_idx(y_max, y_min, num_y_steps, data[i][1])
        z_idx = get_cell_idx(z_max, z_min, num_z_steps, data[i][2])

        current_idx = get_position_idx3d(x_idx, y_idx, z_idx, num_y_steps, num_z_steps)
        result_data[current_idx] = np.append(result_data[current_idx],data[i,:][None,:],axis=0)
        weight_data[current_idx] = np.append(weight_data[current_idx],weights[i])

    return result_data, weight_data


class Vranic_merging_algorithm_parameters:
    """Tolerance is array-like, first -- coordinate tolerance, second -- momentum tolerance"""

    def __init__(self, tolerance_momentum, tolerance_position):
        self.tolerance_momentum = tolerance_momentum
        self.tolerance_position = tolerance_position
        self.dimension = None
        self.mass = None


class Vranic_merging_algorithm:
    """Main algorithm. Parameters is VoronoiMergingAlgorithmParameters """

    def __init__(self, parameters):
        self.parameters = parameters

    def _run(self, data, weigths, dimension):
        self.parameters.dimension = dimension
        mass = self.parameters.mass
        data = np.array(data)
        weigths = np.array(weigths)

        data_sorted_spatially, weights_sorted_spatially = sortParticlesByBins(data, 
                                                 weigths, self.parameters.tolerance_position, dimension)

        result_data = np.zeros((0,data.shape[1]))
        result_weight = np.zeros((0))

        for i in range(len(data_sorted_spatially)):
            if data_sorted_spatially[i].shape[0] == 0:
                continue
            x_segment, y_segment, z_segment = get_start_ranges(data_sorted_spatially[i],dimension)

            sorted_momentum = create_sorted_momentum_list(data_sorted_spatially[i], dimension)
            cells = create_momentum_cells(sorted_momentum, self.parameters.tolerance_momentum, x_segment, y_segment, z_segment)
            data_bin, weights_bin = recount_cells(data_sorted_spatially[i], weights_sorted_spatially[i], cells, mass, dimension, x_segment, y_segment, z_segment)
            print("result_data.shape",result_data.shape)
            result_data = np.append(result_data,data_bin,axis=0)
            result_weight = np.append(result_weight,weights_bin)
        return result_data, result_weight


