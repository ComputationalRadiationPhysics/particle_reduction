import numpy
import math
import statistics
import random


class Voronoi_probabilistic_algorithm_parameters:
    """"""

    def __init__(self, reduction_percent, ratio_left_particles):
        self.reduction_percent = reduction_percent
        self.ratio_left_particles = ratio_left_particles


class Voronoi_probabilistic_algorithm:
    """Main algorithm. Parameters is Voronoi_probabilistic_algorithm_parameters """

    def __init__(self, parameters):
        self.parameters = parameters
        self.dimensions = None

    def _run(self, data, weigths, dict_data_indexes):

        """Points is a collection of Point"""

        print("data len : "+ str(len(data)))
        return _merge(data, weigths, self.parameters, dict_data_indexes)


class _Voronoi_cell:

    """Used to store points in Voronoi cell"""

    def __init__(self, data, weigths, expected_number_of_particles, size_of_divide_particles):


        self.vector = numpy.array(data)
        self.weights = numpy.array(weigths)
        self.size_of_divide_particles = size_of_divide_particles
        self.expected_number_of_particles = expected_number_of_particles



    def get_coeff_var(self, dict_data_indexes):

        """Get max variance coefficient of Voronoi cell"""

        avg_values = []

        range_momentum = range(dict_data_indexes["momentum"][0], dict_data_indexes["momentum"][1])
        range_position = range(dict_data_indexes["position"][0], dict_data_indexes["position"][1])


        for i in range_momentum:
            values_dimension = self.vector[ : , i]

            std = numpy.sqrt(weighted_variance(values_dimension, self.weights))
            avg_values.append(std)


        for i in range_position:
            values_dimension = self.vector[ : , i]

            std = numpy.sqrt(weighted_variance(values_dimension, self.weights))
            avg_values.append(std)


        max_idx, max_avg = get_max_coef(avg_values, range_momentum, range_position)


        return max_idx, max_avg


    def divide(self, devide_hyperlane, reduction_percent):

        """

        Devide Voronoi cell into two Voronoi cells
        max_idx --
        max_key --

        """

        median = statistics.median(self.vector[:, devide_hyperlane])

        first = []
        weights_first = []
        second = []
        weights_second = []

        for i in range(0, len(self.vector)):
            if self.vector[i][devide_hyperlane] > median:
                first.append(self.vector[i])
                weights_first.append(self.weights[i])
            else:
                second.append(self.vector[i])
                weights_second.append(self.weights[i])
        first_number_expected_particles = 0
        second_number_expected_particles = 0
        if len(self.vector) > self.size_of_divide_particles:

            first_number_expected_particles = reduction_percent * len(first)
            second_number_expected_particles = reduction_percent * len(second)

        else:
            b = (self.expected_number_of_particles + len(self.vector)) / 2.0
            first_number_expected_particles = len(first) * b / len(self.vector)
            second_number_expected_particles = len(second) * b / len(self.vector)

        if len(self.vector) < first_number_expected_particles:
            print("len(self.vector)" + str(len(self.vector)))
            print("first_number_expected_particles" + str(first_number_expected_particles))
            assert(False)

        if len(self.vector) < second_number_expected_particles:
            print("len(self.vector)" + str(len(self.vector)))
            print("second_number_expected_particles" + str(second_number_expected_particles))
            assert(False)


        first_cell = _Voronoi_cell(first, weights_first, first_number_expected_particles,
                                  self.size_of_divide_particles)
        second_cell = _Voronoi_cell(second, weights_second, second_number_expected_particles,
                                   self.size_of_divide_particles)

        return first_cell, second_cell

        #"""component - index of coordinate to use for subdivision, this function returns two Voronoi Cells"""

    def merge(self, dict_data_indexes):
        """ Merge Voronoi cell into one point """


        self.vector, self.weights = merge_points(dict_data_indexes, self.vector, self.weights)


def merge_points(dict_data_indexes, vector, weights_vector):
    """ Merge coordinates into one point """

    values_vector = []
    dimension = len(vector[0])
    sum_weights = numpy.sum(weights_vector, dtype=float)
    vector = numpy.array(vector)


    range_momentum = range(dict_data_indexes["momentum"][0], dict_data_indexes["momentum"][1])
    range_position = range(dict_data_indexes["position"][0], dict_data_indexes["position"][1])

    for i in range(0, dimension):
        if i in range_momentum:
            values_vector.append(numpy.average(vector[:, i], weights=weights_vector))
        elif i in range_position:
            values_vector.append(numpy.average(vector[:, i], weights=weights_vector))
        else:
            values_vector.append(vector[0][i])

    return values_vector, sum_weights

def weighted_variance(values, weights):
    """
    Return the weighted variance.

    values, weights -- Numpy ndarrays with the same shape.

    """

    if len(values) == 1:
        return 0

    weighted_average = 0

    for i in range(0, len(values)):
        weighted_average = weighted_average + values[i] * weights[i]

    sum_weights = sum(weights)

    weighted_sq_average = 0
    for i in range(0, len(values)):
        weighted_sq_average = weighted_sq_average + values[i] * values[i] * weights[i]

    weighted_average = weighted_average / sum_weights
    weighted_sq_average = weighted_sq_average / sum_weights
    # Fast and numerically precise:
    variance = weighted_sq_average - weighted_average * weighted_average
    return abs(variance)


def get_max_coef(avg_values, range_momentum, range_position):
    """ Find max coefficient of variance """

    max_value = float("-inf")
    max_idx = -1

    for i in range(0, len(avg_values)):
        if avg_values[i] > max_value:
            max_value = avg_values[i]
            max_idx = i

    if max_idx < len(range_momentum):
        max_idx = range_momentum[0] + max_idx

    else:
        max_idx = max_idx - len(range_momentum) + range_position[0]

    return max_idx, max_value

def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.

    """
    weighted_average = numpy.average(values, weights=weights)
    # Fast and numerically precise:
    variance = numpy.average((values-weighted_average)**2, weights=weights)
    return math.sqrt(variance)


def _merge(data, weights, parameters, dict_data_indexes):
    """
    Merging algorithm:
    points -- original points
    parametes -- input parameters for Voronoi algorithm(tolerances)

    """
    initial_cell = _Voronoi_cell(data, weights, parameters.reduction_percent * len(data),
                                parameters.ratio_left_particles)

    print(parameters.reduction_percent)

    result_vector = []
    result_weights = []
    cells = [initial_cell]
    while len(cells) > 0:
        cell = cells[0]
        max_idx, max_avg = cell.get_coeff_var(dict_data_indexes)
        if len(cell.vector) < 2:
            for i in range(0, len(cell.vector)):
                result_vector.append(cell.vector[i])

            for i in range(0, len(cell.weights)):
                result_weights.append(cell.weights[i])
            cells.remove(cells[0])
            continue

        needs_subdivision = check_statistical_subdivision(cell, parameters.ratio_left_particles)
       # needs_subdivision = check_needs_subdivision(parameters, max_avg, max_idx, dimensions)

        if needs_subdivision:
            first_part_cell, secound_part_cell = cell.divide(max_idx, parameters.reduction_percent)
            if len(first_part_cell.vector) == 0:
                if len(secound_part_cell.vector) != 0:
                    secound_part_cell.merge(dict_data_indexes)
                    result_vector.append(secound_part_cell.vector)
                    result_weights.append(secound_part_cell.weights)

            if len(secound_part_cell.vector) == 0:
                if len(first_part_cell.vector) != 0:
                    first_part_cell.merge()
                    result_vector.append(first_part_cell.vector)
                    result_weights.append(first_part_cell.weights)

            if len(secound_part_cell.vector) != 0 and len(first_part_cell.vector) != 0:
                cells.append(secound_part_cell)
                cells.append(first_part_cell)
        else:

            cell.merge(dict_data_indexes)
            result_vector.append(cell.vector)
            result_weights.append(cell.weights)

        cells.remove(cells[0])

    result_vector = numpy.asarray(result_vector)
    result_weights = numpy.asarray(result_weights)

    return result_vector, result_weights



def check_statistical_subdivision(cell, size_of_divide_particles):


    if len(cell.vector) > size_of_divide_particles:
        return True

    b_value = (cell.expected_number_of_particles + len(cell.vector))/2.

    p_value = (cell.expected_number_of_particles - 1.) / (b_value - 1.)


    if len(cell.vector) == 3:
        p_value = (cell.expected_number_of_particles - 1.) / 2.0

    if len(cell.vector) == 2:
        p_value = cell.expected_number_of_particles - 1.

    random_value = random.uniform(0, 1)

    return random_value < p_value




