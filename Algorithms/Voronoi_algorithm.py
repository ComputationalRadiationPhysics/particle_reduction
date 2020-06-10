import numpy
import math
import statistics


class VoronoiMergingAlgorithmParameters:
    """Tolerance is array-like, first -- coordinate tolerance, second -- momentum tolerance"""

    def __init__(self, tolerance):
        self.tolerance = tolerance


class VoronoiMergingAlgorithm:
    """Main algorithm. Parameters is VoronoiMergingAlgorithmParameters """

    def __init__(self, parameters):
        self.parameters = parameters
        self.dimensions = None

    def _run(self, data, weigths, dimensions):

        """Points is a collection of Point"""
        self.dimensions = dimensions
        return _merge(data, weigths, self.parameters, self.dimensions)


class _VoronoiCell:

    """Used to store points in Voronoi cell"""

    def __init__(self, data, weigths):

        """Points is array of Points"""

        self.vector = numpy.array(data)
        self.weights = numpy.array(weigths)

    def get_coeff_var(self):

        """Get max variance coefficient of Voronoi cell"""

        dimension = len(self.vector[0])

        avg_values = []

        for i in range(0, dimension):
            values_dimension = []
            weights_dimension = []
            for j in range(0, len(self.vector)):
                values_dimension.append(self.vector[j][i])
                weights_dimension.append(self.weights[j])

            std = weighted_std(values_dimension, weights_dimension)
            avg_values.append(std)

        max_idx, max_avg = get_max_coef(avg_values)

        return max_idx, max_avg

    def divide(self, devide_hyperlane):

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

        first_cell = _VoronoiCell(first, weights_first)
        second_cell = _VoronoiCell(second, weights_second)

        return first_cell, second_cell

        #"""component - index of coordinate to use for subdivision, this function returns two Voronoi Cells"""

    def merge(self):
        """ Merge Voronoi cell into one point """

        dimension = len(self.vector[0])

        self.vector, self.weights = merge_points(dimension, self.vector, self.weights)


def merge_points(dimension, vector, weights_vector):
    """ Merge coordinates into one point """

    values_vector = []
    sum_weights = numpy.sum(weights_vector, dtype=float)
    vector = numpy.array(vector)

    for i in range(0, dimension):
        values_vector.append(numpy.average(vector[:, i], weights=weights_vector))

    return values_vector, sum_weights


def get_max_coef(avg_values):
    """ Find max coefficient of variance """

    max_value = float("-inf")
    max_idx = -1

    for i in range(0, len(avg_values)):
        if avg_values[i] > max_value:
            max_value = avg_values[i]
            max_idx = i

    return max_idx, max_value


def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.

    """

    weighted_average = 0

    for i in range(0, len(values)):
        weighted_average = weighted_average + values[i] * weights[i]

    sum_weights = sum(weights)

    weighted_sq_average = 0
    for i in range(0, len(values)):
        weighted_sq_average = weighted_sq_average + values[i] * values[i] *weights[i]

    weighted_average = weighted_average / sum_weights
    weighted_sq_average = weighted_sq_average / sum_weights
    # Fast and numerically precise:
    variance = weighted_sq_average - weighted_average * weighted_average

    return variance


def _merge(data, weights, parameters, dimensions):
    """
    Merging algorithm:
    points -- original points
    parametes -- input parameters for Voronoi algorithm(tolerances)

    """
    initial_cell = _VoronoiCell(data, weights)

    result_vector = []
    result_weights = []
    cells = [initial_cell]
    while len(cells) > 0:
        cell = cells[0]

        max_idx, max_avg = cell.get_coeff_var()
        needs_subdivision = check_needs_subdivision(parameters, max_avg, max_idx, dimensions)

        if needs_subdivision:
            first_part_cell, secound_part_cell = cell.divide(max_idx)
            if len(first_part_cell.vector) == 0:
                if len(secound_part_cell.vector) != 0:
                    secound_part_cell.merge()
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

            cell.merge()
            result_vector.append(cell.vector)
            result_weights.append(cell.weights)

        cells.remove(cells[0])

    result_vector = numpy.asarray(result_vector)
    result_weights = numpy.asarray(result_weights)

    return result_vector, result_weights


def check_needs_subdivision(parameters, max_avg, max_idx, dimensions):
    """
    Check that Voronoi cell need to devide
    parametrs -- subdivision tolerances
    max_avg -- value of max variance
    max_key -- parameter with max variance

    """

    position_tolerance = parameters[0]
    momentum_tolerance = parameters[1]
    position_vector_idx = dimensions.dimension_position


    if max_idx <= position_vector_idx:
        return max_avg > position_tolerance

    if max_idx > position_vector_idx:
        return max_avg > momentum_tolerance
