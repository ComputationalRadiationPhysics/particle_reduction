import math
import statistics


class Voronoi_probabilistic_algorithm_parameters:
    """"""

    def __init__(self, tolerance, reduction_percent, ratio_left_particles):
        self.tolerance = tolerance
        self.reduction_percent = reduction_percent
        self.ratio_left_particles = ratio_left_particles


class Voronoi_probabilistic_algorithm:
    """Main algorithm. Parameters is Voronoi_probabilistic_algorithm_parameters """

    def __init__(self, parameters):
        self.parameters = parameters
        self.dimensions = None

    def _run(self, data, weigths):

        return _merge(data, weigths, self.parameters, self.dimensions)


class _Voronoi_cell:

    """Used to store points in Voronoi cell"""

    def __init__(self, data, weigths, expected_number_of_particles, size_of_divide_particles):


        self.vector = numpy.array(data)
        self.weights = numpy.array(weigths)
        self.size_of_divide_particles = size_of_divide_particles
        self.expected_number_of_particles = expected_number_of_particles

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

    def divide(self, devide_hyperlane, ratio_left_partilces):

        """

        Devide Voronoi cell into two Voronoi cells
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

        if len(self.vector) > self.size_of_divide_particles:
            first_number_expected_particles = ratio_left_partilces * len(first)
            second_number_expected_particles = ratio_left_partilces * len(second)

        else:
            b = (self.expected_number_of_particles + len(self.vector))/ 2.0
            first_number_expected_particles = len(first)/len(self.vector) * b
            second_number_expected_particles = len(second) / len(self.vector) * b


        first_cell = _Voronoi_cell(first, weights_first, first_number_expected_particles,
                                  self.size_of_divide_particles)
        second_cell = _Voronoi_cell(second, weights_second, second_number_expected_particles,
                                   self.size_of_divide_particles)

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

