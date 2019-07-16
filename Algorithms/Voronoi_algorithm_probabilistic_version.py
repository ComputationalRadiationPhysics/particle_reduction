

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
