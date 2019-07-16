

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

