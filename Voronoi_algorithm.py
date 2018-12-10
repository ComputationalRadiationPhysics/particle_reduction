
class Point:
    def __init__(self, coords, weight):
        """A weighted point, coords is array-like, weight is float"""
        self.coords = coords
        self.weight = weight


class Voronoi_merging_algorithm_parametrs:
    def __init__(self, tolerance):
        """..."""
        self.tolerance = tolerance


class Voronoi_merging_algorithm:
    def __init__(self, parameters):
        self.parameters = parameters

    def run(self, points):
        """Points is a collection of Point"""
        return _merge(points, self.parameters)

