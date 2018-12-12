
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


class _Voronoi_cell:
    """ Используется для хранения points in a cell"""

    def __init__(self, points):
        self.points = points

    def get_coeff_var(self):
        first_key = list(self.points.keys())[0]
        print('first == ' + str(first_key))
        demention = len(self.points[first_key][0.coords[0])
        print('demention == ' + str(demention))
    def divide(self, component):
        array = []
        return array
        """component - index of coordinate to use for subdivision, this function returns two Voronoi Cells"""

    def merge(self):
        pass
        """Get one point out of all points in the cell, return Point"""

def merge_coordinates(demention, merged_points):

    size = len(merged_points['position'])

    weights = 0.
    full_values = []
    for i in range(0, demention):

        value = 0
        for point in merged_points['position']:
            value += point.coords[i]
        full_values.append(value/size)

    for point in merged_points['momentum']:
        weights += point.weight

    result = []
    result.append(Point(full_values, weights))
    return result

