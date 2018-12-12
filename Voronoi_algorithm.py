
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
        demention = len(self.points[first_key][0].coords)

        avg_keys = {}

        for keys in self.points.keys():
            avg_values = []
            for i in range(0, demention):
                values = []
                weights = []

                for points in self.points[keys]:
                    values.append(points.coords[i])
                    weights.append(points.weight)

                avg = weighted_avg(values, weights)
                avg_values.append(avg)
            avg_keys[keys] = avg_values

        max_idx, max_key, max_avg = get_max_coef(avg_keys)

        return max_idx, max_key, max_avg

    def merge(self):

        first_key = list(self.points.keys())[0]
        demention = len(self.points[first_key][0].coords)

        for key in self.points.keys():

            if key == 'momentum':
                self.points['momentum'] = merge_momentum(demention, self.points)
            if key == 'position':
                self.points['position'] = merge_coordinates(demention, self.points)
 #       """Get one point out of all points in the cell, return Point"""


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


def merge_momentum(demention, merged_points):

    weights = 0.
    full_values = []
    for i in range(0, demention):
        value = 0
        for point in merged_points['momentum']:
            value += point.coords[i]
        full_values.append(value)

    for point in merged_points['momentum']:
        weights += point.weight

    result = []
    result.append(Point(full_values, weights))
    return result


def get_max_coef(avg_keys):

    max_value = float("-inf")
    max_idx = -1
    max_key = ''
    for key in avg_keys:
        for i in range(0, len(avg_keys[key])):
            if avg_keys[key][i] > max_value:
                max_value = avg_keys[key][i]
                max_idx = i
                max_key = key

    return max_idx, max_key, max_value


def weighted_avg(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = numpy.average(values, weights=weights)
    # Fast and numerically precise:
    variance = numpy.average((values-average)**2, weights=weights)
    return math.sqrt(variance)


def _merge(points, parameters):
    initial_cell = _Voronoi_cell(points)
    max_idx, max_key, max_avg = initial_cell.get_coeff_var()
    initial_cell.merge()
    result = []
    cells = [initial_cell]

    while len(cells) > 0:
        cell = cells[0]
        max_idx, max_key, max_avg = cell.get_coeff_var()
        needs_subdivision = check_needs_subdivision(parameters, max_avg, max_key)
       #  compare coeff_var with parameters
        if needs_subdivision:
            first_part_cell, secound_part_cell = cell.divide(max_idx, max_key)
            cells.append(first_part_cell)
            cells.append(secound_part_cell)
        else:
            new_particle = cell.merge()
            result.append(new_particle)
        del cells[0]

    return result


def check_needs_subdivision(parameters, max_avg, max_key):
    # first -- momentum
    # secound -- coord

    position_tolerance = parameters.tolerance[0]
    momentum_tolerance = parameters.tolerance[1]

    if max_key == 'position':
        return max_avg > position_tolerance

    if max_key == 'momentum':
        return max_avg > momentum_tolerance

