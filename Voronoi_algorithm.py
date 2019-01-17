import numpy
import math
import copy


class Point:

    def __init__(self, coords, weight):
        """A weighted point, coords is array-like, weight is float"""

        self.coords = coords
        self.weight = weight


class VoronoiMergingAlgorithmParameters:
    """Tolerance is array-like, first -- coordinate tolerance, second -- momentum tolerance"""

    def __init__(self, tolerance):
        self.tolerance = tolerance


class VoronoiMergingAlgorithm:
    """Main algorithm. Parameters is VoronoiMergingAlgorithmParameters """

    def __init__(self, parameters):
        self.parameters = parameters

    def run(self, points):

        """Points is a collection of Point"""
        return _merge(points, self.parameters)


class _VoronoiCell:

    """Used to store points in Voronoi cell"""

    def __init__(self, points):

        """Points is array of Points"""

        self.points = points

    def get_coeff_var(self):

        """Get max variance coefficient of Voronoi cell"""


        first_key = list(self.points.keys())[0]
        demention = len(self.points[first_key][0].coords)

        avg_keys = {}

        for key in self.points.keys():

            avg_values = []

            for i in range(0, demention):
                values = []
                weights = []
                for points in self.points[key]:
                    values.append(points.coords[i])
                    weights.append(points.weight)

                avg = weighted_avg(values, weights)
                avg_values.append(avg)
            avg_keys[key] = avg_values

        max_idx, max_key, max_avg = get_max_coef(avg_keys)
        return max_idx, max_key, max_avg

    def divide(self, max_idx, max_key):

        """

        Devide Voronoi cell into two Voronoi cells
        max_idx --
        max_key --

        """

        max_value = float("-inf")
        min_value = float("inf")

        for point in self.points[max_key]:
            if max_value < point.coords[max_idx]:
                max_value = point.coords[max_idx]

            if min_value > point.coords[max_idx]:
                min_value = point.coords[max_idx]

        middle_value = (max_value + min_value)/2.

        array_first = {}
        array_secound = {}

        array_first['position'] = []
        array_first['momentum'] = []
        array_secound['position'] = []
        array_secound['momentum'] = []

        for i in range(0, len(self.points[max_key])):
            if self.points[max_key][i].coords[max_idx] > middle_value:
                for keys in self.points:
                       array_first[keys].append(self.points[keys][i])
            else:
                for keys in self.points:
                    array_secound[keys].append(self.points[keys][i])

        first_cell = _VoronoiCell(array_first)
        secound_cell = _VoronoiCell(array_secound)
        return first_cell, secound_cell
        #"""component - index of coordinate to use for subdivision, this function returns two Voronoi Cells"""

    def merge(self):
        """ Merge Voronoi cell into one point """

        first_key = list(self.points.keys())[0]
        demention = len(self.points[first_key][0].coords)

        for key in self.points.keys():

            if key == 'momentum':
                self.points['momentum'] = merge_momentum(demention, self.points)
            if key == 'position':
                self.points['position'] = merge_coordinates(demention, self.points)


 #       """Get one point out of all points in the cell, return Point"""


def merge_coordinates(demention, merged_points):
    """ Merge coordinates into one point """

    size = len(merged_points['position'])

    weights = 0.
    full_values = []
    for i in range(0, demention):
        value = 0
        for point in merged_points['position']:
            value += point.coords[i]
        full_values.append(value/size)

    for point in merged_points['position']:
        weights += point.weight

    result = []
    result.append(Point(full_values, weights))
    return result


def merge_momentum(demention, merged_points):
    """ Merge momentum into one point """

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
    """ Find max coefficient of variance """

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
    """
    Merging algorithm:
    points -- original points
    parametes -- input parameters for Voronoi algorithm(tolerances)

    """

    initial_cell = _VoronoiCell(points)

    result = []
    cells = [initial_cell]

    while len(cells) > 0:
        cell = cells[0]

        max_idx, max_key, max_avg = cell.get_coeff_var()
        needs_subdivision = check_needs_subdivision(parameters, max_avg, max_key)

        if needs_subdivision:

            first_part_cell, secound_part_cell = cell.divide(max_idx, max_key)
            cells.append(secound_part_cell)
            cells.append(first_part_cell)
        else:
            cell.merge()
            new_particle = copy.deepcopy(cell)
            result.append(new_particle)

        cells.remove(cells[0])

    return result


def check_needs_subdivision(parameters, max_avg, max_key):
    """
    Check that Voronoi cell need to devide
    parametrs -- subdivision tolerances
    max_avg -- value of max variance
    max_key -- parameter with max variance

    """

    position_tolerance = parameters.tolerance[0]
    momentum_tolerance = parameters.tolerance[1]

    if max_key == 'position':
        return max_avg > position_tolerance

    if max_key == 'momentum':
        return max_avg > momentum_tolerance


def run_algorithm(points, tolerances):
    """
    Run main algorithm
    points -- start points
    tolerances -- parameters for algorithm

    """

    parameters = VoronoiMergingAlgorithmParameters(tolerances)
    algorithm = VoronoiMergingAlgorithm(parameters)
    result = algorithm.run(points)

    return result
