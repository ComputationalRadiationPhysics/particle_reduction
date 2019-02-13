import random
import math
import bisect
import copy
import numpy


class RandomThinningAlgorithmParameters:
    def __init__(self, reduction_percent, numParticles, numParticlesOffset):
        """..."""
        self.reduction_percent = reduction_percent
        self.numParticles = numParticles
        self.numParticlesOffset = numParticlesOffset


class RandomThinningAlgorithm:

    def __init__(self, ratio):
        self.ratio = ratio

    def _run(self, data, weigths):
        size = len(data)

        data = numpy.array(data)
        weigths = numpy.array(weigths)

        indices_to_remove = get_indices_to_remove(size, self.ratio)
        all_data_indexes = numpy.array(range(size))

        select = numpy.in1d(range(all_data_indexes.shape[0]), indices_to_remove)

        indices_to_keep = all_data_indexes[~select]
        total_removed_weight = numpy.sum(weigths[indices_to_remove])

        weights_to_keep = weigths[indices_to_keep]
        weight_correction = total_removed_weight / len(weights_to_keep)
        weights_to_keep = weights_to_keep + weight_correction

        return data[indices_to_keep], weights_to_keep


def count_euclidean_distance(point_coordinates, point_momentum):

    sum_coords = point_coordinates.coords[0] * point_coordinates.coords[0]
    sum_coords += point_coordinates.coords[1] * point_coordinates.coords[1]
    sum_coords += point_coordinates.coords[2] * point_coordinates.coords[2]

    sum_momentum = point_momentum.coords[0] * point_momentum.coords[0]
    sum_momentum += point_momentum.coords[1] * point_momentum.coords[1]
    sum_momentum += point_momentum.coords[2] * point_momentum.coords[2]

    return math.sqrt(sum_coords + sum_momentum)


def get_indices_to_remove(size, ratio):

    num_to_remove = int(size * ratio)
    result = random.sample(range(size), num_to_remove)
    result.sort()

    return result






