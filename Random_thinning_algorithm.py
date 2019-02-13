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


def iterate_patch(left_bound, right_bound, weights, reduction_percent):
    print('left_bound '+ str(left_bound))
    print('right_bound  '+ str(right_bound))

    size_of_patch = right_bound - left_bound
    if size_of_patch != 0:
        number_reduction_particles = int(size_of_patch * reduction_percent)
        print('number_reduction_particles  '+ str(number_reduction_particles))

        reduced_points = random.sample(range(left_bound, right_bound), k=number_reduction_particles)
        reduced_points.sort()

        selected_points = weights[reduced_points]

        sum_weights = numpy.sum(selected_points)

        added_weight = sum_weights/(size_of_patch - number_reduction_particles)
        print('added weight  ' + str(added_weight))

      #  print('added_weight' + str(added_weight))
        weights[left_bound:right_bound] = weights[left_bound:right_bound] + added_weight
        print(weights[reduced_points])
        weights[reduced_points] = [0] * len(selected_points)
        print('------------------------------')
        print(weights[reduced_points])


def get_indices_to_remove(size, ratio):

    num_to_remove = int(size * ratio)
    result = random.sample(range(size), num_to_remove)
    result.sort()

    return result
