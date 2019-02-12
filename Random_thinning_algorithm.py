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


def delete_elements_with_null_weight(points):

    result_points = [i for i in points if points[i]. weight != 0]

    return result_points


def count_euclidean_distance(point_coordinates, point_momentum):

    sum_coords = point_coordinates.coords[0] * point_coordinates.coords[0]
    sum_coords += point_coordinates.coords[1] * point_coordinates.coords[1]
    sum_coords += point_coordinates.coords[2] * point_coordinates.coords[2]

    sum_momentum = point_momentum.coords[0] * point_momentum.coords[0]
    sum_momentum += point_momentum.coords[1] * point_momentum.coords[1]
    sum_momentum += point_momentum.coords[2] * point_momentum.coords[2]

    return math.sqrt(sum_coords + sum_momentum)


def weight_distribution(idx_point,  positions, momentum, ranges_patches):

    weight = positions[idx_point].weight
    right_idx = bisect.bisect_left(ranges_patches, idx_point)
    left_idx = right_idx - 1
    size_of_patch = ranges_patches[right_idx] - ranges_patches[left_idx]
    adding_weight = weight/size_of_patch

    for i in range(int(ranges_patches[left_idx]), int(ranges_patches[right_idx])):
        momentum[i].weight += adding_weight
        positions[i].weight += adding_weight


def patches_recount(reduced_points, ranges_patches, parameters):

    reduced_points.sort()
    count_reduced_points = [0] * len(ranges_patches)

    for point in reduced_points:
        patch_idx = bisect.bisect_left(ranges_patches, point) - 1
        count_reduced_points[patch_idx] += 1


    recount_numParticles = []

    for i in range(0, len(parameters.numParticlesOffset)):
        recount_numParticles.append(int(parameters.numParticles[i] - count_reduced_points[i]))

    recount_numParticlesOffset = numpy.cumsum(recount_numParticles, dtype=int)
    recount_numParticlesOffset = numpy.insert(recount_numParticlesOffset, 0, 0)
    recount_numParticlesOffset = numpy.delete(recount_numParticlesOffset, -1)

    return recount_numParticlesOffset, recount_numParticles



def _thinning_ver20():
    data = [[1., 2., 3.], [4., 5., 6.], [7., 1., 5.], [8., 1., 3], [0., 3., 2]]
    weight = [10., 10., 20., 15., 15.]




  #  return weights


def _thinning_ver20(weights, parameters):

    weights = numpy.array(weights)

    ranges_patches = numpy.append(parameters.numParticlesOffset, len(weights))
    for i in range(0, len(ranges_patches) - 1):
        iterate_patch(int(ranges_patches[i]), int(ranges_patches[i + 1]), weights, parameters.reduction_percent)

    return weights


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
