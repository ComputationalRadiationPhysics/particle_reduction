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
    def __init__(self, parameters):
        self.parameters = parameters

    def run(self, points):
        """Points is a collection of Point"""
        return _thinning(points, self.parameters)

        result_points = copy.deepcopy(points)

        weights = []
        sum_points_before = 0.
        for point in result_points:
            weights.append(point.weight)
            sum_points_before += point.weight

        print('weights before == ' + str(sum_points_before))

        weights =_thinning_ver20(weights, self.parameters)
        sum_points = 0.
        for point in weights:
            sum_points += point

        print('weights after == ' + str(sum_points))

        print('error == '+ str((sum_points - sum_points_before)/sum_points))

        return delete_elements_with_null_weight(result_points)


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


def _thinning(points, parameters):

    size_of_points_array = int(len(points['position']))
    number_reduction_particles = int(size_of_points_array * parameters.reduction_percent)
    result_points = copy.deepcopy(points)

    reduced_points = random.sample(range(size_of_points_array), k=number_reduction_particles)

    positions = result_points['position']
    momentum = result_points['momentum']
    ranges_patches = parameters.numParticlesOffset
    ranges_patches = numpy.append(ranges_patches, size_of_points_array)

    for i in range(0, len(reduced_points)):
        weight_distribution(reduced_points[i], positions, momentum, ranges_patches)

    num_particles_offset, num_particles = patches_recount(reduced_points, ranges_patches, parameters)

    positions = numpy.delete(positions, reduced_points)
    momentum = numpy.delete(momentum, reduced_points)
    result = {}
    result['position'] = positions
    result['momentum'] = momentum

    return result, num_particles_offset, num_particles


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

