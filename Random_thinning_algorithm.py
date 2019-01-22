import copy


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

