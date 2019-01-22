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

