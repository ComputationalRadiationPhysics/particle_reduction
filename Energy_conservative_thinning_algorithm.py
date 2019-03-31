import random
import numpy
import collections
import math


class Energy_conservative_thinning_algorithm_parameters:
    def __init__(self, reduction_percent, numParticles, numParticlesOffset):
        """..."""
        self.reduction_percent = reduction_percent
        self.numParticles = numParticles
        self.numParticlesOffset = numParticlesOffset


class Energy_conservative_thinning_algorithm:

    def __init__(self, number_of_k_sample):
        self.number_of_k_sample = number_of_k_sample

    def _run(self, data, weigths, mass):

        size = len(data)
        data = numpy.array(data)
        momentum_values = data[:, 3:6]

        energy_values = calculate_energy_values(momentum_values, mass)
        weigths = numpy.array(weigths)
        sample, sum_weighted_energy = get_random_sample(weigths, energy_values, self.number_of_k_sample)
        indices_to_remove, indexes_to_keep = get_indices_to_remove(sample, size)
        weights_to_keep = recount_weights(weigths, sample, self.number_of_k_sample, indexes_to_keep, energy_values, sum_weighted_energy)
        return data[indexes_to_keep], weights_to_keep


def calculate_energy_from_momentum(momentum, mass):

    c = 299792458.
    len_momentum_vector = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]
    e_ = math.sqrt(len_momentum_vector * c + (mass * mass * c * c) * (mass * mass * c * c))
    return e_


def calculate_energy_values(momentums, mass):

    energy_values = []
    for i in range(0, len(momentums)):
        e = calculate_energy_from_momentum(momentums[i], mass)
        energy_values.append(e)

    return energy_values


def get_random_sample(weights, energy_values, number_of_k_sample):

    weighted_energy = []

    for i in range(0, len(weights)):
        weighted_energy.append(weights[i] * energy_values[i])
    sum_weighted_energy = sum(weighted_energy)
    norm_weights = [x / sum_weighted_energy for x in weighted_energy]
    indexes = list(range(len(weights)))
    sample = numpy.random.choice(indexes, size=number_of_k_sample, replace=True, p=norm_weights)
    return sample, sum_weighted_energy


def recount_weights(weights, sample, number_of_k_sample, indexes_to_keep, energy, weighted_sum_energy):

    sample = sample.tolist()
    result_weights = []
    for i in range(0, len(indexes_to_keep)):
        size_of_values = sample.count(indexes_to_keep[i])
        new_weights = size_of_values * weighted_sum_energy/(number_of_k_sample * weights[indexes_to_keep[i]] * energy[indexes_to_keep[i]])
        result_weights.append(new_weights)
    return result_weights


def get_indices_to_remove(sample, size):

    indexes = numpy.array(list(range(size)))
    indexes_to_keep = [item for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)
    indices_to_remove = select[~select]
    return indices_to_remove, indexes_to_keep

