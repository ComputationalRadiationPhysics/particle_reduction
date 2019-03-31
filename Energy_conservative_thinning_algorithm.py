import random
import numpy
import collections
import math



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
    indexes = list[range(len(weights))]
    sample = numpy.random.choice(indexes, size=number_of_k_sample, replace=True, p=norm_weights)
    return sample, sum_weighted_energy


def recount_weights(weights, sample, number_of_k_sample, indexes_to_keep, energy, weighted_sum_energy):

    result_weights = []
    for i in range(0, len(indexes_to_keep)):
        size_of_values = sample.count(indexes_to_keep[i])
        new_weights = size_of_values * weighted_sum_energy/(number_of_k_sample * weights[indexes_to_keep[i]] * energy[indexes_to_keep[i]])
        result_weights.append(new_weights)
    return result_weights


def get_indices_to_remove(sample, size):

    indexes = list[range(len(size))]
    indexes_to_keep = [item for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)
    indices_to_remove = select[~select]
    return indices_to_remove, indexes_to_keep

