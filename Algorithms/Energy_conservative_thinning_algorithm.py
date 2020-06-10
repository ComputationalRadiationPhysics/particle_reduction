import numpy
import collections
import math
import time


class Energy_conservative_thinning_algorithm_parameters:
    def __init__(self, reduction_percent):
        """..."""
        self.reduction_percent = reduction_percent


class Energy_conservative_thinning_algorithm:

    def __init__(self, ratio, mass):
        self.ratio = ratio
        self.mass = mass
        self.dimensions = None

    def _run(self, data, weigths, dimensions):

        self.dimensions = dimensions

        size = len(data)
        number_of_k_sample = int((1 - self.ratio) * size)
        data = numpy.array(data)
        end_of_dimensions = len(data[0])
        momentum_values = data[:, self.dimensions.dimension_momentum:end_of_dimensions]

        energy_values = calculate_energy_values(momentum_values, self.mass)
        weigths = numpy.array(weigths)

        sample, sum_weighted_energy = get_random_sample(weigths, energy_values, number_of_k_sample)
        indices_to_remove, indexes_to_keep = get_indices_to_remove(sample, size)

        num_indexes_to_keep = numpy.array(indexes_to_keep)[:, 0]

        weights_to_keep = recount_weights(weigths, sample, number_of_k_sample, indexes_to_keep, energy_values, sum_weighted_energy)
        end_time = time.time()


        result_data = numpy.asarray(data[num_indexes_to_keep])
        result_weights = numpy.asarray(weights_to_keep)

        return result_data, result_weights


def calculate_energy_from_momentum(momentum, mass):

    c = 299792458.
    len_momentum_vector = 0.
    if len(momentum) == 3:
        len_momentum_vector = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]
    elif len(momentum) == 2:
        len_momentum_vector = momentum[0] * momentum[0] + momentum[1] * momentum[1]
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
        new_weights = indexes_to_keep[i][1] * weighted_sum_energy/(number_of_k_sample * weights[indexes_to_keep[i][0]] * energy[indexes_to_keep[i][0]])
        result_weights.append(new_weights)
    return result_weights


def get_indices_to_remove(sample, size):

    indexes = numpy.array(list(range(size)))
    indexes_to_keep = [(item, count) for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)
    indices_to_remove = select[~select]
    return indices_to_remove, indexes_to_keep

