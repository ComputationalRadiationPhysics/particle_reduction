import random
import numpy
import collections


class Number_conservative_thinning_algorithm_parameters:
    def __init__(self, reduction_percent):
        """..."""
        self.reduction_percent = reduction_percent


class Number_conservative_thinning_algorithm:

    def __init__(self, reduction_percent):
        self.reduction_percent = reduction_percent
        self.dimensions = None

    def _run(self, data, weigths):

        size = len(data)
        number_of_k_sample = int((1 - self.reduction_percent) * size)
        data = numpy.array(data)
        weigths = numpy.array(weigths)
        sample = get_random_sample(weigths, number_of_k_sample)
        indices_to_remove, indexes_to_keep = get_indices_to_remove(sample, size)
        weights_to_keep = recount_weights(weigths, sample, number_of_k_sample, indexes_to_keep)
        return data[indexes_to_keep], weights_to_keep


def get_random_sample(weights, number_of_k_sample):

    sum_weights = sum(weights)
    norm_weights = [x / sum_weights for x in weights]
    indexes = list(range(0, len(weights)))
    sample = numpy.random.choice(indexes, size=number_of_k_sample, replace=True, p=norm_weights)
    return sample


def recount_weights(weights, sample, number_of_k_sample, indexes_to_keep):

    sample = sample.tolist()
    result_weights = []
    for i in range(0, len(indexes_to_keep)):
        size_of_values = sample.count(indexes_to_keep[i])
        new_weights = size_of_values * weights[indexes_to_keep[i]]/number_of_k_sample
        result_weights.append(new_weights)
    return result_weights


def get_indices_to_remove(sample, size):

    indexes = numpy.array(list(range(0, size)))
    sample.tolist()
    indexes_to_keep = [item for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)
    indices_to_remove = select[~select]
    return indices_to_remove, indexes_to_keep
