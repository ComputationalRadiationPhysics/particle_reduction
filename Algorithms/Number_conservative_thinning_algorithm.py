import random
import numpy
import collections
import time

class Number_conservative_thinning_algorithm_parameters:
    def __init__(self, reduction_percent):
        """..."""
        self.reduction_percent = reduction_percent


class Number_conservative_thinning_algorithm:

    def __init__(self, reduction_percent):
        self.reduction_percent = reduction_percent
        self.dimensions = None

    def _run(self, data, weigths, dimensions):

        size = len(data)
        if size == 0:
            return [], []
        number_of_k_sample = int((1 - self.reduction_percent) * size)
        data = numpy.array(data)
        weigths = numpy.array(weigths)
        sample = get_random_sample(weigths, number_of_k_sample)
        indices_to_remove, indexes_to_keep = get_indices_to_remove(sample, size)
        num_indexes_to_keep = numpy.array(indexes_to_keep)[:, 0]
        weights_to_keep = recount_weights(weigths, sample, number_of_k_sample, indexes_to_keep)

        result_data = numpy.asarray(data[num_indexes_to_keep])
        result_weights = numpy.asarray(weights_to_keep)

        return result_data, result_weights


def get_random_sample(weights, number_of_k_sample):

    sum_weights = sum(weights)
    norm_weights = [x / sum_weights for x in weights]
    indexes = list(range(0, len(weights)))
    sample = numpy.random.choice(indexes, size=number_of_k_sample, replace=True, p=norm_weights)
    return sample


def recount_weights(weights, sample, number_of_k_sample, indexes_to_keep):

    result_weights = []
    for i in range(0, len(indexes_to_keep)):
       # size_of_values = sample.count(indexes_to_keep[i])
        new_weights = indexes_to_keep[i][1] * weights[indexes_to_keep[i][0]]/number_of_k_sample
        result_weights.append(new_weights)

    return result_weights


def get_indices_to_remove(sample, size):

    indexes = numpy.array(list(range(0, size)))

    sample.tolist()
    indexes_to_keep = [(item, count) for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)

    indices_to_remove = select[~select]

    return indices_to_remove, indexes_to_keep
