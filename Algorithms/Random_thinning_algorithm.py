import random
import numpy


class Random_thinning_algorithm_parameters:
    def __init__(self, reduction_percent):
        """..."""
        self.reduction_percent = reduction_percent


class Random_thinning_algorithm:

    def __init__(self, reduction_percent):
        self.reduction_percent = reduction_percent
        self.dimensions = None

    def _run(self, data, weigths):
        size = len(data)

        data = numpy.array(data)
        weigths = numpy.array(weigths)

        indices_to_remove = get_indices_to_remove(size, self.reduction_percent)
        all_data_indexes = numpy.array(range(size))

        select = numpy.in1d(range(all_data_indexes.shape[0]), indices_to_remove)

        indices_to_keep = all_data_indexes[~select]

        total_removed_weight = numpy.sum(weigths[indices_to_remove])
        empty_array = []
        if len(indices_to_keep) == 0:
            return empty_array, empty_array

        weights_to_keep = weigths[indices_to_keep]
        weight_correction = total_removed_weight / len(weights_to_keep)
        weights_to_keep = weights_to_keep + weight_correction

        return data[indices_to_keep], weights_to_keep


def get_indices_to_remove(size, ratio):

    num_to_remove = int(size * ratio)
    result = random.sample(range(size), num_to_remove)
    result.sort()

    return result






