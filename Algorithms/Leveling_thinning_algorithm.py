import numpy


class Leveling_thinning_algorithm_parameters:
    def __init__(self, awg_weight_coef):
        """..."""
        self.awg_weight_coef = awg_weight_coef


class Leveling_thinning_algorithm:

    def __init__(self, awg_weight_coef):
        self.awg_weight_coef = awg_weight_coef
        self.dimensions = None

    def _run(self, data, weigths, dimensions):
        size = len(data)

        data = numpy.array(data)
        weigths = numpy.array(weigths)
        avg_weight = sum(weigths)/size
        weight_value_of_reduced_particle = self.awg_weight_coef * avg_weight

        indices_to_remove = get_indices_to_remove(weight_value_of_reduced_particle, weigths)
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


def get_indices_to_remove(weight_value_of_reduced_particle, weigths):

    result = []

    for i in range(0, len(weigths)):
        if weigths[i] < weight_value_of_reduced_particle:
            result.append(i)

    return result

