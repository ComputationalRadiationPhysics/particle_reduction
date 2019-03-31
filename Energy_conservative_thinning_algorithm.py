

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

