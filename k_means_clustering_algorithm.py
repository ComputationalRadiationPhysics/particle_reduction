import numpy


def merge_points(dimension, vector, weights_vector):
    """ Merge coordinates into one point """

    values_vector = []
    sum_weights = numpy.sum(weights_vector, dtype=float)
    vector = numpy.array(vector)

    for i in range(0, dimension):
        values_vector.append(numpy.average(vector[:, i], weights=weights_vector))

    return values_vector, sum_weights


def recount_data(dimension, num_to_keep, labels, data, weights):
    """ Recalculation of the initial values: if there is one vector in the cell,
    we leave it unchanged, if there are several, we replace it with the average vector.
    dimension -- dimension of source data
    num_to_keep --
    labels -- array indexes in kmeans
    data -- source data
    weights -- source weights

    """

    get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
    result_data = []
    result_weights = []
    for i in range(0, num_to_keep):
        indexes = get_indexes(i, labels)
        selected_weights = weights[indexes]
        selected_data = data[indexes]

        if len(indexes) > 1:
            merged_data, merged_weight = merge_points(dimension, selected_data, selected_weights)
            result_data.append(merged_data)
            result_weights.append(merged_weight)
        else:
            result_data.append(selected_data)
            result_weights.append(selected_weights)
    return result_data, result_weights

