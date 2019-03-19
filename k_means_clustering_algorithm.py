import numpy


def merge_points(dimension, vector, weights_vector):
    """ Merge coordinates into one point """

    values_vector = []
    sum_weights = numpy.sum(weights_vector, dtype=float)
    vector = numpy.array(vector)

    for i in range(0, dimension):
        values_vector.append(numpy.average(vector[:, i], weights=weights_vector))

    return values_vector, sum_weights

