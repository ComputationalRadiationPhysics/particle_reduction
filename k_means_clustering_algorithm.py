from sklearn.cluster import KMeans
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


class K_means_clustering_algorithm_Parameters:

    """ Parametrs of k means clustering algorithm
        reduction_percent -- percent of reduced particles
    """

    def __init__(self, reduction_percent):

        self.reduction_percent = reduction_percent


class K_means_clustering_algorithm:
    """ k means reduction algorithm
        iterative optimization procedure for controlling particle populations in particle-in-
        cell (PIC) codes via merging and splitting of computational macro-particles
        based on article https://arxiv.org/abs/1504.03849
    """

    def __init__(self, parameters):
        self.parameters = parameters

    def _run(self, data, weigths):

        dimension = len(data[0])
        size = len(data)
        data = numpy.array(data)
        weights = numpy.array(weigths)
        num_to_remove = int(size * self.parameters.reduction_percent)
        num_to_keep = size - num_to_remove
        kmeans = KMeans(n_clusters=num_to_keep, random_state=0, max_iter=30, tol=0.1).fit(data,
                                                                                          sample_weight=weights)
        result_data, result_weights = recount_data(dimension, num_to_keep, kmeans.labels_, data, weights)

        return result_data, result_weights

