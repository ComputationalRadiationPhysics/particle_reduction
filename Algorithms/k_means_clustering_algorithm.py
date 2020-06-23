from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
import numpy
from Algorithms import K_means_divisions



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
        if len(indexes) == 0:
            continue

        selected_weights = weights[indexes]
        selected_data = data[indexes]

        if len(indexes) > 1:
            merged_data, merged_weight = merge_points(dimension, selected_data, selected_weights)
            result_data.append(merged_data)
            result_weights.append(merged_weight)
        elif len(indexes) > 0:
            result_data.append(selected_data[0])
            result_weights.append(selected_weights[0])
    return result_data, result_weights


class K_means_clustering_algorithm_parameters:

    """ Parametrs of k means clustering algorithm
        reduction_percent -- percent of reduced particles
    """

    def __init__(self, reduction_percent, max_iterations=30, tolerance=0.1):
        self.reduction_percent = reduction_percent
        self.max_iterations = max_iterations
        self.tolerance = tolerance


class K_means_clustering_algorithm:
    """ k means reduction algorithm
        iterative optimization procedure for controlling particle populations in particle-in-
        cell (PIC) codes via merging and splitting of computational macro-particles
        based on article https://arxiv.org/abs/1504.03849
    """

    def __init__(self, reduction_percent, max_iterations, tolerance, divisions):
        self.reduction_percent = reduction_percent
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.dimensions = None
        self.divisions = divisions
        self.min_max_values = []

    def _run(self, data, weights,dict_data_indexes):

        range_position = [dict_data_indexes["position"][0], dict_data_indexes["position"][1]]

        self.dict_data_indexes = dict_data_indexes
        if range_position[1] - range_position[1] == 2:
            self.divisions = self.divisions[0:2]

        if len(data) == 0:
            return [], []
        data = numpy.array(data)

        coordinates = data[:, range_position[0]:range_position[1]]
        num_particles_offset, num_particles, moved_values, moved_weights\
            = K_means_divisions.handle_particle_data(coordinates, data, self.divisions, weights)

        result_data = []
        result_weights = []

        for i in range(0, len(num_particles_offset) - 1):
            current_data, current_weights = self.iterate_cell(moved_values[num_particles_offset[i]:num_particles_offset[i + 1]],
                              moved_weights[num_particles_offset[i]:num_particles_offset[i + 1]])

            for point in current_data:
                result_data.append(point)

            for weight in current_weights:
                result_weights.append(weight)

        result_data = numpy.asarray(result_data)
        result_weights = numpy.asarray(result_weights)

        return result_data, result_weights

    def iterate_cell(self, data, weights):
        if len(data) == 0:
            return [], []

        dimension = len(data[0])
        size = len(data)
        data = numpy.array(data)
        weights = numpy.array(weights)
        num_to_remove = int(size * self.reduction_percent)
        num_to_keep = size - num_to_remove

        data_copy = self.normalize_array(data)

      #  kmeans = KMeans(n_clusters=num_to_keep, random_state=0, max_iter=self.max_iterations, tol=self.tolerance).fit(
       #     data,
        #    sample_weight=weights)

        kmeans = MiniBatchKMeans(n_clusters=num_to_keep, max_iter=self.max_iterations,
                                 batch_size=num_to_keep * 3, tol=self.tolerance)\
            .fit(data_copy)


        result_data, result_weights = recount_data(dimension, num_to_keep, kmeans.labels_, data, weights)

        return result_data, result_weights

    def normalize_array(self, data):

        range_momentum = range(self.dict_data_indexes["momentum"][0], self.dict_data_indexes["momentum"][1])
        range_position = range(self.dict_data_indexes["position"][0], self.dict_data_indexes["position"][1])

        normalized_values = []
        dimensions = len(data[0])
        for i in range(0, dimensions):
            current_values = data[:, i]
            if i in range_momentum or i in range_position:
                current_values = normalize_values(current_values)
                normalized_values.append(current_values)



        normalized_values = numpy.transpose(normalized_values)
        return normalized_values


def normalize_values(vector):

    normalize = []
    max_value = max(vector)
    min_value = min(vector)

    diff = max_value - min_value
    if diff == 0:
        normalize = numpy.zeros(shape=(1, len(vector)))[0]
    else:
        normalize = (vector - min_value) / diff

    return normalize


