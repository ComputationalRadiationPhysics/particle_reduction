from sklearn.cluster import KMeans
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
        selected_weights = weights[indexes]
        selected_data = data[indexes]

        if len(indexes) > 1:
            merged_data, merged_weight = merge_points(dimension, selected_data, selected_weights)
            result_data.append(merged_data)
            result_weights.append(merged_weight)
        else:
            result_data.append(selected_data[0])
            result_weights.append(selected_weights[0])
    return result_data, result_weights


class K_means_merge_average_algorithm_parameters:

    """ Parametrs of k means clustering algorithm
        reduction_percent -- percent of reduced particles
    """

    def __init__(self, reduction_percent, max_iterations=30, tolerance=0.1):
        self.reduction_percent = reduction_percent
        self.max_iterations = max_iterations
        self.tolerance = tolerance


class K_means_merge_average_algorithm:
    """
    """

    def __init__(self, reduction_percent, max_iterations, tolerance, divisions):
        self.reduction_percent = reduction_percent
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.dimensions = None
        self.divisions = divisions

    def _run(self, data, weights):

        if len(data) == 0:
            return [], []
        data = numpy.array(data)

        coordinates = data[:, 0:self.dimensions.dimension_position]
        num_particles_offset, num_particles, moved_values, moved_weights \
            = K_means_divisions.handle_particle_data(coordinates, data, self.divisions, weights)

        result_data = []
        result_weights = []

        for i in range(0, len(num_particles_offset) - 1):
            current_data, current_weights = self.iterate_cell(
                moved_values[num_particles_offset[i]:num_particles_offset[i + 1]],
                moved_weights[num_particles_offset[i]:num_particles_offset[i + 1]])

            for point in current_data:
                result_data.append(point)

            for weight in current_weights:
                result_weights.append(weight)

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
        kmeans = KMeans(n_clusters=num_to_keep, random_state=0, max_iter=self.max_iterations, tol=self.tolerance).fit(
            data,
            sample_weight=weights)
        result_data, result_weights = recount_data(dimension, num_to_keep, kmeans.labels_, data, weights)

        return result_data, result_weights
