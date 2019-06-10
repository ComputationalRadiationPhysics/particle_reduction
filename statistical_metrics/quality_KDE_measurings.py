import sys
import numpy
import math


def norm_array01(vector, max_value, min_value):
    array_norm = []
    for point in vector:
        value = (point - min_value) /(max_value - min_value)
        array_norm.append(value)
    return array_norm



def compute_kernel_2d(array_x, array_y, weights_values):

    copy_weights = copy.deepcopy(weights_values)

    x_min = numpy.min(array_x)
    x_max = numpy.max(array_x)
    array_x_norm = norm_array01(array_x, x_min, x_max)

    y_min = numpy.min(array_y)
    y_max = numpy.max(array_y)
    array_y_norm = norm_array01(array_y, y_min, y_max)

    x_min = 0.
    x_max = 1.

    y_min = 0.
    y_max = 1.

    X, Y = numpy.mgrid[x_min:x_max:50j, y_min:y_max:50j]
    lover_values = [x_min, y_min]
    upper_values = [x_max, y_max]

    positions = numpy.vstack([X.ravel(), Y.ravel()])
    values = numpy.vstack([array_x_norm, array_y_norm])
    kernel = stats.gaussian_kde(values, weights=copy_weights)
    value_intergrate = kernel.integrate_box(lover_values, upper_values)
    array_kernell = kernel(positions)

    return array_kernell
