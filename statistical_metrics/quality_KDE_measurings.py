import sys


def norm_array01(vector, max_value, min_value):
    array_norm = []
    for point in vector:
        value = (point - min_value) /(max_value - min_value)
        array_norm.append(value)
    return array_norm

