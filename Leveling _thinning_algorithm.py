
def get_indices_to_remove(weight_value_of_reduced_particle, weigths):

    result = []

    for i in range(0, len(weigths)):
        if weigths[i] < weight_value_of_reduced_particle:
            result.append(i)

    return result