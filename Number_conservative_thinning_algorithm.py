




def get_indices_to_remove(sample, size):

    indexes = numpy.array(list(range(0, size)))
    sample.tolist()
    indexes_to_keep = [item for item, count in collections.Counter(sample).items() if count > 1]
    select = numpy.in1d(range(indexes.shape[0]), indexes_to_keep)
    indices_to_remove = select[~select]
    return indices_to_remove, indexes_to_keep
