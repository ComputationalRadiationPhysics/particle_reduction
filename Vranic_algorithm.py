import numpy as np

class Segment:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def is_in_segment(self, value):
        if value >=self.start and value <= self.end:
            return True
        else:
            return False


class Momentum_cell:

    def __init__(self):
        self.segment_x = 0.
        self.segment_y = 0.
        self.segment_z = 0.

        self.momentums = []

    def get_idixes(self):
        self.momentums = np.array(self.momentums)
        indexes = np.array(self.momentums[:,3], dtype=('int'))
        return indexes


def create_sorted_momentum_list(data):
    
    size = len(data)
    indexes = [i for i in range(len(data))]
    result = []
    for i in range(0, size):
        result.append([data[i][3], data[i][4], data[i][5], indexes[i]])

    sort_list = sorted(result, key=lambda e: (e[0], e[1], e[2]))

    return sort_list


def get_segment(momentum_min, idx, step):
    start = momentum_min + (idx - 1) * step
    end = momentum_min + (idx) * step
    segment = Segment(start, end)
    return segment


def get_position_idx3d(x_patch, y_patch, z_patch, size_y, size_z):
    return (x_patch * size_y + y_patch) * size_z + z_patch


def get_cell_idx(max_coord, min_coord, separator, x_current):
    """ Get name of particles group """
    lenght = max_coord - min_coord
    return max(0, min(int((x_current - min_coord) * separator / lenght), separator - 1))


def create_momentum_cells(sorted_momentum_list, tolerance_momentum, segment_x, segment_y, segment_z):

    x_momentum_min = segment_x.start
    x_momentum_max = segment_x.end

    y_momentum_min = segment_y.start
    y_momentum_max = segment_y.end

    z_momentum_min = segment_z.start
    z_momentum_max = segment_z.end

    num_x_steps = int((x_momentum_max - x_momentum_min)/ tolerance_momentum) + 1
    num_y_steps = int((y_momentum_max - y_momentum_min) / tolerance_momentum) + 1
    num_z_steps = int((z_momentum_max - z_momentum_min) / tolerance_momentum) + 1

    size_cells = num_x_steps * num_y_steps * num_z_steps
    cells = [Momentum_cell() for i in range(size_cells)]

    for i in range(0, len(sorted_momentum_list)):
        x_idx = get_cell_idx(x_momentum_max, x_momentum_min, num_x_steps, sorted_momentum_list[i][0])

        y_idx = get_cell_idx(y_momentum_max, y_momentum_min, num_y_steps, sorted_momentum_list[i][1])
        z_idx = get_cell_idx(z_momentum_max, z_momentum_min, num_z_steps, sorted_momentum_list[i][2])

        current_idx = get_position_idx3d(x_idx, y_idx, z_idx, num_y_steps, num_z_steps)
        cells[current_idx].momentums.append(sorted_momentum_list[i])

    return cells


def get_start_ranges(data):

    x_min = min(data[:,3])
    x_max = max(data[:,3])

    x_segment = Segment(x_min, x_max)

    y_min = min(data[:,4])
    y_max = max(data[:,4])

    y_segment = Segment(y_min, y_max)

    z_min = min(data[:,5])
    z_max = max(data[:,5])

    z_segment = Segment(z_min, z_max)

    return x_segment, y_segment, z_segment

