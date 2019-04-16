import numpy as np

def get_ranges(data):
    ranges = []
    dimension = len(data[0])
    for i in range(0, dimension):
        ranges.append((min(data[:, i]), max(data[:, i])))
    return ranges


def count_points_idx(data, splitting_sizes):

    ranges = get_ranges(data)
    size_array = len(data)

    patch_data = None

    if size_array != 0:
        patch_data = Cells_data(data, splitting_sizes, ranges)

    patches_sizes = patch_data.get_size_split()

    list_number_particles_in_parts, links_to_array = \
        points_to_patches(patch_data)

    amount_particles_in_patches, final_size = count_patches_sizes(size_array, patches_sizes, list_number_particles_in_parts,
                                                  links_to_array)

    return amount_particles_in_patches, final_size, list_number_particles_in_parts


def move_values(data, final_size, new_indexes):

    size = len(data)
    dimension = len(data[0])
    moved_values = []
    for i in range(0, size):
        moved_values.append( np.zeros(dimension))

    for i in range(0, len(final_size) - 1):
        for j in range(int(final_size[i]), int(final_size[i + 1])):
            moved_values[j] = data[int(new_indexes[j])]

    return moved_values


def handle_particle_data(coordinates, data, splitting_sizes):

    new_patches_indexes, final_size, num_particles\
        = count_points_idx(coordinates, splitting_sizes)

    moved_values = move_values(data, final_size, new_patches_indexes)

    return final_size, num_particles, moved_values


class Cells_data():

    def __init__(self, data, splitting_number, ranges):
        self.data = data
        self.splitting_number = splitting_number
        self.ranges = ranges
        self.dimension = len(data[0])

    def get_size_split(self):
        size = 0
        if len(self.splitting_number) == 2:
            size = self.splitting_number[0] * self.splitting_number[1]
        else:
            size = self.splitting_number[0] * self.splitting_number[1] * self.splitting_number[2]
        return size

    def get_array_lenght(self):
        return len(self.data)

    def get_patch_x(self, i):
        return get_positon(self.ranges[0][1], self.ranges[0][0], self.splitting_number[0], self.data[i][0])

    def get_patch_y(self, i):
        return get_positon(self.ranges[1][1], self.ranges[1][0], self.splitting_number[1], self.data[i][1])

    def get_patch_z(self, i):
        return get_positon(self.ranges[2][1], self.ranges[2][0], self.splitting_number[2], self.data[i][2])

    def get_position_idx2d(self, x_patch, y_patch):
        return x_patch * self.splitting_number[1] + y_patch

    def get_position_idx3d(self, x_patch, y_patch, z_patch):
        return (x_patch * self.splitting_number[1] + y_patch) * self.splitting_number[2] + z_patch

    def get_position_idx(self, i):

        if self.dimension == 2:
            x_patch = self.get_patch_x(i)
            y_patch = self.get_patch_y(i)
            particle_idx = self.get_position_idx2d(x_patch, y_patch)
        else:
            x_patch = self.get_patch_x(i)
            y_patch = self.get_patch_y(i)
            z_patch = self.get_patch_z(i)
            particle_idx = self.get_position_idx3d(x_patch, y_patch, z_patch)
        return particle_idx


def get_positon(max_coord, min_coord, separator, x_current):

    lenght = max_coord - min_coord
    return max(0, min(int((x_current - min_coord) * separator / lenght), separator - 1))


def count_cells_sizes(links_to_array, final_size, size_indexes, size_array):

    counter_indexes = np.zeros(size_indexes)
    amount_particles_in_patches = np.zeros(max(size_indexes, size_array))

    for i in range(0, len(links_to_array)):
        xy_idx = links_to_array[i]
        start_size = final_size[xy_idx]
        adding_counter = counter_indexes[xy_idx]
        amount_particles_in_patches[int(start_size + adding_counter)] = i
        counter_indexes[xy_idx] = adding_counter + 1
    return amount_particles_in_patches


def points_to_patches(patch_data):

    list_number_particles_in_parts = np.zeros(patch_data.get_size_split() + 1, dtype=int)
    links_to_array = []
    for i in range(0, patch_data.get_array_lenght()):
        particle_idx = patch_data. get_position_idx(i)
        sum_links = list_number_particles_in_parts[particle_idx]
        list_number_particles_in_parts[particle_idx] = sum_links + 1
        links_to_array.append(particle_idx)
    return list_number_particles_in_parts, links_to_array


def count_patches_sizes(size_array, size_indexes, list_number_particles_in_parts, links_to_array):
    final_size = np.cumsum(list_number_particles_in_parts, dtype=int)
    final_size = np.insert(final_size, 0, 0)

    amount_particles_in_patches = count_cells_sizes(links_to_array, final_size, size_indexes, size_array)
    return amount_particles_in_patches, final_size

