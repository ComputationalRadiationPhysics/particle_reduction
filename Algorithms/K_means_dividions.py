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
