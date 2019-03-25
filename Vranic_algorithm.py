
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

