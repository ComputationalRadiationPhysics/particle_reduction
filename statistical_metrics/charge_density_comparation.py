
import numpy
import math
class get_fields():
    """

    Collect particles groups from hdf file
    particles_name -- name of main partcles group

    """

    def __init__(self):
        self.fields = []

    def __call__(self, name, node):
        if node.name.endswith('fields'):
            self.fields = node
        return None

