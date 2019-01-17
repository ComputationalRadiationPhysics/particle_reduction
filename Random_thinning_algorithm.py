

class Random_thinning_alorithm_parametrs:
    def __init__(self, reduction_percent):
        """..."""
        self.reduction_percent = reduction_percent

class Random_thinning_alorithm:
    def __init__(self, parameters):
        self.parameters = parameters

    def run(self, points):
        """Points is a collection of Point"""
        return _merge(points, self.parameters)

