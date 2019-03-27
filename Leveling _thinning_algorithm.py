

class Leveling_Thinning_Algorithm_Parameters:
    def __init__(self, awg_weight_coef, numParticles, numParticlesOffset):
        """..."""
        self.awg_weight_coef = awg_weight_coef
        self.numParticles = numParticles
        self.numParticlesOffset = numParticlesOffset


def get_indices_to_remove(weight_value_of_reduced_particle, weigths):

    result = []

    for i in range(0, len(weigths)):
        if weigths[i] < weight_value_of_reduced_particle:
            result.append(i)

    return result