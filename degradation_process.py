
class Degradation:

    def __init__(self):
        pass

    def D_deg(self, V, k_deg, Z):
        '''
        Species-specific degradation is not considered here.

        :param V: volumn of the compartment
        :param k_deg: degradation rate in 1/day in the compartment
        :param Z_bulk: Z value for neutral species in the compartment (mol/m3)

        :return:
        '''

        # unit: 1/day * m3 * unitless = m3/day
        D_deg_val = k_deg * V * Z

        return D_deg_val