
class Diffusion:
    def __init__(self, molar_mass, Kaw_n, wind_speed):
        '''
        Source: A multimedia activity model for ionizable compounds:
            validation study with 2,4-dichlorophenoxyacetic acid, aniline, and trimethoprim.
        Franco, 2009, SI
        (https://www.ncbi.nlm.nih.gov/pubmed/20821507)

        '''
        self.molar_mass = molar_mass
        self.Kaw_n = Kaw_n
        self.wind_speed = wind_speed # m/s


    # diffusion process between compartments
    # only applied for neural species
    def D_diffu_comp1_comp2(self, P_comp1_comp2, area):
        # m/day * m2 = m3/day
        D_comp1_comp2_val = (P_comp1_comp2 * 24.0) * area

        return D_comp1_comp2_val


    """ Permeability P is a mass transfer coefficient or velocity (m/s) for diffusion process """
    def P_air_n(self):
        # P_air_n is the permeability of the neutral molecule in air, m/h
        # m/h = m/s * kg/mol? - the units not make sense here???
        P_air_n_val = (0.01 * (0.3 + 0.2 * self.wind_speed) * ((18.0/(self.molar_mass*1000))**0.335)) * 3600 # mami, m/h
        # P_air_n_val = ((0.43/(24*3600))/0.00475)*3600
        # print P_air_n_val/3600, 'P_air_n_val'
        return P_air_n_val


    def P_Wair(self):
        # P_water_n is the permeability of the neutral molecule in water
        # assume permeabilities of the neutral and the ionic molecule are equal
        # m/h = m/s * kg/mol? - the units not make sense here???
        P_Wair_val = (0.01 * (0.0004 + 0.00004 * (self.wind_speed**2)) * ((32/(self.molar_mass*1000))**0.25)) * 3600 # mami
        return P_Wair_val


    def P_water_air(self, Fr_n_water):
        P_air_n_val = self.P_air_n()
        P_Wair_val = self.P_Wair()
        P_water_air_val = 1/(Fr_n_water/P_Wair_val + 1/(P_air_n_val * self.Kaw_n))
        return P_water_air_val


    def P_soil_air(self, Fr_n_soil):
        # the addtion of the resistances of the soil boundary layer and the air boundary layer
        P_air_n_val = 0.43/(24*0.00475)
        # P_air_n_val = self.P_air_n() # mami
        P_soilAir_n_val = 0.02 # m/h, mami
        # P_soilAir_n_val = 4.3e-9 * 3600
        P_soilWater_val = 2e-6 # m/h, mami

        try:
            part1 = 1/(P_soilWater_val/Fr_n_soil + P_soilAir_n_val * self.Kaw_n)
            part2 = 1/(P_air_n_val *  self.Kaw_n)
            P_soil_air_val = 1/(part1 + part2)
        except:
            P_soil_air_val = 0

        return P_soil_air_val

    # P_waterSed_val = 2.78e-6 # simplebox, fixed value
    # P_sedWater_val = 2.78d-8 # simplebox, fixed value
    # P_waterSed_val = 0.01, P_sedWater_val = 0.0001, mami
    def P_water_sed(self, P_waterSed_val = 2.78e-6*3600, P_sedWater_val = 2.78e-8*3600):
        # the permeability between water and sediment-pore water
        P_water_sed_val = 1/(1/P_waterSed_val + 1/P_sedWater_val)
        return P_water_sed_val


