from __future__ import division

class Diffusion:

    def __init__(self):

        pass


    def diffusion_air_water (self, airMTC, waterMTC, cross_sectional_area, zAirSub, zWaterSub):
        # unit: m/day * m^2 * mol/m^3-Pa = mol/Pa-day
        air_part = airMTC * cross_sectional_area * zAirSub
        water_part = waterMTC * cross_sectional_area * zWaterSub
        try:
            air_water_diffusion = 1 / (1/air_part + 1/water_part)
        except:
            air_water_diffusion = 0
        return air_water_diffusion



    def diffusion_air_soil (self, airSoilMTC, soilAirMTC, soilWaterMTC, cross_sectional_area, zAirSub, zWaterSub):
        # unit: m/day * m^2 * mol/m^3-Pa = mol/Pa-day
        soilAirBound_part = airSoilMTC * cross_sectional_area * zAirSub
        soilAir_part = cross_sectional_area * soilAirMTC * zAirSub
        soilWater_part = cross_sectional_area * soilWaterMTC * zWaterSub
        try:
            air_soil_diffusion = 1/(1/soilAirBound_part + 1/(soilAir_part + soilWater_part))
        except:
            air_soil_diffusion = 0
        return air_soil_diffusion


    def diffusion_sediment_water (self, sedmtWaterMTC, cross_sectional_area, zWaterSub):
        # unit: m/day * m^2 * mol/m^3-Pa = mol/Pa-day
        waterSedmMTC = 0.01 * 24.0 # Mackay P178, m/h
        water_part = waterSedmMTC * cross_sectional_area * zWaterSub
        sediment_part = sedmtWaterMTC * cross_sectional_area * zWaterSub
        try:
            sediment_water_diffusion = 1/(1/water_part + 1/sediment_part)
        except:
            sediment_water_diffusion = 0
        return sediment_water_diffusion



class MTC:

    def __init__(self, molar_volume, tempK, MW):
        self.molar_volume = molar_volume # cm3/mol
        self.tempK = tempK
        self.MW = MW


    def airMD(self):
        # air diffusivity in cm2/second
        # reference temperature 25 C: refer_temp is the room temperature
        # reference: Handbook of Chemical Mass Transport in the Environment, P78
        refer_temp = 298.15  # K
        airMD_value = 2.35 * (self.molar_volume ** (-0.73)) * ((self.tempK / refer_temp)) ** 1.75
        # cm2/second * 8.64 = m2/day
        airMD_value = airMD_value * 8.64
        return airMD_value


    def waterMD(self):
        # water diffusivity in cm2/s
        # reference: Schwarzenbach, Environmental Organic Chemistry, 1993
        waterMD = 8.64 * 0.000274 * (self.MW ** -0.71) # m2/day
        return waterMD


    def airWaterMTC(self, airMD):
        # reference: <chemical mass transport in the environment> P78, equation 5.18)
        # air boundary layer thickness = 5 mm = 0.005 m (Mackay, P197)
        diffusion_pathLength = 0.005
        # m2/day / m = m/day
        airWaterMTCair = airMD / diffusion_pathLength
        return airWaterMTCair


    def waterAirMTC(self, waterMD):
        # reference: <chemical mass transport in the environment> P81, equation 5.24)
        # still water layer thickness = 0.5 mm = 0.0005 m (Mackay, P152)
        diffusion_pathLength = 0.0005
        waterAirMTC = waterMD / diffusion_pathLength
        return waterAirMTC


    def airSoilMTC(self, airMD):
        # reference: <chemical mass transport in the environment> P177)
        diffusion_pathLength = 0.005
        airSoilMTCair = airMD / diffusion_pathLength
        return airSoilMTCair


    def soilAirMTC(self, airMD, soilAVf, soilWVf):
        diffusion_pathLength = 0.025
        # calcuate the effective diffusivity of soil air
        # Millington-Quirk equation (Mackay P197)
        soilAirMDeff = airMD * (soilAVf ** (10/3))/((soilAVf + soilWVf) ** 2)
        soilAirMTCsoil = soilAirMDeff/diffusion_pathLength
        return soilAirMTCsoil


    def soilWaterMTC(self, waterMD, soilAVf, soilWVf):
        diffusion_pathLength = 0.05
        soilWaterMDeff = waterMD * (soilWVf ** (10/3))/((soilAVf + soilWVf) ** 2)
        soilWaterMTCsoil = soilWaterMDeff / diffusion_pathLength
        return soilWaterMTCsoil


    def sedmtWaterMTC(self, waterMD, fsedpercSolid):
        diffusion_pathLength = 0.005 # diffusion path length in sediment P178 Mackay
        waterMDeff = waterMD * ((1-fsedpercSolid)**1.5) # P160 Mackay
        waterParticleMTC = waterMDeff/diffusion_pathLength
        return waterParticleMTC

