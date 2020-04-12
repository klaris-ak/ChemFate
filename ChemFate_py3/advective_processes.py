import numpy as np

"""
Advective processes include:
    - Advection
    - Deposition (wet and dry, for aerosol and suspended sediment - sedimentation)
    - Resuspension
    - Burial
    - Runoff
    - Leaching
    - Infiltration
    - Erosion
    - Wind Erosion
    - Aerosol Resuspension (Marine)

"""


class AdvectiveProcess:

    def __init__(self):
        pass

    ''' Advection '''

    def G_advec_air(self, wind_speed, area_air, height_air):
        '''
        :param wind_speed: (m/day)
        :param area_air: area of air (m2)
        :param height_air: height of air (m)
        :return:
        '''
        # unit: m/day * m * m = m3/day
        G_air = wind_speed * np.sqrt(area_air) * height_air  # m3/day

        return G_air

    def G_advec_aer(self, wind_speed, area_air, height_air, aerVf):
        '''
        :param wind_speed: (m/day)
        :param area_air: area of air (m2)
        :param height_air: height of air (m)
        :return:
        '''
        # unit: m/day * m * m = m3/day
        G_air = wind_speed * np.sqrt(area_air) * height_air * aerVf  # m3/day

        return G_air

    def G_advec_susSed(self, waterflow, susSed_Vf):
        '''
        :param wind_speed: (m/day)
        :param area_air: area of air (m2)
        :param height_air: height of air (m)
        :return:
        '''
        # unit: m/day * m * m = m3/day
        G_susSed = waterflow * susSed_Vf  # m3/day

        return G_susSed

    def D_advec_air(self, wind_speed, area_air, height_air, Z_air):
        '''
        :param wind_speed: (m/day)
        :param area_air: area of air (m2)
        :param height_air: height of air (m)
        :return:
        '''
        # unit: m/day * m * m * mol/m3-Pa = mol/Pa-day
        D_air = wind_speed * np.sqrt(area_air) * height_air * Z_air

        return D_air

    def D_advec_aer(self, wind_speed, area_air, height_air, aer_Vf, Z_aer):
        '''
        :param wind_speed: (m/day)
        :param area_air: area of air (m2)
        :param height_air: height of air (m)
        :return:
        '''
        # unit: m/day * m * m * mol/m3-Pa = mol/Pa-day
        D_aer = wind_speed * np.sqrt(area_air) * height_air * aer_Vf * Z_aer  # m3/day

        return D_aer

    def D_advec_water(self, waterflow, Z_water):
        '''
        :param waterflow: (m3/day)
        :return:
        '''
        # unit: m/day * m * m * mol/m3-Pa = mol/Pa-day
        D_fw = waterflow * Z_water

        return D_fw

    def D_advec_susSed(self, waterflow, susSed_Vf, Z_susSed):
        '''
        :param waterflow: (m3/day)
        :return:
        '''
        # unit: m/day * m * m * mol/m3-Pa = mol/Pa-day
        D_fw = waterflow * susSed_Vf * Z_susSed  # m3/day

        return D_fw

    """ Deposition """

    def k_dep_dry(self, particleP, envP, viscosity, radius):
        '''
        Stoke's Law is used to calculate deposition rate of particles
        :param particleP: particle density, unit: kg/m3
        :param envP: environment density, such as airP, unit kg/m3
        :param viscosity: dynamic viscosity of air or water in (kg/m*s)
        :param radius: particle radius env['radiusParticlesAer'] or 'radiusParticleWater' in m

        :return:
        '''

        g = 9.8  # gravitation acceleration in 9.8 m/s^2
        # unit: (kg/m^3 * m/s^2 * m2^)/(kg/(m*s)) = m/s
        k_dep_second = (particleP - envP) * g * (radius ** 2.0) * 4.0 / (18.0 * viscosity)
        # second to day conversion (m/s) to (m/day)
        k_dep_day = k_dep_second * 86400.0  # m/day

        return k_dep_day

    """ air or aerosol wet deposition """

    def D_dep_wet(self, k_dep_wet, scavenging, area, Z):
        # unit: m/day * unitless * m2 * unitless = m3/day
        D_rain_val = k_dep_wet * scavenging * area * Z
        return D_rain_val

    """ air or aerosol rain dissolution """

    def D_rain_diss(self, k_dep_wet, area, Z):
        # unit: m/day * unitless * m2 * unitless = m3/day
        D_rain_val = k_dep_wet * area * Z
        return D_rain_val

    """ air or suspended sediment dry deposition """

    def D_dep_dry(self, k_dep_dry, area, Z):
        # area need to multiply the volumn fraction to indicate area fraciton of those particles
        # unit:  m/day * m2 * unitless = m3/day
        D_dep_dry_val = k_dep_dry * area * Z

        return D_dep_dry_val

    ''' Sediment Resuspension '''

    def D_sedResusp(self, k_resusp, area, Z):
        # unit: m3/m2-hr * 24hr/day * m2 * unitless = m3/day
        D_resusp_val = (k_resusp * 24.0) * area * Z

        return D_resusp_val

    ''' Sediment Burial '''

    def D_burial(self, k_burial, area, Z):
        # unit: m3/m2-hr * 24hr/day * m2 * unitless = m3/day
        D_burial_val = (k_burial * 24.0) * area * Z

        return D_burial_val

    ''' Aerosol Resuspension (Marine) '''

    def D_aeroResusp(self, windspeed_s, coastalA, enrichFactor, seawD, Z_solid):
        # enrichment factor values (EF) estimated from Eisenreich 1980, Sievering et al., 1980 and Eisenreich 1982
        # transport equation adapted from Eisenreich 1982
        # assumes bubbles can only be formed in the top foot (1 m of the marine environment) 1 m/(depth m)

        bubbleRate = self.form_bubble(windspeed_s)  # m/day
        # unitless * m/day * m2 * (1m/m)  = m3/day
        if seawD == 0:
            D_aeroResusp_val = 0
        else:
            D_aeroResusp_val = enrichFactor * bubbleRate * coastalA * (1.0 / seawD) * Z_solid

        return D_aeroResusp_val


    def form_bubble(self, windspeed_s):
        # droplet volume flux (bubbleRate) taken from Kerman 1986, dependent on
        # windspeed (1*10^-9 cm/s at windspeeds of 6 m/s, and (9*10^-4 cm/s at windspeeds greater than 12 m/s)

        if windspeed_s < 6:  # m/s
            bubbleRate = 0
        elif windspeed_s > 12:  # m/s
            # rate from Wu et al. 1984
            # convert from cm/s to m/day
            bubbleRate = (9 * 10 ** -4) * 86400 / 100.0
        else:
            # rate from blanchard 1963
            # convert from cm/s to m/day
            bubbleRate = (1 * 10 ** -9) * 86400 / 100.0

        return bubbleRate

    ''' Water Runoff '''

    def D_runoff(self, precip_mm, CN, soilA, Z_water):
        # precip in mm/day
        # CN - NRCS curve number, unitless
        # Reference:  # http://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1082989.pdf
        # conversion factor from mm/d to inches/d
        # 1 mm = 0.039 inch
        Qthreshold = 0
        if (precip_mm * 0.0393701) > (0.2 * CN):
            Qthreshold = (((precip_mm * 0.0393701) - 0.2 * CN) ** 2) / (((precip_mm * 0.0393701) + 0.8 * CN))

        # 1 in/hr = 1 cfs/acre
        # ft3/s-acre * M3/ft3 * s/d * acre/m2 * m2 * unitless = m3/day
        D_runoff_val = Qthreshold * (1 / 35.312) * (86400) * (1 / 4046.86) * soilA * Z_water

        return D_runoff_val


    ''' Soil Erosion '''

    def D_erosion(self, precip_mm, slope, Kfact, cropManageFactor, supportFactor, soilA, soilP, Z_soilSolid):
        '''
        :param precip: mm/day
        :param slope: soil slope in percentage %
        :param Kfact: soil erodibility factor in hundreds-feet-ton-in/acre-hr-year
        :param cropManageFactor: unitless
        :param supportFactor: unitless
        :param soilA:
        :param soilP: soil density
        :param Z_soilSolid:
        :return:
        Rainfall erosivity factor
        Reference: http://www.sciencedirect.com.proxy.library.ucsb.edu:2048/science/article/pii/S004896971500011X
        '''

        lenslope = self.isFactor(slope)
        er = 0.29 * (1 - 0.72 * np.exp(-0.05 * precip_mm))  # MJ/mm-ha-day
        EI = er * precip_mm * precip_mm * (1.0 / 24.0)  # (MJ-mm/ha-hr-day)
        # unit: MJ-mm/ha-h-day * ha/2.471 acre * in/25.4 mm * 7375.62 hundred ft-lbs/MJ * tons/2000lbs = hundred ft-ton-in/acre-hr-day
        R_unit = EI * (1 / 2.471) * (1 / 25.4) * (7375.62 / 2000)  # % converts R to SI units

        # unit: hundred ft-ton-in/acre-hr-day to kg/m^2-day
        RUSLE = R_unit * Kfact * lenslope * cropManageFactor * supportFactor * (907.185 / 4046.86)

        # unit: kg/(m^2-day) * m^2 * m^3/kg * unitless = m3/day
        D_soilLoss = RUSLE * soilA / soilP * Z_soilSolid

        return D_soilLoss

    ''' Soil Wind Erosion'''

    def D_windErosion(self, windspeed_s, precip_mm, roughness, Kconstant, airP, soilA, A, TSV, TSVmin, z_wind, percWind,
                      windConstant, percUncovered, percSuspended, soilP, Z_soilSolid):
        '''
        :param windspeed_s:
        :param precip:
        :param roughness: in m
        :param Kconstant: conversion constant between horizontal and vertical fluxes in 1/m
        :param airP: density of air, kg/m3
        :param soilA:
        :param A: Saltation fitting parameter
        :param TSV: threshold wind velocity to cause erosion soil in m/s
        :param TSVmin:
        :param z_wind: height at which wind measurements were taken (in US, typically 1.5 m)
        :param percWind:
        :param windConstant:
        :param percUncovered: percent land uncovered and available for wind erosion
        :param percSuspended: percent particles that remain suspended based on average soil size distribution
        :param soilP:
        :param soilV:
        :return:
        '''
        # annualFlux per land use type is given in kg/m2 year
        # we assume that the starting soil moisture is currently 0
        # we assume that if there was no rain that day erosion can occur.
        # we use the saltation equation, given in Kelly et al. 2004, and the
        # vertical flux conversion to estimate total transport of soil between soil and aerosols
        # wind profile and shear rates are taken from Gillette 1978
        # minimum wind speed needed to cause erosion
        # if there was precipitation recently, then the minimum threshold is
        # increased to 30 m/s windspeed, cause not much wind erosion when soil is saturated

        if precip_mm > 10:
            TSV = TSVmin  # % m/s
        # Uz=wind; % m/s
        # k=0.41; % von Karmen constant, unitless

        # u* is the wind shear velocity in m/s for each day based on the existing wind speed
        # m/s / m/m = m/s
        ustar = (windspeed_s * 0.41) / np.log(z_wind / roughness)
        # g=9.81; acceleration due to gravity in m/s2
        # calculate horizontal Qtot and vertical Fa fluxes

        Fa = 0
        if ustar > TSV:
            # unitless * (kg/m3)/(m/s2) * (m/s * ((m/s)^2 - (m/s)^2)) * 1/m = kg/m2-s
            # kg/m-s * m-1 = kg/m2 s
            Fa = (A * (airP / 9.81) * (ustar * ((ustar ** 2) - (TSV ** 2)))) * Kconstant

        # convert to kg/m2-day
        # if wind constantly throughout the day Fa*86400
        # but it doesn't so Fa*86400*0.3 because wind dies down at night (~10
        # hours) and is not consistant throughout the day (10% of the time)
        Fa = Fa * 86400.0 * percWind * windConstant
        # convert to m3/day
        # kg/m2-day * m2 / (kg/m3) = m3/day
        D_windErosion_val = ((Fa * soilA * percUncovered * percSuspended) / soilP) * Z_soilSolid

        return D_windErosion_val


    ''' Soil Infiltration'''
    # https://swat.tamu.edu/media/99192/swat2009-theory.pdf (percolation calculation)
    def D_infiltra(self, precip_mm, CN, evap_mm, FC, soilWC, soilV, soilA, Z_water):
        # FC = field capacity(m^3/m^3): https://stormwater.pca.state.mn.us/index.php?title=Soil_water_storage_properties
        runoff_mm = 0
        if (precip_mm * 0.0393701) > (0.2 * CN):
            runoff_mm = (((precip_mm * 0.0393701) - 0.2 * CN) ** 2) / (((precip_mm * 0.0393701) + 0.8 * CN))

        # unit: mm/day
        infil_mm = precip_mm - runoff_mm - evap_mm

        if infil_mm <= 0:
            k_infil, D_infiltra_val = 0, 0
        else:
            # unit: unitless * m3 + mm/day * m/1000mm * m2 = m3
            soil_water = soilWC * soilV + infil_mm * 0.001 * soilA

            # unit: unitless * m3 = m3
            if soil_water >= FC * soilV:
                # unit: m3/day
                k_infil = soil_water - FC * soilV
            else:
                k_infil = 0

            # unit: m3/day * unitless = m3/day
            D_infiltra_val = k_infil * Z_water

        return D_infiltra_val, k_infil


    ''' Soil Leaching '''
    # assume leaching rate = infiltration rate
    def D_leach(self, k_infil, Z_water):
        # unit: m3/day * unitless = m3/day
        D_leach_val = k_infil * Z_water

        return D_leach_val


    def isFactor(self, slope):
        lenslope = 0.0
        # calculate lenght of slope
        if slope <= 0.2:
            lenslope = 0.06
        elif slope <= 0.5:
            lenslope = 0.1
        elif slope <= 1:
            lenslope = 0.2
        elif slope <= 2:
            lenslope = 0.47
        elif slope <= 3:
            lenslope = 0.8
        elif slope <= 4:
            lenslope = 1.19
        elif slope <= 5:
            lenslope = 1.63
        elif slope <= 6:
            lenslope = 2.11
        elif slope <= 8:
            lenslope = 3.15
        elif slope <= 10:
            lenslope = 4.56
        elif slope <= 12:
            lenslope = 6.28
        elif slope <= 14:
            lenslope = 8.11
        elif slope <= 16:
            lenslope = 10.02
        elif slope <= 20:
            lenslope = 13.99
        elif slope <= 25:
            lenslope = 19.13
        elif slope <= 30:
            lenslope = 24.31
        elif slope <= 40:
            lenslope = 34.48
        elif slope <= 50:
            lenslope = 44.02
        elif slope <= 60:
            lenslope = 52.7

        return lenslope
