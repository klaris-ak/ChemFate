from __future__ import division


##########################################################################
#
#   CLiCC Organic Chemical F&T Model Developed by
#   Mengya Tao, Rucha Thakar, Dillon Elsbury, Kendra Garner, Arturo Keller
#
#   Updated August 31, 2017
#
##########################################################################


class zValue:
    # unit of z value = mol/(Pa-m^3)
    def __init__(self, tempK, Kaw, Kp, aerP, Koc):
        self.tempK = tempK
        self.Kaw = Kaw
        self.Kp = Kp
        self.aerP = aerP
        self.Koc = Koc

    def zAirSub (self):
        R = 8.3144598 # gas constant in Pa-m^3/mol-K
        zAirSub = 1 / (self.tempK * R)
        return zAirSub


    def zAerSub (self, zAirSub):
        # Kp - aerosol-air partition coefficient in m^3-air/ug-aer
        # aerP - aerosol density in kg-aer/m^3-aer
        # unit: mol/m^3-Pa * m^3-air/ug-aer * 10^9 ug/kg * kg-aer/m^3-aer = mol/m^3-Pa
        zAerSub = zAirSub * self.Kp * (10**9) * self.aerP
        # zAerSub = zAirSub * self.Kp
        return zAerSub


    def zWaterSub (self, zAirSub):
        zWaterSub = zAirSub/ self.Kaw
        return zWaterSub


    def zWaterSusSedSub (self, zWaterSub, Kssw_unitless):
        zWaterSusSedSub = zWaterSub * Kssw_unitless
        return zWaterSusSedSub

    def zWaterSedSolidSub (self, zWaterSub, Kbsw_unitless):
        # dWaterSedS - density of freshwater/seawater sediment solid, kg/m^3
        # sedWaterOC - freshwater/seawater sediment organic carbon content
        zWaterSedSolidSub = zWaterSub * Kbsw_unitless
        return zWaterSedSolidSub

    def zSoilSolidSub (self, zWaterSub, Kd_unitless):
        zSoilSolidSub = zWaterSub * Kd_unitless
        return zSoilSolidSub

    def zDeepS (self, zWaterSub, Kd_d_unitless):
        # dsoilP - soil density in kg/m^3
        zDeepS = zWaterSub * Kd_d_unitless
        return zDeepS

    def zAirBulk (self, aerVf, zAirSub, zAerSub):
        zAirBulk = aerVf * zAerSub + (1-aerVf) * zAirSub
        return zAirBulk

    def zWaterBulk (self, susSedVf, zWaterSusSedSub, zWaterSub):
        zWaterBulk = susSedVf * zWaterSusSedSub + (1-susSedVf) * zWaterSub
        return zWaterBulk

    def zSedimentBulk (self, sedSolidVf, zWaterSub, zWaterSedSolidSub):
        zSedimentBulk = sedSolidVf * zWaterSedSolidSub + (1-sedSolidVf) * zWaterSub
        return zSedimentBulk

    def zSoilBulk (self, soilAirVf, soilWaterVf, zAirSub, zWaterSub, zSoilSolidSub):
        zSoilBulk = soilAirVf * zAirSub + soilWaterVf * zWaterSub + (1-soilAirVf-soilWaterVf) * zSoilSolidSub
        return zSoilBulk
