from __future__ import division
import xlrd
from collections import OrderedDict
from xlrd import open_workbook
from datetime import datetime
from numpy import *
import pandas as pd


class LoadData:

    def __init__(self, chem_type, chem_file, region_file, release_file, start_date, end_date):
        self.chem_type = chem_type
        self.chem_file = chem_file
        self.release_file = release_file
        self.region_file = region_file
        self.start_date = start_date
        self.end_date = end_date
        self.region_workbook = open_workbook(self.region_file)
        self.release_workbook = open_workbook(self.release_file)


    def load_date(self):
        # load the start rows and end rows from the climate dataset
        original_start_date = datetime.strptime('2005 1 1', "%Y %m %d")
        start_day = datetime.strptime(self.start_date, "%Y %m %d")
        start_row = (start_day - original_start_date).days + 1
        end_day = datetime.strptime(self.end_date, "%Y %m %d")
        sim_days = (end_day - start_day).days
        end_row = start_row + sim_days + 1
        return start_row, end_row


    def get_Koc_acid(self, smiles, cas):
        # if the chemical is organic acid, this parameter would be used to calculate Kd_i in soil
        Koc_acid = None
        df = pd.read_excel('./IonizableChem_DB.xlsx', sheet_name='Koc_organicAcid')
        # check if smiles in the SMILES column
        # if contains, a row of values would return
        # if not contain, an empty dataframe would return
        if smiles is not None:
            df2 = df[df['SMILES'] == smiles]
        else:
            df2 = df[df['SMILES'] == '--']

        if cas is not None:
            df3 = df[df['CAS'].str.contains(cas)]
        else:
            df3 = df[df['CAS'].str.contains('--')]

        if df2['SMILES'].empty:
            # check if cas column contains the cas
            if df3['CAS'].empty:
                return Koc_acid
            else:
                Koc_acid = df3['Koc_i'].values[0]
        else:
            Koc_acid = df2['Koc_i'].values[0]
        return Koc_acid


    def get_infil_rate(self, soil_type, slope):
        infil_rate = None
        slope_col = None
        df = pd.read_excel('./IonizableChem_DB.xlsx', sheetname='infiltrationRate', index_col=0)

        if slope <= 4:
            slope_col = '0-4%'
        elif slope <= 8:
            slope_col = '5-8%'
        elif slope <= 12:
            slope_col = '8-12%'
        elif slope <= 16:
            slope_col = '12-16%'
        else:
            slope_col = 'over 16%'

        infil_rate = df.loc[soil_type, slope_col]
        return infil_rate


    def load_chemParams(self, chem_type, env):
        # load chemical properties
        chem_wb = xlrd.open_workbook(self.chem_file)
        chem_ws = chem_wb.sheet_by_name('Sheet1')

        chem_code = chem_ws.col_values(1, start_rowx=1, end_rowx=None)
        chem_value = chem_ws.col_values(2, start_rowx=1, end_rowx=None)

        chem_params = {}

        for i in range(len(chem_code)):
            chem_params[chem_code[i]] = chem_value[i]

        # get the Koc_acid value
        if chem_type == 'IonizableOrganic':
            chem_params['Koc_acid'] = self.get_Koc_acid(chem_params['smiles'], chem_params['cas'])
        
        if chem_type == 'NonionizableOrganic':
            # soil/water partition coefficient Kd = Koc * foc
            # sorbed concentration (mg/kg) / dissolved concentration (mg/L) = L/kg
            # unit of Koc is equal to the unit of Kd: L/kg, divide by 1000, 1000 L = 1 m^3
            # Kd in m^3-water/kg-soil
            chem_params['Kd1'] = (chem_params['Koc_n'] * env['soilOC1']) / 1000
            chem_params['Kd2'] = (chem_params['Koc_n'] * env['soilOC2']) / 1000
            chem_params['Kd3'] = (chem_params['Koc_n'] * env['soilOC3']) / 1000
            chem_params['Kd4'] = (chem_params['Koc_n'] * env['soilOC4']) / 1000

            chem_params['Kd1_d'] = (chem_params['Koc_n'] * env['dsoilOC1']) / 1000
            chem_params['Kd2_d'] = (chem_params['Koc_n'] * env['dsoilOC1']) / 1000
            chem_params['Kd3_d'] = (chem_params['Koc_n'] * env['dsoilOC1']) / 1000
            chem_params['Kd4_d'] = (chem_params['Koc_n'] * env['dsoilOC1']) / 1000
            # convert Kd to unitless, multiply by soil density
            # m3-water/kg-soil * kg-soil/m3-soil = m3-water/m3-soil
            chem_params['Kd1_unitless'] = chem_params['Kd1'] * env['dSS1']
            chem_params['Kd2_unitless'] = chem_params['Kd2'] * env['dSS2']
            chem_params['Kd3_unitless'] = chem_params['Kd3'] * env['dSS3']
            chem_params['Kd4_unitless'] = chem_params['Kd4'] * env['dSS4']

            chem_params['Kd1_d_unitless'] = chem_params['Kd1_d'] * env['deepsP1']
            chem_params['Kd2_d_unitless'] = chem_params['Kd2_d'] * env['deepsP2']
            chem_params['Kd3_d_unitless'] = chem_params['Kd3_d'] * env['deepsP3']
            chem_params['Kd4_d_unitless'] = chem_params['Kd4_d'] * env['deepsP4']

            # sediment/water partition coefficient Kssw (suspended sediment - water) and Kbsw (bottom sediment - water)
            chem_params['Kssfw'] = (chem_params['Koc_n'] * env['freshssOC']) / 1000
            chem_params['Ksssw'] = (chem_params['Koc_n'] * env['seassOC']) / 1000
            chem_params['Kbsfw'] = (chem_params['Koc_n'] * env['sedFWOC']) / 1000
            chem_params['Kbssw'] = (chem_params['Koc_n'] * env['sedSWOC']) / 1000
            # convert Kss and Kbs to unitless, multiple suspended sediment and sediment's density
            chem_params['Kssfw_unitless'] = chem_params['Kssfw'] * env['freshssP']
            chem_params['Ksssw_unitless'] = chem_params['Ksssw'] * env['seassP']
            chem_params['Kbsfw_unitless'] = chem_params['Kbsfw'] * env['dFWSedS']
            chem_params['Kbssw_unitless'] = chem_params['Kbssw'] * env['dSWSedS']

            # aerosol-air partition coefficient Kp in m^3-air/ug-aer
            # m^3/ug * 10^9 ug/kg = 10^9 m^3/kg
            # convert Kp to unitless, multiply by its density
            # chem_params['Kp_unitless'] = chem_params['Kp_n'] * (10 ** 9) * env['aerP']
            chem_params['Kp_unitless'] = 0.54 * (chem_params['Kow_n']/chem_params['Kaw_n']) * env['aerOC'] * (env['aerP']/1000)
            # air-aerosol partiton coefficient Kairaer
            try:
                chem_params['Kairaer'] = 1 / chem_params['Kp_unitless']
            except:
                chem_params['Kairaer'] = 0

        # degradation rate: k = 0.693/(halflife/24)
        # transform the units from hours to days
        if self.chem_type != 'Metal':
            chem_params['kDeg_air_n'] = 24.0 * log(2.0) / chem_params['HL_air_n']
            chem_params['kDeg_aer_n'] = 24.0 * log(2.0) / chem_params['HL_aer_n']
            chem_params['kDeg_fw_n'] = 24.0 * log(2.0) / chem_params['HL_fWater_n']
            chem_params['kDeg_fSS_n'] = 24.0 * log(2.0) / chem_params['HL_fSS_n']
            chem_params['kDeg_fSedW_n'] = 24.0 * log(2.0) / chem_params['HL_fSedW_n']
            chem_params['kDeg_fSedS_n'] = 24.0 * log(2.0) / chem_params['HL_fSedS_n']
            chem_params['kDeg_sw_n'] = 24.0 * log(2.0) / chem_params['HL_sWater_n']
            chem_params['kDeg_sSS_n'] = 24.0 * log(2.0) / chem_params['HL_sSS_n']
            chem_params['kDeg_sSedW_n'] = 24.0 * log(2.0) / chem_params['HL_sSedW_n']
            chem_params['kDeg_sSedS_n'] = 24.0 * log(2.0) / chem_params['HL_sSedS_n']
            chem_params['kDeg_soilA1_n'] = 24.0 * log(2.0) / chem_params['HL_soilA1_n']
            chem_params['kDeg_soilW1_n'] = 24.0 * log(2.0) / chem_params['HL_soilW1_n']
            chem_params['kDeg_soilS1_n'] = 24.0 * log(2.0) / chem_params['HL_soilS1_n']
            chem_params['kDeg_deepS1_n'] = 24.0 * log(2.0) / chem_params['HL_soilDeep1_n']
            chem_params['kDeg_soilA2_n'] = 24.0 * log(2.0) / chem_params['HL_soilA2_n']
            chem_params['kDeg_soilW2_n'] = 24.0 * log(2.0) / chem_params['HL_soilW2_n']
            chem_params['kDeg_soilS2_n'] = 24.0 * log(2.0) / chem_params['HL_soilS2_n']
            chem_params['kDeg_deepS2_n'] = 24.0 * log(2.0) / chem_params['HL_soilDeep2_n']
            chem_params['kDeg_soilA3_n'] = 24.0 * log(2.0) / chem_params['HL_soilA3_n']
            chem_params['kDeg_soilW3_n'] = 24.0 * log(2.0) / chem_params['HL_soilW3_n']
            chem_params['kDeg_soilS3_n'] = 24.0 * log(2.0) / chem_params['HL_soilS3_n']
            chem_params['kDeg_deepS3_n'] = 24.0 * log(2.0) / chem_params['HL_soilDeep3_n']
            chem_params['kDeg_soilA4_n'] = 24.0 * log(2.0) / chem_params['HL_soilA4_n']
            chem_params['kDeg_soilW4_n'] = 24.0 * log(2.0) / chem_params['HL_soilW4_n']
            chem_params['kDeg_soilS4_n'] = 24.0 * log(2.0) / chem_params['HL_soilS4_n']
            chem_params['kDeg_deepS4_n'] = 24.0 * log(2.0) / chem_params['HL_soilDeep4_n']

        if chem_type == 'IonizableOrganic':
            chem_params['kDeg_aer_i'] = 24.0 * log(2.0) / chem_params['HL_aer_i']
            chem_params['kDeg_fw_i'] = 24.0 * log(2.0) / chem_params['HL_fWater_i']
            chem_params['kDeg_fSS_i'] = 24.0 * log(2.0) / chem_params['HL_fSS_i']
            chem_params['kDeg_fSedW_i'] = 24.0 * log(2.0) / chem_params['HL_fSedW_i']
            chem_params['kDeg_fSedS_i'] = 24.0 * log(2.0) / chem_params['HL_fSedS_i']
            chem_params['kDeg_sw_i'] = 24.0 * log(2.0) / chem_params['HL_sWater_i']
            chem_params['kDeg_sSS_i'] = 24.0 * log(2.0) / chem_params['HL_sSS_i']
            chem_params['kDeg_sSedW_i'] = 24.0 * log(2.0) / chem_params['HL_sSedW_i']
            chem_params['kDeg_sSedS_i'] = 24.0 * log(2.0) / chem_params['HL_sSedS_i']
            chem_params['kDeg_soilA1_i'] = 24.0 * log(2.0) / chem_params['HL_soilA1_i']
            chem_params['kDeg_soilW1_i'] = 24.0 * log(2.0) / chem_params['HL_soilW1_i']
            chem_params['kDeg_soilS1_i'] = 24.0 * log(2.0) / chem_params['HL_soilS1_i']
            chem_params['kDeg_deepS1_i'] = 24.0 * log(2.0) / chem_params['HL_soilDeep1_i']
            chem_params['kDeg_soilA2_i'] = 24.0 * log(2.0) / chem_params['HL_soilA2_i']
            chem_params['kDeg_soilW2_i'] = 24.0 * log(2.0) / chem_params['HL_soilW2_i']
            chem_params['kDeg_soilS2_i'] = 24.0 * log(2.0) / chem_params['HL_soilS2_i']
            chem_params['kDeg_deepS2_i'] = 24.0 * log(2.0) / chem_params['HL_soilDeep2_i']
            chem_params['kDeg_soilA3_i'] = 24.0 * log(2.0) / chem_params['HL_soilA3_i']
            chem_params['kDeg_soilW3_i'] = 24.0 * log(2.0) / chem_params['HL_soilW3_i']
            chem_params['kDeg_soilS3_i'] = 24.0 * log(2.0) / chem_params['HL_soilS3_i']
            chem_params['kDeg_deepS3_i'] = 24.0 * log(2.0) / chem_params['HL_soilDeep3_i']
            chem_params['kDeg_soilA4_i'] = 24.0 * log(2.0) / chem_params['HL_soilA4_i']
            chem_params['kDeg_soilW4_i'] = 24.0 * log(2.0) / chem_params['HL_soilW4_i']
            chem_params['kDeg_soilS4_i'] = 24.0 * log(2.0) / chem_params['HL_soilS4_i']
            chem_params['kDeg_deepS4_i'] = 24.0 * log(2.0) / chem_params['HL_soilDeep4_i']

        if chem_type == 'Metal':
            # assign the fraction value when user enters 0, make it to 1e-20
            for i in range(len(chem_code)):
                if chem_value[i] == 0:
                    chem_value[i] = 1e-20
                chem_params[chem_code[i]] = chem_value[i]

        # g/mol / (g/cm3) = cm3/mol
        chem_params['molar_mass'] = chem_params['MW'] / 1000.0  # kg/mol
        if self.chem_type == 'NonionizableOrganic':
            chem_params['molar_volume'] = chem_params['MW'] / chem_params['MD']  # cm3/mol

        return chem_params


    def load_compart_presence(self):
        # load presence of each compartment
        presence_worksheet = self.region_workbook.sheet_by_name('Presence')
        presence_code = presence_worksheet.col_values(1, start_rowx=1, end_rowx=None)
        presence_value = presence_worksheet.col_values(2, start_rowx=1, end_rowx=None)
        presence_loading = zip(presence_code, presence_value)
        presence = OrderedDict()
        for name, value in presence_loading:
            presence[name] = value

        return presence


    def load_env_params(self):
        # load the environmental parameters
        env_ws = self.region_workbook.sheet_by_name('Environment')
        env_code = env_ws.col_values(1, start_rowx=1, end_rowx=None)
        env_value = env_ws.col_values(2, start_rowx=1, end_rowx=None)
        env_loading = zip(env_code, env_value)
        env = OrderedDict()
        for name, value in env_loading:
            env[name] = value

        # area calculation
        env['area'] = env['freshwA'] + env['seawA'] + env['soilA1'] + env['soilA2'] + env['soilA3'] + env['soilA4']
        env['fwA'] = env['freshwA']
        env['swA'] = env['seawA']
        env['airA'] = env['area']
        env['sedFWA'] = env['freshwA']
        env['sedSWA'] = env['seawA']
        env['deepsA1'] = env['soilA1']
        env['deepsA2'] = env['soilA2']
        env['deepsA3'] = env['soilA3']
        env['deepsA4'] = env['soilA4']
        # volume calculation
        env['areaV'] = env['area'] * env['airH']
        env['fWaterV'] = env['freshwA'] * env['freshwD']
        env['sWaterV'] = env['seawA'] * env['seawD']
        # kg-aer/m3-air * m3-air / (kg-aer/m3-aer) = m3 aer
        if env['aerP'] == 0:
            env['aerV'] = 0
        else:
            env['aerV'] = env['aerC'] * (env['areaV'] / env['aerP'])

        if env['freshssP'] == 0:
            env['fSSV'] = 0
        else:
            env['fSSV'] = env['freshssC'] * (env['fWaterV'] / env['freshssP'])

        if env['seassP'] == 0:
            env['sSSV'] = 0
        else:
            env['sSSV'] = env['seassC'] * (env['sWaterV'] / env['seassP'])

        env['airV'] = env['areaV'] - env['aerV']
        env['fwV'] = env['fWaterV'] - env['fSSV']
        env['swV'] = env['sWaterV'] - env['sSSV']
        # freshwater sediment volume
        env['sedFWV'] = env['sedFWA'] * env['sedFWD']
        env['sedSWV'] = env['sedSWA'] * env['sedSWD']
        env['fSedWV'] = env['sedFWV'] * (1 - env['fsedpercSolid'])
        env['fSedSV'] = env['sedFWV'] * env['fsedpercSolid']
        env['sSedWV'] = env['sedSWV'] * (1 - env['ssedpercSolid'])
        env['sSedSV'] = env['sedSWV'] * env['ssedpercSolid']
        # soil commpartments volume (m^3)
        env['soilV1'] = env['soilA1'] * env['soilD1']
        env['soilV2'] = env['soilA2'] * env['soilD2']
        env['soilV3'] = env['soilA3'] * env['soilD3']
        env['soilV4'] = env['soilA4'] * env['soilD4']
        env['soilAV1'] = env['soilA1'] * env['soilD1'] * env['soilAC1']
        env['soilAV2'] = env['soilA2'] * env['soilD2'] * env['soilAC2']
        env['soilAV3'] = env['soilA3'] * env['soilD3'] * env['soilAC3']
        env['soilAV4'] = env['soilA4'] * env['soilD4'] * env['soilAC4']
        # surface soil water volume
        env['soilWV1'] = env['soilA1'] * env['soilD1'] * env['soilWC1']
        env['soilWV2'] = env['soilA2'] * env['soilD2'] * env['soilWC2']
        env['soilWV3'] = env['soilA3'] * env['soilD3'] * env['soilWC3']
        env['soilWV4'] = env['soilA4'] * env['soilD4'] * env['soilWC4']
        # surface soil solid volume
        env['soilSV1'] = env['soilA1'] * env['soilD1'] * (1 - env['soilWC1'] - env['soilAC1'])
        env['soilSV2'] = env['soilA2'] * env['soilD2'] * (1 - env['soilWC2'] - env['soilAC2'])
        env['soilSV3'] = env['soilA3'] * env['soilD3'] * (1 - env['soilWC3'] - env['soilAC3'])
        env['soilSV4'] = env['soilA4'] * env['soilD4'] * (1 - env['soilWC4'] - env['soilAC4'])
        # env['soilSV2'] = env['soilA2'] * env['soilD2'] * env['soilSC2']
        # env['soilSV1'] = env['soilA1'] * env['soilD1'] * env['soilSC1']
        # env['soilSV3'] = env['soilA3'] * env['soilD3'] * env['soilSC3']
        # env['soilSV4'] = env['soilA4'] * env['soilD4'] * env['soilSC4']

        # deep soil volume (m^3)
        env['deepSV1'] = env['soilA1'] * env['deepsD1']
        env['deepSV2'] = env['soilA2'] * env['deepsD2']
        env['deepSV3'] = env['soilA3'] * env['deepsD3']
        env['deepSV4'] = env['soilA4'] * env['deepsD4']
        # volumn percentage calculation
        if env['fWaterV'] == 0:
            env['fSSVf'] = 0
        else:
            env['fSSVf'] = (env['fSSV'] / env['fWaterV'])

        if env['sWaterV'] == 0:
            env['sSSVf'] = 0
        else:
            env['sSSVf'] = (env['sSSV'] / env['sWaterV'])

        env['fwVf'] = 1 - env['fSSVf']
        env['swVf'] = 1 - env['sSSVf']
        env['aerVf'] = (env['aerV'] / env['areaV'])
        env['airVf'] = 1 - env['aerVf']
        env['soilSC1'] = 1 - env['soilWC1'] - env['soilAC1']
        env['soilSC2'] = 1 - env['soilWC2'] - env['soilAC2']
        env['soilSC3'] = 1 - env['soilWC3'] - env['soilAC3']
        env['soilSC4'] = 1 - env['soilWC4'] - env['soilAC4']

        # density calculation
        # soil bulk density (kg/m3)
        env['soilP1'] = env['dSS1'] * env['soilSC1'] + env['freshwP'] * env['soilWC1'] + env['airP'] * env['soilAC1']
        env['soilP2'] = env['dSS2'] * env['soilSC2'] + env['freshwP'] * env['soilWC2'] + env['airP'] * env['soilAC2']
        env['soilP3'] = env['dSS3'] * env['soilSC3'] + env['freshwP'] * env['soilWC3'] + env['airP'] * env['soilAC3']
        env['soilP4'] = env['dSS4'] * env['soilSC4'] + env['freshwP'] * env['soilWC4'] + env['airP'] * env['soilAC4']

        # CN values will be used in soilRunoff.py
        env['CN1'] = 1000.0 / env['CN1'] - 10.0
        env['CN2'] = 1000.0 / env['CN2'] - 10.0
        env['CN3'] = 1000.0 / env['CN3'] - 10.0
        env['CN4'] = 1000.0 / env['CN4'] - 10.0

        return env


    def load_climate(self):
        # load climate parameters
        climate_ws = self.region_workbook.sheet_by_name('Climate')
        start_row, end_row = self.load_date()
        climate_month = climate_ws.col_values(0, start_rowx=start_row, end_rowx=end_row)
        climate_day = climate_ws.col_values(1, start_rowx=start_row, end_rowx=end_row)
        climate_year = climate_ws.col_values(2, start_rowx=start_row, end_rowx=end_row)
        # precipitation unit: mm/day
        climate_precip = climate_ws.col_values(3, start_rowx=start_row, end_rowx=end_row)
        # windspeed unit: m/second
        climate_windspeed = climate_ws.col_values(4, start_rowx=start_row, end_rowx=end_row)
        # water flow unit: m3/second
        climate_flow = climate_ws.col_values(5, start_rowx=start_row, end_rowx=end_row)
        # temperature unit: C
        climate_temp = climate_ws.col_values(6, start_rowx=start_row, end_rowx=end_row)
        # evaporation unit: mm
        climate_evap = climate_ws.col_values(7, start_rowx=start_row, end_rowx=end_row)

        new_datetime = []
        new_month = [int(i) for i in climate_month]
        new_day = [int(i) for i in climate_day]
        new_year = [int(i) for i in climate_year]
        date = zip(new_year, new_month, new_day)

        # unit conversion
        climate = OrderedDict()
        climate['dates'] = date
        climate['precip_mm'] = climate_precip  # mm/day
        climate['precip_m'] = [(x / 1000.0) for x in climate['precip_mm']]  # m/day
        climate['windspeed_s'] = climate_windspeed  # m/second
        climate['windspeed_d'] = [(x * 86400.0) for x in climate['windspeed_s']]  # m/day
        climate['waterflow_s'] = climate_flow  # m3/s
        climate['waterflow_d'] = [(x * 86400.0) for x in climate['waterflow_s']]  # m^3/day
        climate['temp_C'] = climate_temp  # C - celcius
        climate['temp_K'] = [(x + 273.15) for x in climate['temp_C']]  # K
        climate['evap_mm'] = climate_evap

        return climate


    def load_bg_conc(self, chem_params):
        # load background concentration
        bgConc_ws = self.release_workbook.sheet_by_name('bgConc')
        bgConc_code = bgConc_ws.col_values(1, start_rowx=2, end_rowx=None)
        bgConc_value = bgConc_ws.col_values(2, start_rowx=2, end_rowx=None)

        bgConc = {}
        for i in range(len(bgConc_code)):
            # unit conversion from kg/m^3 to mol/m^3
            # kg/m3 / kg/mol = mol/m3
            bgConc[bgConc_code[i]] = bgConc_value[i] / chem_params['molar_mass']

        return bgConc


    def load_release(self, chem_params, presence):
        # load release data
        release_ws = self.release_workbook.sheet_by_name('Release')
        start_row, end_row = self.load_date()
        # to take into account of the row of release scenario
        start_row = start_row + 1
        end_row = end_row + 1
        release_scenario = release_ws.cell_value(rowx=0, colx=1)
        release_month = release_ws.col_values(0, start_rowx=start_row, end_rowx=end_row)
        release_day = release_ws.col_values(1, start_rowx=start_row, end_rowx=end_row)
        release_year = release_ws.col_values(2, start_rowx=start_row, end_rowx=end_row)
        release_air = release_ws.col_values(3, start_rowx=start_row, end_rowx=end_row)
        release_fw = release_ws.col_values(4, start_rowx=start_row, end_rowx=end_row)
        release_fSS = release_ws.col_values(5, start_rowx=start_row, end_rowx=end_row)
        release_fwSed = release_ws.col_values(6, start_rowx=start_row, end_rowx=end_row)
        release_sw = release_ws.col_values(7, start_rowx=start_row, end_rowx=end_row)
        release_sSS = release_ws.col_values(8, start_rowx=start_row, end_rowx=end_row)
        release_swSed = release_ws.col_values(9, start_rowx=start_row, end_rowx=end_row)
        release_soil1 = release_ws.col_values(10, start_rowx=start_row, end_rowx=end_row)
        release_dsoil1 = release_ws.col_values(11, start_rowx=start_row, end_rowx=end_row)
        release_soil2 = release_ws.col_values(12, start_rowx=start_row, end_rowx=end_row)
        release_dsoil2 = release_ws.col_values(13, start_rowx=start_row, end_rowx=end_row)
        release_soil3 = release_ws.col_values(14, start_rowx=start_row, end_rowx=end_row)
        release_dsoil3 = release_ws.col_values(15, start_rowx=start_row, end_rowx=end_row)
        release_soil4 = release_ws.col_values(16, start_rowx=start_row, end_rowx=end_row)
        release_dsoil4 = release_ws.col_values(17, start_rowx=start_row, end_rowx=end_row)
        # Create Datetime Objects
        new_datetime = []
        new_month = [int(i) for i in release_month]
        new_day = [int(i) for i in release_day]
        new_year = [int(i) for i in release_year]
        dt = zip(new_year, new_month, new_day)

        release = {}
        release['dates'] = dt

        # mol/day
        release['air'] = [(x / chem_params['molar_mass'])  for x in release_air]
        release['fw'] = [(x / chem_params['molar_mass']) for x in release_fw]
        release['fSS'] = [(x / chem_params['molar_mass']) for x in release_fSS]
        release['fwSed'] = [(x / chem_params['molar_mass']) for x in release_fwSed]
        release['sw'] = [(x / chem_params['molar_mass'])  for x in release_sw]
        release['sSS'] = [(x / chem_params['molar_mass']) for x in release_sSS]
        release['swSed'] = [(x / chem_params['molar_mass']) for x in release_swSed]
        release['soil1'] = [(x / chem_params['molar_mass']) for x in release_soil1]
        release['dsoil1'] = [(x / chem_params['molar_mass']) for x in release_dsoil1]
        release['soil2'] = [(x / chem_params['molar_mass']) for x in release_soil2]
        release['dsoil2'] = [(x / chem_params['molar_mass']) for x in release_dsoil2]
        release['soil3'] = [(x / chem_params['molar_mass'])  for x in release_soil3]
        release['dsoil3'] = [(x / chem_params['molar_mass']) for x in release_dsoil3]
        release['soil4'] = [(x / chem_params['molar_mass']) for x in release_soil4]
        release['dsoil4'] = [(x / chem_params['molar_mass']) for x in release_dsoil4]

        return release, release_scenario


    def run_loadData(self):
        # run the functions above to load all of the data
        presence = self.load_compart_presence()
        env = self.load_env_params()
        chem_params = self.load_chemParams(self.chem_type, env)
        climate = self.load_climate()
        bgConc = self.load_bg_conc(chem_params)
        release, release_scenario = self.load_release(chem_params, presence)

        return chem_params, presence, env, climate, bgConc, release, release_scenario






