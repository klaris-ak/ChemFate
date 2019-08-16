#!/usr/bin/env python
import xlrd
import collections
from collections import OrderedDict
from datetime import datetime
import numpy as np
from advective_processes_nano import lsFactor

#################################################################
#
#   CLiCC Nanomaterials F&T Model Developed by Dr. Garner
#	Model original devloped in MATLAB
#   Date: August 16th, 2016
#   Converted to Python by Dillon Elsbury
#   Updated by Dr. Garner on July 26, 2017
#
#################################################################


def load_bgConc(filename, sheetname, presence):
    chem_workbook = xlrd.open_workbook(filename)
    bgValues_worksheet = chem_workbook.sheet_by_name(sheetname)
    bgValues_code = bgValues_worksheet.col_values(1, start_rowx=2, end_rowx=None)
    bgValues_value = bgValues_worksheet.col_values(2, start_rowx=2, end_rowx=None)
    # convert to units of ug/m3
    ## always converts units so need to make sure that inputs are always kg
    bgValues_value2 = [x * 10 ** 9 for x in bgValues_value]

    bgValues = OrderedDict(zip(bgValues_code, bgValues_value2))

    if presence['air'] == 0:
        bgValues['air'] = 0
        bgValues['aer'] = 0
        bgValues['gairc_n'] = 0
    if presence['fw'] == 0:
        bgValues['fw'] = 0
        bgValues['fSS'] = 0
        bgValues['fSedS'] = 0
        bgValues['fwdis'] = 0
        bgValues['fwSeddis'] = 0
    if presence['sw'] == 0:
        bgValues['sw'] = 0
        bgValues['sSS'] = 0
        bgValues['sSedS'] = 0
        bgValues['swdis'] = 0
        bgValues['swSeddis'] = 0
    if presence['soil1'] == 0:
        bgValues['soilS1'] = 0
        bgValues['soilW1'] = 0
        bgValues['dsoil1'] = 0
        bgValues['soilW1dis'] = 0
    if presence['soil2'] == 0:
        bgValues['soilS2'] = 0
        bgValues['soilW2'] = 0
        bgValues['dsoil2'] = 0
        bgValues['soilW2dis'] = 0
    if presence['soil3'] == 0:
        bgValues['soilS3'] = 0
        bgValues['soilW3'] = 0
        bgValues['dsoil3'] = 0
        bgValues['soilW3dis'] = 0
    if presence['soil4'] == 0:
        bgValues['soilS4'] = 0
        bgValues['soilW4'] = 0
        bgValues['dsoil4'] = 0
        bgValues['soilW4dis'] = 0

    bgConc = collections.OrderedDict()
    bgConc['A'] = bgValues['air']
    bgConc['Aer'] = bgValues['aer']
    bgConc['fW'] = bgValues['fw']
    bgConc['fSS'] = bgValues['fSS']
    bgConc['fwSed'] = bgValues['fSedS']
    bgConc['sW'] = bgValues['sw']
    bgConc['sSS'] = bgValues['sSS']
    bgConc['swSed'] = bgValues['sSedS']
    bgConc['S1'] = bgValues['soilS1']
    bgConc['soilW1'] = bgValues['soilW1']
    bgConc['S2'] = bgValues['soilS2']
    bgConc['soilW2'] = bgValues['soilW2']
    bgConc['S3'] = bgValues['soilS3']
    bgConc['soilW3'] = bgValues['soilW3']
    bgConc['S4'] = bgValues['soilS4']
    bgConc['soilW4'] = bgValues['soilW4']
    bgConc['fWdis'] = bgValues['fwdis']
    bgConc['fWSeddis'] = bgValues['fwSeddis']
    bgConc['Swdis'] = bgValues['swdis']
    bgConc['swSeddis'] = bgValues['swSeddis']
    bgConc['soilW1dis'] = bgValues['soilW1dis']
    bgConc['soilW2dis'] = bgValues['soilW2dis']
    bgConc['soilW3dis'] = bgValues['soilW3dis']
    bgConc['soilW4dis'] = bgValues['soilW4dis']
    bgConc['dsoil1'] = bgValues['dsoil1']
    bgConc['dsoil2'] = bgValues['dsoil2']
    bgConc['dsoil3'] = bgValues['dsoil3']
    bgConc['dsoil4'] = bgValues['dsoil4']
    bgConc['gairc'] = bgValues['gairc_n']

    return bgConc


def load_climate(filename, sheetname, start_row, end_row):
    env_workbook = xlrd.open_workbook(filename)
    climate_worksheet = env_workbook.sheet_by_name(sheetname)

    # Climate Parameter Loading
    climate_month = climate_worksheet.col_values(0, start_rowx=start_row, end_rowx=end_row)
    climate_day = climate_worksheet.col_values(1, start_rowx=start_row, end_rowx=end_row)
    climate_year = climate_worksheet.col_values(2, start_rowx=start_row, end_rowx=end_row)
    climate_precip = climate_worksheet.col_values(3, start_rowx=start_row, end_rowx=end_row)
    climate_windspeed = climate_worksheet.col_values(4, start_rowx=start_row, end_rowx=end_row)
    climate_flow = climate_worksheet.col_values(5, start_rowx=start_row, end_rowx=end_row)
    climate_temp = climate_worksheet.col_values(6, start_rowx=start_row, end_rowx=end_row)
    climate_evap = climate_worksheet.col_values(7, start_rowx=start_row, end_rowx=end_row)

    names = climate_worksheet.row_values(0, start_colx=3, end_colx=None)
    headers = []

    headers.append('dates')
    for val in names:
        headers.append(val)

    # Create Datetime Objects
    new_datetime = []
    new_month = [int(i) for i in climate_month]
    new_day = [int(i) for i in climate_day]
    new_year = [int(i) for i in climate_year]
    dt = zip(new_year, new_month, new_day)
    for val in dt:
        mystring = ' '.join(map(str, val))
        dt = datetime.strptime(mystring, "%Y %m %d")
        new_datetime.append(dt)
    climate_loading = zip(new_datetime, climate_precip, climate_windspeed, climate_flow, climate_temp)

    # climate = {}
    climate = OrderedDict()
    climate['dates'] = new_datetime
    climate['precip'] = climate_precip
    climate['windspeed'] = climate_windspeed
    climate['flow'] = climate_flow
    climate['temp'] = climate_temp
    climate['evap'] = climate_evap

    return climate


def load_ENM(filename, sheetname, presence):
    chem_workbook = xlrd.open_workbook(filename)
    ENM_worksheet = chem_workbook.sheet_by_name(sheetname)

    ENM_code = ENM_worksheet.col_values(1, start_rowx=2, end_rowx=None)
    ENM_value = ENM_worksheet.col_values(2, start_rowx=2, end_rowx=None)
    ENM_loading = zip(ENM_code, ENM_value)

    # ENM = {}
    ENM = OrderedDict()
    for name, value in ENM_loading:
        ENM[name] = value

    if presence['air'] == 0:
        ENM['khetA'] = 0
    if presence['fw'] == 0:
        ENM['kdisFW'] = 0
        ENM['kdisFWsed'] == 0
        ENM['ksedFW'] = 0
        ENM['khetFW'] = 0
    if presence['sw'] == 0:
        ENM['kdisSW'] = 0
        ENM['kdisSWsed'] = 0
        ENM['ksedSW'] = 0
        ENM['khetSW'] = 0
        ENM['enrichFactor'] = 0
    if presence['soil1'] == 0:
        ENM['kdisS1'] = 0
        ENM['elutionS1'] = 0
    if presence['soil2'] == 0:
        ENM['kdisS2'] = 0
        ENM['elutionS2'] = 0
    if presence['soil3'] == 0:
        ENM['kdisS3'] = 0
        ENM['elutionS3'] = 0
    if presence['soil4'] == 0:
        ENM['kdisS4'] = 0
        ENM['elutionS4'] = 0

    ENM['density'] = ENM['density'] * 10 ** 9
    ENM['khetA'] = np.true_divide(ENM['khetA'], 10 ** 9)
    ENM['khetFW'] = np.true_divide(ENM['khetFW'], 10 ** 9)
    ENM['khetSW'] = np.true_divide(ENM['khetSW'], 10 ** 9)

    return ENM


def load_env(filename, sheetname, presence):
    env_workbook = xlrd.open_workbook(filename)
    env_worksheet = env_workbook.sheet_by_name(sheetname)

    env_code = env_worksheet.col_values(1, start_rowx=1, end_rowx=None)
    env_value = env_worksheet.col_values(2, start_rowx=1, end_rowx=None)
    env_loading = zip(env_code, env_value)

    env = OrderedDict()
    for name, value in env_loading:
        env[name] = value

    # %% convert all units from kg to pg
    # % concentration and density units are all in kg
    env['airP'] = env['airP'] * 10 ** 9
    env['dynViscAir'] = env['dynViscAir'] * 10 ** 9
    env['aerP'] = env['aerP'] * 10 ** 9
    env['aerC'] = env['aerC'] * 10 ** 9
    env['freshwP'] = env['freshwP'] * 10 ** 9
    env['dynViscFW'] = env['dynViscFW'] * 10 ** 9
    env['freshssP'] = env['freshssP'] * 10 ** 9
    env['freshssC'] = env['freshssC'] * 10 ** 9
    env['dFWSedS'] = env['dFWSedS'] * 10 ** 9
    env['seawP'] = env['seawP'] * 10 ** 9
    env['dynViscSW'] = env['dynViscSW'] * 10 ** 9
    env['seassP'] = env['seassP'] * 10 ** 9
    env['seassC'] = env['seassC'] * 10 ** 9
    env['dSWSedS'] = env['dSWSedS'] * 10 ** 9
    env['dSS1'] = env['dSS1'] * 10 ** 9
    env['dSS2'] = env['dSS2'] * 10 ** 9
    env['dSS3'] = env['dSS3'] * 10 ** 9
    env['dSS4'] = env['dSS4'] * 10 ** 9

    # Air
    # Area of Air (m^2)
    env['area'] = env['freshwA'] + env['seawA'] + env['soilA1'] + env['soilA2'] + env['soilA3'] + env['soilA4']
    # Volume of air (m^3)
    env['airV'] = env['area'] * env['airH']
    # Volume fraction of aerosols
    env['aerVf'] = np.true_divide(env['aerC'], env['aerP'])
    # Aerosols volume (m^3)
    env['aerV'] = env['aerC'] * np.true_divide(env['airV'], env['aerP'])
    # Seawater
    # Volume of seawater (m^3)
    env['seawV'] = env['seawA'] * env['seawD']
    # Seawater suspended sediment volume (m^3)
    env['seassV'] = env['seassC'] * np.true_divide(env['seawV'], env['seassP'])
    # Area of marine sediment (m^2)
    env['sedSWA'] = env['seawA']
    # Volume of marine sediment (m^3)
    env['sedSWV'] = env['sedSWA'] * env['sedSWD']
    # Density of Seawater Sediment (kg/m3)
    env['sedSWP'] = env['dSWSedS'] * env['ssedpercSolid'] + env['seawP'] * (1 - env['ssedpercSolid'])
    # Freshwater
    # Volume of freshwater (m^3)
    env['freshwV'] = env['freshwA'] * env['freshwD']
    # Freshwater suspended sediment volume (m^3)
    env['freshssV'] = env['freshssC'] * np.true_divide(env['freshwV'], env['freshssP'])
    # Area of freshwater sediment (m^2)
    env['sedFWA'] = env['freshwA']
    # Volume of freshwater sediment (m^3)
    env['sedFWV'] = env['sedFWA'] * env['sedFWD']
    # Density of Freshwater Sediment (kg/m3)
    env['sedFWP'] = env['dFWSedS'] * env['fsedpercSolid'] + env['freshwP'] * (1 - env['fsedpercSolid'])
    # Soil 1
    # Soil Volume 1 (m^3)
    env['soilV1'] = env['soilA1'] * env['soilD1'] * (1 - env['soilWC1'] - env['soilAC1'])
    # Soil Density 1 (kg/m3)
    env['soilP1'] = env['dSS1'] * (1 - env['soilWC1'] - env['soilAC1']) + env['freshwP'] * env['soilWC1'] + env[
        'airP'] * env['soilAC1']
    # Soil Water Volume 1 (m^3)
    env['soilwV1'] = env['soilA1'] * env['soilD1'] * env['soilWC1']
    # Soil Air Volume 1 (m^3)
    env['soilaV1'] = env['soilA1'] * env['soilD1'] * env['soilAC1']
    # Length slope factor
    env['lenslope1'] = lsFactor(env['slope1'])
    # CN calculations of S
    env['CN1'] = np.true_divide(1000, env['CN1']) - 10
    # Deep Soil Volume 1 (m^3)
    env['deepsV1'] = env['soilA1'] * env['deepsD1']
    # Soil 2
    # Soil Volume 2 (m^3)
    env['soilV2'] = env['soilA2'] * env['soilD2'] * (1 - env['soilWC2'] - env['soilAC2'])
    # Soil Density 2 (kg/m3)
    env['soilP2'] = env['dSS2'] * (1 - env['soilWC2'] - env['soilAC2']) + env['freshwP'] * env['soilWC2'] + env[
        'airP'] * env['soilAC2']
    # Soil Water Volume 2 (m^3)
    env['soilwV2'] = env['soilA2'] * env['soilD2'] * env['soilWC2']
    # Soil Air Volume 2 (m^3)
    env['soilaV2'] = env['soilA2'] * env['soilD2'] * env['soilAC2']
    # Length slope factor
    env['lenslope2'] = lsFactor(env['slope2'])
    # CN calculations of S
    env['CN2'] = np.true_divide(1000, env['CN2']) - 10
    # Deep Soil Volume 2 (m^3)
    env['deepsV2'] = env['soilA2'] * env['deepsD2']
    # Soil 3
    # Soil Volume 3 (m^3)
    env['soilV3'] = env['soilA3'] * env['soilD3'] * (1 - env['soilWC3'] - env['soilAC3'])
    # Soil Density 3 (kg/m3)
    env['soilP3'] = env['dSS3'] * (1 - env['soilWC3'] - env['soilAC3']) + env['freshwP'] * env['soilWC3'] + env[
        'airP'] * env['soilAC3']
    # Soil Water Volume 3 (m^3)
    env['soilwV3'] = env['soilA3'] * env['soilD3'] * env['soilWC3']
    # Soil Air Volume 3 (m^3)
    env['soilaV3'] = env['soilA3'] * env['soilD3'] * env['soilAC3']
    # Length slope factor
    env['lenslope3'] = lsFactor(env['slope3'])
    # CN calculations of S
    env['CN3'] = np.true_divide(1000, env['CN3']) - 10
    # Deep Soil Volume 3 (m^3)
    env['deepsV3'] = env['soilA3'] * env['deepsD3']
    # Soil 4
    # Soil Volume 4 (m^3)
    env['soilV4'] = env['soilA4'] * env['soilD4'] * (1 - env['soilWC4'] - env['soilAC4'])
    # Soil Density 1 (kg/m3)
    env['soilP4'] = env['dSS4'] * (1 - env['soilWC4'] - env['soilAC4']) + env['freshwP'] * env['soilWC4'] + env[
        'airP'] * env['soilAC4']
    # Soil Water Volume 4 (m^3)
    env['soilwV4'] = env['soilA4'] * env['soilD4'] * env['soilWC4']
    # Soil Air Volume 4 (m^3)
    env['soilaV4'] = env['soilA4'] * env['soilD4'] * env['soilAC4']
    # Length slope factor
    env['lenslope4'] = lsFactor(env['slope4'])
    # CN calculations of S
    env['CN4'] = np.true_divide(1000, env['CN4']) - 10
    # Deep Soil Volume 4 (m^3)
    env['deepsV4'] = env['soilA4'] * env['deepsD4']

    if presence['air'] == 0:
        env['airA'] = 0
        env['airV'] = 0
        env['airH'] = 0
        env['dynViscAir'] = 0
        env['scavengingENM'] = 0
    if presence['aer'] == 0:
        env['aerVf'] = 0
        env['aerV'] = 0
        env['aerP'] = 0
        env['aerC'] = 0
        env['radiusParclesAer'] = 0
        env['scavenging'] = 0
    if presence['fw'] == 0:
        env['freshwA'] = 0
        env['freshwD'] = 0
        env['freshwV'] = 0
        env['freshwP'] = 0
        env['freshwpH'] = 0
        env['dynViscFW'] = 0
        env['freshssP'] = 0
        env['freshssC'] = 0
        env['freshssV'] = 0
        env['radiusParticlesFW'] = 0
        env['freshssOC'] = 0
        env['sedFWD'] = 0
        env['sedFWA'] = 0
        env['sedFWV'] = 0
        env['dFWSedS'] = 0
        env['fsedpercSolid'] = 0
        env['burialRateFW'] = 0
        env['resuspensionRateFW'] = 0
        env['fwadvfrac'] = 0
    if presence['fw'] == 0:
        env['freshwA'] = 0
        env['freshwD'] = 0
        env['freshwV'] = 0
        env['freshwP'] = 0
        env['freshwpH'] = 0
        env['dynViscFW'] = 0
        env['fwadvfrac'] = 0
    if presence['fSS'] == 0:
        env['freshssP'] = 0
        env['freshssC'] = 0
        env['freshssV'] = 0
        env['radiusParticlesFW'] = 0
        env['freshssOC'] = 0
    if presence['fSed'] == 0:
        env['sedFWD'] = 0
        env['sedFWA'] = 0
        env['sedFWV'] = 0
        env['dFWSedS'] = 0
        env['fsedpercSolid'] = 0
        env['burialRateFW'] = 0
        env['resuspensionRateFW'] = 0
    if presence['sw'] == 0:
        env['seawA'] = 0
        env['seawD'] = 0
        env['seawV'] = 0
        env['seawP'] = 0
        env['seawpH'] = 0
        env['dynViscSW'] = 0
        env['coastalA'] = 0
        env['seassP'] = 0
        env['seassC'] = 0
        env['seassV'] = 0
        env['radiusParticlesSW'] = 0
        env['marinessOC'] = 0
        env['sedSWD'] = 0
        env['sedSWA'] = 0
        env['sedSWV'] = 0
        env['dSWSedS'] = 0
        env['ssedpercSolid'] = 0
        env['sedSWOC'] = 0
        env['burialRateSW'] = 0
        env['resuspensionRateSW'] = 0
        env['swadvfrac'] = 0
    if presence['sw'] == 0:
        env['seawA'] = 0
        env['seawD'] = 0
        env['seawV'] = 0
        env['seawP'] = 0
        env['seawpH'] = 0
        env['dynViscSW'] = 0
        env['coastalA'] = 0
        env['swadvfrac'] = 0
    if presence['sSS'] == 0:
        env['seassP'] = 0
        env['seassC'] = 0
        env['seassV'] = 0
        env['radiusParticlesSW'] = 0
        env['marinessOC'] = 0
    if presence['sSed'] == 0:
        env['sedSWD'] = 0
        env['sedSWA'] = 0
        env['sedSWV'] = 0
        env['dSWSedS'] = 0
        env['ssedpercSolid'] = 0
        env['sedSWOC'] = 0
        env['burialRateSW'] = 0
        env['resuspensionRateSW'] = 0
    if presence['soil1'] == 0:
        env['soilA1'] = 0
        env['soilD1'] = 0
        env['soilV1'] = 0
        env['soilP1'] = 0
        env['dSS1'] = 0
        env['soilwV1'] = 0
        env['soilaV1'] = 0
        env['soilWpH1'] = 0
        env['deepsV1'] = 0
        env['deepswV1'] = 0
        env['deepsaV1'] = 0
        env['lenslope1'] = 0
        env['CN1'] = 0
        env['deepsD1'] = 0
        env['A1'] = 0
        env['TSV1'] = 0
        env['z_wind1'] = 0
        env['roughness1'] = 0
        env['Kconstant1'] = 0
        env['percWind1'] = 0
        env['windConstant1'] = 0
        env['percUncovered1'] = 0
        env['percSuspended1'] = 0
        env['Kfact1'] = 0
        env['slope1'] = 0
        env['cropManageFactor1'] = 0
        env['supportFactor1'] = 0
        env['leachingR1'] = 0
    if presence['soil2'] == 0:
        env['soilA2'] = 0
        env['soilD2'] = 0
        env['soilV2'] = 0
        env['soilP2'] = 0
        env['dSS2'] = 0
        env['soilwV2'] = 0
        env['soilaV2'] = 0
        env['soilWpH2'] = 0
        env['deepsV2'] = 0
        env['deepswV2'] = 0
        env['deepsaV2'] = 0
        env['lenslope2'] = 0
        env['CN2'] = 0
        env['deepsD2'] = 0
        env['A2'] = 0
        env['TSV2'] = 0
        env['z_wind2'] = 0
        env['roughness2'] = 0
        env['Kconstant2'] = 0
        env['percWind2'] = 0
        env['windConstant2'] = 0
        env['percUncovered2'] = 0
        env['percSuspended2'] = 0
        env['Kfact2'] = 0
        env['slope2'] = 0
        env['cropManageFactor2'] = 0
        env['supportFactor2'] = 0
        env['leachingR2'] = 0
    if presence['soil3'] == 0:
        env['soilA3'] = 0
        env['soilD3'] = 0
        env['soilV3'] = 0
        env['soilP3'] = 0
        env['dSS3'] = 0
        env['soilwV3'] = 0
        env['soilaV3'] = 0
        env['soilWpH3'] = 0
        env['deepsV3'] = 0
        env['deepswV3'] = 0
        env['deepsaV3'] = 0
        env['lenslope3'] = 0
        env['CN3'] = 0
        env['deepsD3'] = 0
        env['A3'] = 0
        env['TSV3'] = 0
        env['z_wind3'] = 0
        env['roughness3'] = 0
        env['Kconstant3'] = 0
        env['percWind3'] = 0
        env['windConstant3'] = 0
        env['percUncovered3'] = 0
        env['percSuspended3'] = 0
        env['Kfact3'] = 0
        env['slope3'] = 0
        env['cropManageFactor3'] = 0
        env['supportFactor3'] = 0
        env['leachingR3'] = 0
    if presence['soil4'] == 0:
        env['soilA4'] = 0
        env['soilD4'] = 0
        env['soilV4'] = 0
        env['soilP4'] = 0
        env['dSS4'] = 0
        env['soilwV4'] = 0
        env['soilaV4'] = 0
        env['soilWpH4'] = 0
        env['deepsV4'] = 0
        env['deepswV4'] = 0
        env['deepsaV4'] = 0
        env['lenslope4'] = 0
        env['CN4'] = 0
        env['deepsD4'] = 0
        env['A4'] = 0
        env['TSV4'] = 0
        env['z_wind4'] = 0
        env['roughness4'] = 0
        env['Kconstant4'] = 0
        env['percWind4'] = 0
        env['windConstant4'] = 0
        env['percUncovered4'] = 0
        env['percSuspended4'] = 0
        env['Kfact4'] = 0
        env['slope4'] = 0
        env['cropManageFactor4'] = 0
        env['supportFactor4'] = 0
        env['leachingR4'] = 0

    return env


def load_presence(filename, sheetname):
    env_workbook = xlrd.open_workbook(filename)
    presence_worksheet = env_workbook.sheet_by_name(sheetname)
    presence_code = presence_worksheet.col_values(1, start_rowx=1, end_rowx=None)
    presence_value = presence_worksheet.col_values(2, start_rowx=1, end_rowx=None)
    presence_loading = zip(presence_code, presence_value)

    # presence = {}
    presence = OrderedDict()
    for name, value in presence_loading:
        presence[name] = value

    str_list = list(filter(None, presence_value))

    # if the user ever wants to make unique changes to the presence/absence data
    # beyond simply the bulk compartments.  They would do it here for individual runs.
    # For example, the user might want freshwater without suspended sediment. By changing
    # line 51 set to 0 always, you eliminate the suspended sediment compartment.
    if presence['air'] == 1:
        presence['aer'] = 1
    else:
        presence['aer'] = 0
    if presence['fw'] == 1:
        presence['fw'] = 1
        presence['fSS'] = 1
        presence['fSed'] = 1
        presence['fSedS'] = 1
        presence['fSedW'] = 1
    else:
        presence['fw'] = 0
        presence['fSS'] = 0
        presence['fSed'] = 0
        presence['fSedS'] = 0
        presence['fSedW'] = 0
    if presence['sw'] == 1:
        presence['sw'] = 1
        presence['sSS'] = 1
        presence['sSed'] = 1
        presence['sSedS'] = 1
        presence['sSedW'] = 1
    else:
        presence['sw'] = 0
        presence['sSS'] = 0
        presence['sSed'] = 0
        presence['sSedS'] = 0
        presence['sSedW'] = 0
    if presence['soil1'] == 1:
        presence['soilS1'] = 1
        presence['soilW1'] = 1
        presence['soilDeep1'] = 1
    else:
        presence['soilS1'] = 0
        presence['soilW1'] = 0
        presence['soilDeep1'] = 0
    if presence['soil2'] == 1:
        presence['soilS2'] = 1
        presence['soilW2'] = 1
        presence['soilDeep2'] = 1
    else:
        presence['soilS2'] = 0
        presence['soilW2'] = 0
        presence['soilDeep2'] = 0
    if presence['soil3'] == 1:
        presence['soilS3'] = 1
        presence['soilW3'] = 1
        presence['soilDeep3'] = 1
    else:
        presence['soilS3'] = 0
        presence['soilW3'] = 0
        presence['soilDeep3'] = 0
    if presence['soil4'] == 1:
        presence['soilS4'] = 1
        presence['soilW4'] = 1
        presence['soilDeep4'] = 1
    else:
        presence['soilS4'] = 0
        presence['soilW4'] = 0
        presence['soilDeep4'] = 0

    return presence


def load_release(filename, sheetname, start_row, end_row, presence):
    chem_workbook = xlrd.open_workbook(filename)
    release_ws = chem_workbook.sheet_by_name(sheetname)

    release_scenario = release_ws.cell_value(rowx=0, colx=1)

    start_row = start_row + 1
    end_row = end_row + 1

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
    for val in dt:
        mystring = ' '.join(map(str, val))
        dt = datetime.strptime(mystring, "%Y %m %d")
        new_datetime.append(dt)

    # release = {}
    release = OrderedDict()
    release['dates'] = new_datetime
    release['air'] = [x * 10 ** 9 for x in release_air]
    release['fw'] = [x * 10 ** 9 for x in release_fw]
    release['fSS'] = [x * 10 ** 9 for x in release_fSS]
    release['fwSed'] = [x * 10 ** 9 for x in release_fwSed]
    release['sw'] = [x * 10 ** 9 for x in release_sw]
    release['sSS'] = [x * 10 ** 9 for x in release_sSS]
    release['swSed'] = [x * 10 ** 9 for x in release_swSed]
    release['soil1'] = [x * 10 ** 9 for x in release_soil1]
    release['dsoil1'] = [x * 10 ** 9 for x in release_dsoil1]
    release['soil2'] = [x * 10 ** 9 for x in release_soil2]
    release['dsoil2'] = [x * 10 ** 9 for x in release_dsoil2]
    release['soil3'] = [x * 10 ** 9 for x in release_soil3]
    release['dsoil3'] = [x * 10 ** 9 for x in release_dsoil3]
    release['soil4'] = [x * 10 ** 9 for x in release_soil4]
    release['dsoil4'] = [x * 10 ** 9 for x in release_dsoil4]

    if presence['air'] == 0:
        release['Air'] = [x * 0 for x in release_air]
    if presence['fw'] == 0:
        release['fw'] = [x * 0 for x in release_fw]
    if presence['sw'] == 0:
        release['sw'] = [x * 0 for x in release_sw]
    if presence['soil1'] == 0:
        release['Soil1'] = [x * 0 for x in release_soil1]
    if presence['soil2'] == 0:
        release['Soil2'] = [x * 0 for x in release_soil2]
    if presence['soil3'] == 0:
        release['Soil3'] = [x * 0 for x in release_soil3]
    if presence['soil4'] == 0:
        release['Soil4'] = [x * 0 for x in release_soil4]

    return release, release_scenario


def load_data(env_filename, enmConc_filename, enm_filename, start_date, end_date):
    # original start date will change if user does custom datasets...
    original_start_date = datetime.strptime('2005 1 1', "%Y %m %d")
    start_day = datetime.strptime(start_date, "%Y %m %d")
    start_row = (start_day - original_start_date).days + 1
    end_day = datetime.strptime(end_date, "%Y %m %d")
    time = (end_day - start_day).days + 1
    end_row = start_row + time

    presence = load_presence(env_filename, 'Presence')
    env = load_env(env_filename, 'Environment', presence)
    climate = load_climate(env_filename, 'Climate', start_row, end_row)
    bgConc = load_bgConc(enmConc_filename, 'bgConc', presence)
    ENM = load_ENM(enm_filename, 'Sheet1', presence)
    release, release_scenario = load_release(enmConc_filename, 'Release', start_row, end_row, presence)

    return time, presence, env, climate, bgConc, ENM, release, release_scenario