import numpy as np
from datetime import datetime, timedelta
import math
from scipy.integrate import ode
import json

from Y_ion import Y_Value
from Z_non_ion import zValue
from eqDissolution import eqDissolution

from ode_non_ion import org_ode
from ode_ion import ion_ode
from ode_metal import metal_ode
from ode_nano import ode_nano

from ode_ion_process import ion_process
from ode_metal_process import metal_process
from ode_non_ion_process import org_process
from ode_nano_process import nano_process


def org_solver(start_date, time, presence, env, climate, chemParams, bgConc, release):
    V_bulk = [env['areaV'], env['fwV'], env['sedFWV'], env['swV'], env['sedSWV'], env['soilV1'], env['deepSV1'], env['soilV2'],
              env['deepSV2'], env['soilV3'], env['deepSV3'], env['soilV4'], env['deepSV4']]

    V_sub = [env['airV'], env['aerV'], env['fWaterV'], env['fSSV'], env['fSedWV'], env['fSedSV'],
             env['sWaterV'], env['sSSV'], env['sSedWV'], env['sSedSV'], env['soilAV1'], env['soilWV1'],
             env['soilSV1'], env['deepSV1'], env['soilAV2'], env['soilWV2'], env['soilSV2'], env['deepSV2'],
             env['soilAV3'], env['soilWV2'], env['soilSV3'], env['deepSV3'],
             env['soilAV4'], env['soilWV2'], env['soilSV3'], env['deepSV4']]

    # print ["%E" % e for e in V_bulk]

    # initialize the list to store fugacity values in each compartment in Pa
    funF = np.zeros((time, len(V_bulk)))
    # initialize the list to store the concentration value in each compartment in kg/m^3
    funC_bulk_kg = np.zeros((time, len(V_bulk)))
    funC_bulk_mol = np.zeros((time, len(V_bulk)))
    funC_sub_kg = np.zeros((time, len(V_sub)))
    funC_sub_mol = np.zeros((time, len(V_sub)))

    # initialize the list to store the mass value in kg
    funM_bulk_kg = np.zeros((time, len(V_bulk)))
    funM_sub_kg = np.zeros((time, len(V_sub)))

    zV = zValue(climate['temp_K'][0], chemParams['Kaw_n'], chemParams['Kp_n'], env['aerP'], chemParams['Koc_n'])
    zAirSub = zV.zAirSub()
    zAerSub = zV.zAerSub(zAirSub)
    zWaterSub = zV.zWaterSub(zAirSub)
    zFWSusSedSub = zV.zWaterSusSedSub(zWaterSub, chemParams['Kssfw_unitless'])
    zSWSusSedSub = zV.zWaterSusSedSub(zWaterSub, chemParams['Ksssw_unitless'])
    zFSedSSub = zV.zWaterSedSolidSub(zWaterSub, chemParams['Kbsfw_unitless'])
    zSSedSSub = zV.zWaterSedSolidSub(zWaterSub, chemParams['Kbssw_unitless'])
    zS1SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd1_unitless'])
    zS2SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd2_unitless'])
    zS3SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd3_unitless'])
    zS4SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd4_unitless'])
    zS1DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd1_d_unitless'])
    zS2DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd2_d_unitless'])
    zS3DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd3_d_unitless'])
    zS4DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd4_d_unitless'])
    zAirBulk = zV.zAirBulk(env['aerVf'], zAirSub, zAerSub)
    zFWBulk = zV.zWaterBulk(env['fSSVf'], zFWSusSedSub, zWaterSub)
    zSWBulk = zV.zWaterBulk(env['sSSVf'], zSWSusSedSub, zWaterSub)
    zFWSedimentBulk = zV.zSedimentBulk(env['fsedpercSolid'], zWaterSub, zFSedSSub)
    zSWSedimentBulk = zV.zSedimentBulk(env['ssedpercSolid'], zWaterSub, zSSedSSub)
    zSoil1Bulk = zV.zSoilBulk(env['soilAC1'], env['soilWC1'], zAirSub, zWaterSub, zS1SolidSub)
    zSoil2Bulk = zV.zSoilBulk(env['soilAC2'], env['soilWC2'], zAirSub, zWaterSub, zS2SolidSub)
    zSoil3Bulk = zV.zSoilBulk(env['soilAC3'], env['soilWC3'], zAirSub, zWaterSub, zS3SolidSub)
    zSoil4Bulk = zV.zSoilBulk(env['soilAC4'], env['soilWC4'], zAirSub, zWaterSub, zS4SolidSub)

    Z_bulk = [zAirBulk, zFWBulk, zFWSedimentBulk, zSWBulk, zSWSedimentBulk, zSoil1Bulk, zS1DeepSSub,
              zSoil2Bulk, zS2DeepSSub, zSoil3Bulk, zS3DeepSSub, zSoil4Bulk, zS4DeepSSub]

    # initial conditions for solver step 1
    bgConcNames = ['air', 'fw', 'fSedS', 'sw', 'sSedS', 'soilS1', 'dsoil1', 'soilS2', 'dsoil2',
                   'soilS3', 'dsoil3', 'soilS4', 'dsoil4']

    for i in range(len(V_bulk)):
        # mol/m3
        funC_bulk_mol[0, i] = bgConc[bgConcNames[i]]

    for i in range(len(V_bulk)):
        try:
            # mol/m3 / mol/(Pa-m^3) = Pa
            funF[0, i] = funC_bulk_mol[0,i]/Z_bulk[i]  # fugacity values from concentration and Z
        except:
            funF[0, i] = 0

    # initialize the first solution of fugacity
    f = [funF[0]]
    date_array = []
    process_array = np.zeros((time, 85))
    start_day = datetime.strptime(start_date, "%Y %m %d")

    for i in range(time):
        print (i)
        r = ode(org_ode).set_integrator('vode', method='bdf', order=5, with_jacobian=True,
                                        nsteps= 5000, rtol=1e-6, atol=1e-14)
        r.set_initial_value(f[-1], 0)
        r.set_f_params(i, presence,env,climate,chemParams,release, bgConc)
        soln = r.integrate(1)
        f.append(soln)

        funF[i] = f[-1]
        date = (start_day + timedelta(days = i)).strftime('%Y %m %d')
        date_array.append(date)

        process_org = org_process(funF[i], i, env, climate, chemParams, bgConc)
        for j in range(0, len(process_org)):
            process_array[i, j] = process_org[j]

        zV = zValue(climate['temp_K'][i], chemParams['Kaw_n'], chemParams['Kp_n'], env['aerP'], chemParams['Koc_n'])
        zAirSub = zV.zAirSub()
        zAerSub = zV.zAerSub(zAirSub)
        zWaterSub = zV.zWaterSub(zAirSub)
        zFWSusSedSub = zV.zWaterSusSedSub(zWaterSub, chemParams['Kssfw_unitless'])
        zSWSusSedSub = zV.zWaterSusSedSub(zWaterSub, chemParams['Ksssw_unitless'])
        zFSedSSub = zV.zWaterSedSolidSub(zWaterSub, chemParams['Kbsfw_unitless'])
        zSSedSSub = zV.zWaterSedSolidSub(zWaterSub, chemParams['Kbssw_unitless'])
        zS1SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd1_unitless'])
        zS2SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd2_unitless'])
        zS3SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd3_unitless'])
        zS4SolidSub = zV.zSoilSolidSub(zWaterSub, chemParams['Kd4_unitless'])
        zS1DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd1_d_unitless'])
        zS2DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd2_d_unitless'])
        zS3DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd3_d_unitless'])
        zS4DeepSSub = zV.zDeepS(zWaterSub, chemParams['Kd4_d_unitless'])
        zAirBulk = zV.zAirBulk(env['aerVf'], zAirSub, zAerSub)
        zFWBulk = zV.zWaterBulk(env['fSSVf'], zFWSusSedSub, zWaterSub)
        zSWBulk = zV.zWaterBulk(env['sSSVf'], zSWSusSedSub, zWaterSub)
        zFWSedimentBulk = zV.zSedimentBulk(env['fsedpercSolid'], zWaterSub, zFSedSSub)
        zSWSedimentBulk = zV.zSedimentBulk(env['ssedpercSolid'], zWaterSub, zSSedSSub)
        zSoil1Bulk = zV.zSoilBulk(env['soilAC1'], env['soilWC1'], zAirSub, zWaterSub, zS1SolidSub)
        zSoil2Bulk = zV.zSoilBulk(env['soilAC2'], env['soilWC2'], zAirSub, zWaterSub, zS2SolidSub)
        zSoil3Bulk = zV.zSoilBulk(env['soilAC3'], env['soilWC3'], zAirSub, zWaterSub, zS3SolidSub)
        zSoil4Bulk = zV.zSoilBulk(env['soilAC4'], env['soilWC4'], zAirSub, zWaterSub, zS4SolidSub)

        Z_bulk = [zAirBulk, zFWBulk, zFWSedimentBulk, zSWBulk, zSWSedimentBulk, zSoil1Bulk, zS1DeepSSub,
                  zSoil2Bulk, zS2DeepSSub, zSoil3Bulk, zS3DeepSSub, zSoil4Bulk, zS4DeepSSub]
        Z_sub = [zAirSub, zAerSub, zWaterSub, zFWSusSedSub, zWaterSub, zFSedSSub, zWaterSub, zSWSusSedSub,
                 zWaterSub, zSSedSSub, zAirSub, zWaterSub, zS1SolidSub, zS1DeepSSub, zAirSub, zWaterSub, zS2SolidSub, zS2DeepSSub,
                 zAirSub, zWaterSub, zS3SolidSub, zS3DeepSSub, zAirSub, zWaterSub, zS4SolidSub, zS4DeepSSub]

        # multiply fugacity*Zvalue to get the concentration in each compartment
        for j in range(len(V_bulk)):
            # unit: Pa * mol/m^3-Pa * kg/mol = kg/m^3
            funC_bulk_kg[i, j] = f[-1][j] * Z_bulk[j] * chemParams['molar_mass']
            funC_bulk_mol[i, j] = f[-1][j] * Z_bulk[j]
            funM_bulk_kg[i, j] = funC_bulk_kg[i, j] * V_bulk[j]


        # calculate for the subcompartments concentration and mass
        funC_sub_mol[i, 0] = funF[i][0] * Z_sub[0]  # air
        funC_sub_mol[i, 1] = funF[i][0] * Z_sub[1]  # aerosol
        funC_sub_mol[i, 2] = funF[i][1] * Z_sub[2]  # fw
        funC_sub_mol[i, 3] = funF[i][1] * Z_sub[3]  # fw sus sed
        funC_sub_mol[i, 4] = funF[i][2] * Z_sub[4]  # fw sediment water
        funC_sub_mol[i, 5] = funF[i][2] * Z_sub[5]  # fw sediment sediment
        funC_sub_mol[i, 6] = funF[i][3] * Z_sub[6]  # sw
        funC_sub_mol[i, 7] = funF[i][3] * Z_sub[7]  # sw sus sed
        funC_sub_mol[i, 8] = funF[i][4] * Z_sub[8]  # sw sediment water
        funC_sub_mol[i, 9] = funF[i][4] * Z_sub[9]  # sw sediment sediment
        funC_sub_mol[i, 10] = funF[i][5] * Z_sub[10]  # soil 1 air
        funC_sub_mol[i, 11] = funF[i][5] * Z_sub[11]  # soil 1 water
        funC_sub_mol[i, 12] = funF[i][5] * Z_sub[12]  # soil 1 solid
        funC_sub_mol[i, 13] = funF[i][6] * Z_sub[13]  # soil 1 deep
        funC_sub_mol[i, 14] = funF[i][7] * Z_sub[14]  # soil 2 air
        funC_sub_mol[i, 15] = funF[i][7] * Z_sub[15]  # soil 2 water
        funC_sub_mol[i, 16] = funF[i][7] * Z_sub[16]  # soil 2 solid
        funC_sub_mol[i, 17] = funF[i][8] * Z_sub[17]  # soil 2 deep
        funC_sub_mol[i, 18] = funF[i][9] * Z_sub[18]  # soil 3 air
        funC_sub_mol[i, 19] = funF[i][9] * Z_sub[19]  # soil 3 water
        funC_sub_mol[i, 20] = funF[i][9] * Z_sub[20]  # soil 3 solid
        funC_sub_mol[i, 21] = funF[i][10] * Z_sub[21]  # soil 3 deep
        funC_sub_mol[i, 22] = funF[i][11] * Z_sub[22]  # soil 4 air
        funC_sub_mol[i, 23] = funF[i][11] * Z_sub[23]  # soil 4 water
        funC_sub_mol[i, 24] = funF[i][11] * Z_sub[24]  # soil 4 solid
        funC_sub_mol[i, 25] = funF[i][12] * Z_sub[25]  # soil 4 deep

        for j in range(len(V_sub)):
            # unit: mol/m3 * kg/mol = kg/m^3
            funC_sub_kg[i, j] = funC_sub_mol[i, j] * chemParams['molar_mass']
            funM_sub_kg[i, j] = funC_sub_kg[i, j] * V_sub[j]

        output_array = [funC_bulk_kg, funC_sub_kg, funM_bulk_kg, funM_sub_kg]
        output_array = remove_floating_values(output_array)

    return date_array, process_array, output_array[0], output_array[1], output_array[2], output_array[3]


def ion_solver(chem_type, start_date, time, presence, env, climate, chemParams, bgConc, release):

    with open('./IonizableChem_helper.json') as f:
        data = json.load(f)

    subcompart_map = data['subcompart_map']

    subcompart_list = ['air', 'aer', 'fw', 'fSS', 'fSedW', 'fSedS', 'sw', 'sSS', 'sSedW', 'sSedS',
                       'soilA1', 'soilW1', 'soilS1', 'deepS1', 'soilA2', 'soilW2', 'soilS2', 'deepS2',
                       'soilA3', 'soilW3', 'soilS3', 'deepS3', 'soilA4', 'soilW4', 'soilS4', 'deepS4']

    compart_list = ['air', 'fw', 'fwSed', 'sw', 'swSed', 'soil1', 'deepS1', 'soil2', 'deepS2',
                    'soil3', 'deepS3', 'soil4', 'deepS4']

    Y_val = Y_Value(chem_type, chemParams, env)

    Z_ij_dict, Z_ij_dict_sub = Y_val.Z_ij()
    Y_ij_dict = Y_val.Y_ij()
    X_ij_dict = Y_val.X_ij()
    Z_i_dict = Y_val.Z_i()

    # initial conditions for solver step 1
    bgConcNames = ['air', 'fw', 'fSedS', 'sw', 'sSedS', 'soilS1', 'dsoil1', 'soilS2', 'dsoil2',
                   'soilS3', 'dsoil3', 'soilS4', 'dsoil4']

    # initialize the list to store aquivalence values in each compartment
    funF = np.zeros((time, len(compart_list)))
    # initialize the list to store the concentration value in each compartment in kg/m^3
    # for ionizable organic, 1 - neutral, 2 - ionic
    # for metal, 1 - particle, 2 - colloidal, 3 - dissolved
    funC_mol_1 = np.zeros((time, len(compart_list)))
    funC_mol_2 = np.zeros((time, len(compart_list)))
    funC_mol_3 = np.zeros((time, len(compart_list)))
    funC_mol_1_sub = np.zeros((time, len(subcompart_list)))
    funC_mol_2_sub = np.zeros((time, len(subcompart_list)))
    funC_mol_3_sub = np.zeros((time, len(subcompart_list)))
    funC_kg_1 = np.zeros((time, len(compart_list)))
    funC_kg_2 = np.zeros((time, len(compart_list)))
    funC_kg_3 = np.zeros((time, len(compart_list)))
    funM_kg_1 = np.zeros((time, len(compart_list)))
    funM_kg_2 = np.zeros((time, len(compart_list)))
    funM_kg_3 = np.zeros((time, len(compart_list)))

    funC_kg_1_sub = np.zeros((time, len(subcompart_list)))
    funC_kg_2_sub = np.zeros((time, len(subcompart_list)))
    funC_kg_3_sub = np.zeros((time, len(subcompart_list)))
    funM_kg_1_sub = np.zeros((time, len(subcompart_list)))
    funM_kg_2_sub = np.zeros((time, len(subcompart_list)))
    funM_kg_3_sub = np.zeros((time, len(subcompart_list)))


    for i in range(len(bgConcNames)):
        compart = compart_list[i]
        # mol/m3 / unitless = mol/m3
        funF[0, i] = bgConc[bgConcNames[i]] / Z_i_dict[compart]

    # initialize the first solution of aquivalency
    f = [funF[0]]
    date_array = []

    if chem_type == 'IonizableOrganic':
        process_array = np.zeros((time, 86))
    else:
        process_array = np.zeros((time, 54))

    start_day = datetime.strptime(start_date, "%Y %m %d")


    V_bulk = [env['areaV'], env['fwV'], env['sedFWV'], env['swV'], env['sedSWV'], env['soilV1'], env['deepSV1'],
              env['soilV2'], env['deepSV2'], env['soilV3'], env['deepSV3'], env['soilV4'], env['deepSV4']]

    V_sub = [env['airV'], env['aerV'], env['fWaterV'], env['fSSV'], env['fSedWV'], env['fSedSV'],
             env['sWaterV'], env['sSSV'], env['sSedWV'], env['sSedSV'], env['soilAV1'], env['soilWV1'],
             env['soilSV1'], env['deepSV1'], env['soilAV2'], env['soilWV2'], env['soilSV2'], env['deepSV2'],
             env['soilAV3'], env['soilWV2'], env['soilSV3'], env['deepSV3'],
             env['soilAV4'], env['soilWV2'], env['soilSV3'], env['deepSV4']]

    for i in range(time):
        print (i)
        if chem_type == 'IonizableOrganic':
            r = ode(ion_ode).set_integrator('vode', method='bdf', order=5, with_jacobian=True,
                                            nsteps=5000, rtol=1e-6, atol=1e-14)
            r.set_initial_value(f[-1], 0)
            r.set_f_params(i, presence, env, chemParams, climate, release, bgConc,
                           Z_ij_dict, Z_ij_dict_sub, Y_ij_dict, X_ij_dict, Z_i_dict)
            soln = r.integrate(1)
            f.append(soln)
            funF[i] = f[-1]

            process_ion = ion_process(funF[i], i, env, chemParams, climate, bgConc, Z_ij_dict, Z_ij_dict_sub, Y_ij_dict, X_ij_dict)
            for j in range(0, len(process_ion)):
                process_array[i,j] = process_ion[j]

        elif chem_type == 'Metal':
            r = ode(metal_ode).set_integrator('vode', method='bdf', order=5, with_jacobian=True,
                                            nsteps=1000, rtol=1e-9, atol=1e-10)
            r.set_initial_value(f[-1], 0)
            r.set_f_params(i, presence, env, chemParams, climate, release, bgConc,
                           Z_ij_dict, Y_ij_dict, Z_i_dict)
            soln = r.integrate(1)
            f.append(soln)
            funF[i] = f[-1]

            process_metal = metal_process(funF[i], i, env, chemParams, climate, bgConc, Z_ij_dict, Y_ij_dict)
            for j in range(0, len(process_metal)):
                process_array[i, j] = process_metal[j]


        date = (start_day + timedelta(days=i)).strftime('%Y %m %d')
        date_array.append(date)

        # multiply aquavalency*Zvalue to get the concentration in each compartment
        # Cij = Qij*Zij = Qit*Yij*Zij

        for j in range(len(compart_list)):
            compart = compart_list[j]
            if chem_type == 'IonizableOrganic':
                funC_mol_1[i, j] = funF[i][j] * Y_ij_dict[compart][0] * Z_ij_dict[compart][0]
                funC_mol_2[i, j] = funF[i][j] * Y_ij_dict[compart][1] * Z_ij_dict[compart][1]
                funC_kg_1[i, j] = funC_mol_1[i, j] * chemParams['molar_mass']
                funC_kg_2[i, j] = funC_mol_2[i, j] * chemParams['molar_mass']
                funM_kg_1[i, j] = funC_kg_1[i, j] * V_bulk[j]
                funM_kg_2[i, j] = funC_kg_2[i, j] * V_bulk[j]
            elif chem_type == 'Metal':
                # for particulate, need to adjust the conc. to particle volumes
                funC_mol_1[i, j] = funF[i][j] * Y_ij_dict[compart][0] * Z_ij_dict[compart][0]
                funC_mol_2[i, j] = funF[i][j] * Y_ij_dict[compart][1] * Z_ij_dict[compart][1]
                funC_mol_3[i, j] = funF[i][j] * Y_ij_dict[compart][2] * Z_ij_dict[compart][2]

                # mol/m3 * kg/mol = kg/m3
                funC_kg_1[i, j] = funC_mol_1[i, j] * chemParams['molar_mass']
                funC_kg_2[i, j] = funC_mol_2[i, j] * chemParams['molar_mass']
                funC_kg_3[i, j] = funC_mol_3[i, j] * chemParams['molar_mass']
                funM_kg_1[i, j] = funC_kg_1[i, j] * V_bulk[j]
                funM_kg_2[i, j] = funC_kg_2[i, j] * V_bulk[j]
                funM_kg_3[i, j] = funC_kg_3[i, j] * V_bulk[j]

        # calculate for the subcompartments concentration
        for j in range(len(subcompart_list)):
            subcompart = subcompart_list[j]
            if chem_type == 'IonizableOrganic':
                funC_mol_1_sub[i, j] = funF[i][subcompart_map[subcompart][0]] * Y_ij_dict[subcompart_map[subcompart][1]][0] * Z_ij_dict_sub[subcompart][0]
                funC_mol_2_sub[i, j] = funF[i][subcompart_map[subcompart][0]] * Y_ij_dict[subcompart_map[subcompart][1]][1] * Z_ij_dict_sub[subcompart][1]
                funC_kg_1_sub[i, j] = funC_mol_1_sub[i, j] * chemParams['molar_mass']
                funC_kg_2_sub[i, j] = funC_mol_2_sub[i, j] * chemParams['molar_mass']
                funM_kg_1_sub[i, j] = funC_kg_1_sub[i, j] * V_sub[j]
                funM_kg_2_sub[i, j] = funC_kg_2_sub[i, j] * V_sub[j]

        output_array = [funC_kg_1, funC_kg_2, funC_kg_3, funC_kg_1_sub, funC_kg_2_sub, funC_kg_3_sub,
                        funM_kg_1, funM_kg_2, funM_kg_3, funM_kg_1_sub, funM_kg_2_sub, funM_kg_3_sub]
        # output_array = remove_floating_values(output_array)


    return date_array, process_array, output_array[0], output_array[1], output_array[2], output_array[3], \
           output_array[4], output_array[5], output_array[6], output_array[7], output_array[8], output_array[9], \
           output_array[10], output_array[11]


def nano_solver(start_date, time, presence, env, climate, ENM, bgConc, release):
    # %   Nano solver function solves the giant differential equation over time
    # %   in a for loop where the coefficients are dependent on the previous solution from the
    # %   previous time step
    # %   Inputs include simulation time, presence of compartments, the
    # %   environment, the climate, the ENM, the background starting
    # %   concentrations, and the releases

    # %% Volume vector
    # %  Needed for calculations
    # V = [env['airV'], env['aerV'], env['freshwV'], env['freshssV'], env['sedFWV'], env['seawV'], env['seassV'],
    #      env['sedSWV'], env['soilV1'], env['soilwV1'], env['soilV2'], env['soilwV2'], env['soilV3'], env['soilwV3'],
    #      env['soilV4'], env['soilwV4'], env['freshwV'], env['sedFWV'], env['seawV'], env['sedSWV'], env['soilwV1'],
    #      env['soilwV2'], env['soilwV3'], env['soilwV4'], env['deepsV1'], env['deepsV2'], env['deepsV3'], env['deepsV4']]

    # for solids, divide the volume of bulk compartment instead of just solids?
    V = [env['airV'], env['aerV'], env['freshwV'], env['freshwV'], env['sedFWV'], env['seawV'], env['seawV'],
         env['sedSWV'], env['soilV1'], env['soilV1'], env['soilV2'], env['soilV2'], env['soilV3'], env['soilV3'],
         env['soilV4'], env['soilwV4'], env['freshwV'], env['sedFWV'], env['seawV'], env['sedSWV'], env['soilwV1'],
         env['soilwV2'], env['soilwV3'], env['soilwV4'], env['deepsV1'], env['deepsV2'], env['deepsV3'], env['deepsV4']]

    # %% Set initial conditions/background concentration for the model to solve
    bgConcNames = list(bgConc.keys())
    funC = np.zeros((time, len(V)))
    funM = np.zeros((time, len(V)))

    # nanoFate keep tracks of three ENM states: 1) free nanoparticles in water; 2) ENM particles with solids; 3) ENM dissolved
    funC_kg_1 = np.zeros((time, 13)) # aer, fwSS, fwSed, swSS, swSed, soil1, deepS1, soil2, deepS2, soil3, deepS3, soil4, deepS4
    funC_kg_2 = np.zeros((time, 13)) # air, fw, sw, soil1w, soil2w, soil3w, soil4w
    funC_kg_3 = np.zeros((time, 13)) # fwDis, fwSedDis, swDis, swSedDis, soil1wDis, soil2wDis, soil3wDis, soil4wDis

    funM_kg_1 = np.zeros((time, 13))
    funM_kg_2 = np.zeros((time, 13))
    funM_kg_3 = np.zeros((time, 13))

    funC_2_array = [0, 1, 3, 5, 7, 9, 11]
    funC_3_array = [1, 2, 3, 4, 5, 7, 9, 11]

    particle_array = [1, 3, 4, 6, 7, 8, 24, 10, 25, 12, 26, 14, 27]
    freeNano_array = [0, 2, 5, 9, 11, 13, 15]
    dissolved_array = [16, 17, 18, 19, 20, 21, 22, 23]

    date_array = []
    start_day = datetime.strptime(start_date, "%Y %m %d")
    process_array = np.zeros((time, 73))

    # % initial conditions for solver step 1 found in first row
    for i in range(len(V)):
        funC[0, i] = bgConc[bgConcNames[i]]

    for i in range(len(V)):
        funM[0, i] = funC[0, i] * V[i]  # % mass values from concentration and volume
        if math.isnan(funM[0, i]):
            funM[0, i] = 0

    # %% Dissolution equilibrium across range of concentrations and pHs
    # % calculates the maximum possible dissolution in the system given the ENM
    # % metal and the pH of the water system.  Dissolution cannot exceed this
    # % value
    DIS = eqDissolution(ENM['ENM'], env['freshwpH'], env['seawpH'], env['soilWpH1'], env['soilWpH2'], env['soilWpH3'],
                        env['soilWpH4'], presence)

    f = [funM[0]]

    # matched tolerance to matlab, can't go lower and still get a match and run matlab
    for i in range(time):
        print (i)
        # -9 and -10 are a statistical match to matlab
        r = ode(ode_nano).set_integrator('vode', method='bdf', with_jacobian=True,
                                         nsteps=5000, rtol=1e-6, atol=1e-14)
        r.set_initial_value(f[-1], 0)
        r.set_f_params(i, V, presence, env, climate, ENM, release, bgConc, DIS, time)
        soln = r.integrate(1)
        f.append(soln)

        # % output is mass values
        funM[i] = f[-1]  # % mass values for time=i
        date = (start_day + timedelta(days=i)).strftime('%Y %m %d')
        date_array.append(date)

        process_nano = nano_process(funM[i],i,V,presence,env,climate,ENM,DIS)

        for j in range(0, len(process_nano)):
            process_array[i, j] = process_nano[j]

        for j in range(len(V)):
            with np.errstate(divide='ignore', invalid='ignore'):
                c = np.true_divide(f[-1][j], V[j])  # % concentration values for time=i
                funC[i, j] = np.nan_to_num(c)

        for j in range(0, 13):
            funC_kg_1[i, j] = funC[i, particle_array[j]]
            funM_kg_1[i, j] = funM[i, particle_array[j]]
            if j in funC_2_array:
                index = funC_2_array.index(j)
                funC_kg_2[i, j] = funC[i, freeNano_array[index]]
                funM_kg_2[i, j] = funM[i, freeNano_array[index]]
            if j in funC_3_array:
                index = funC_3_array.index(j)
                funC_kg_3[i, j] = funC[i, dissolved_array[index]]
                funM_kg_3[i, j] = funM[i, dissolved_array[index]]

    bulk_M, bulk_C = bulkCalculator(time, V, funM, presence)

    # convert units back to kg/m3
    bulk_M = bulk_M / (10 ** 9)
    bulk_C = bulk_C / (10 ** 9)
    funC = funC / (10 ** 9)
    funM = funM / (10 ** 9)
    funC_kg_1 = funC_kg_1 / (10 ** 9)
    funC_kg_2 = funC_kg_2 / (10 ** 9)
    funC_kg_3 = funC_kg_3 / (10 ** 9)
    funM_kg_1 = funM_kg_1 / (10 ** 9)
    funM_kg_2 = funM_kg_2 / (10 ** 9)
    funM_kg_3 = funM_kg_3 / (10 ** 9)

    output_array = [bulk_C, funC, bulk_M, funM, funC_kg_1, funC_kg_2, funC_kg_3, funM_kg_1, funM_kg_2, funM_kg_3]
    output_array = remove_floating_values(output_array)

    return date_array, process_array, output_array[0], output_array[1], output_array[2], output_array[3], \
           output_array[4], output_array[5], output_array[6], output_array[7], output_array[8], output_array[9]


def bulkCalculator(time, V, funM, presence):
    bulk_C = np.zeros((time, 9))
    bulk_M = np.zeros((time, 9))
    for i in range(time):
        # BULK CONCENTRATIONS DO NOT INCLUDE DISSOLVED
        # bulk air
        if presence['air'] == 1:
            bulk_M[i, 0] = funM[i, 0] + funM[i, 1]
            bulk_C[i, 0] = (funM[i, 0] + funM[i, 1]) / (V[0] + V[1])
        else:
            bulk_C[i, 0] = 0
        # bulk freshwater column
        if presence['fw'] == 1:
            bulk_M[i, 1] = funM[i, 2] + funM[i, 3]
            bulk_C[i, 1] = (funM[i, 2] + funM[i, 3]) / (V[2] + V[3])
            # bulk freshwater sediment
            bulk_M[i, 2] = funM[i, 4]
            bulk_C[i, 2] = (funM[i, 4] / V[4])
        else:
            bulk_C[i, 1] = 0
            bulk_C[i, 2] = 0
        # bulk marine column
        if presence['sw'] == 1:
            bulk_M[i, 3] = funM[i, 5] + funM[i, 6]
            bulk_C[i, 3] = (funM[i, 5] + funM[i, 6]) / (V[5] + V[6])
            # bulk marine sediment
            bulk_M[i, 4] = funM[i, 7]
            bulk_C[i, 4] = (funM[i, 7] / V[7])
        else:
            bulk_C[i, 3] = 0
            bulk_C[i, 4] = 0
        # bulk undeveloped soil
        if presence['soil1'] == 1:
            bulk_M[i, 5] = funM[i, 8] + funM[i, 9]
            bulk_C[i, 5] = (funM[i, 8] + funM[i, 9]) / (V[8] + V[9])
        else:
            bulk_C[i, 5] = 0
        # bulk urban soil
        if presence['soil2'] == 1:
            bulk_M[i, 6] = funM[i, 10] + funM[i, 11]
            bulk_C[i, 6] = (funM[i, 10] + funM[i, 11]) / (V[10] + V[11])
        else:
            bulk_C[i, 6] = 0
        # bulk agricultural soil
        if presence['soil3'] == 1:
            bulk_M[i, 7] = funM[i, 12] + funM[i, 13]
            bulk_C[i, 7] = (funM[i, 12] + funM[i, 13]) / (V[12] + V[13])
        else:
            bulk_C[i, 7] = 0
        # bulk biosolids agricultural soil
        if presence['soil4'] == 1:
            bulk_M[i, 8] = funM[i, 14] + funM[i, 15]
            bulk_C[i, 8] = (funM[i, 14] + funM[i, 15]) / (V[14] + V[15])
        else:
            bulk_C[i, 8] = 0

    return bulk_M, bulk_C


def remove_floating_values(array_list):
    new_array_list = []
    for array in array_list:
        array = array.clip(min=0)
        new_array_list.append(array)
    return new_array_list


# def remove_floating_values(array_list):
#     for array in array_list:
#         row_count = len(array)
#         col_count = len(array[0])
#         for row in range(0, row_count):
#             for col in range(0, col_count):
#                 # if -1e-20 < array[row][col] < 1e-20:
#                 #     array[row][col] = 0
#                 if array[row][col] < 0:
#                     array[row][col] = 0
#
#     return array_list

