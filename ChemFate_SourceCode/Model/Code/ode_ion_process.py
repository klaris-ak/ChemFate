from degradation_process import Degradation
from advective_processes import AdvectiveProcess
from diffusion_process_ion import Diffusion

import numpy as np

def ion_process(Q, i, env, chemParams, climate, bgConc, Z_ij, Z_ij_sub, Y_ij, X_ij):
    # Q is the total aquivalence by each subcompartment, also Q_t, unit: mol/m3
    # Q_n = Q_t * Y_n, and Q_i = Q_t * Y_i,
    # D unit: m3/day, N unit: mol/day, Q unit: mol/m3
    # N = D*Q

    deg = Degradation()
    adv = AdvectiveProcess()
    diff = Diffusion(chemParams['molar_mass'], chemParams['Kaw_n'], climate['windspeed_s'][i])

    # Air - air, Q[0]
    Q_air_n = Q[0] * Y_ij['air'][0]
    Q_air_i = Q[0] * Y_ij['air'][1]

    # within air: 1) degradation (n), 2) advection in (n), 3) advection out (n)
    air_deg_n = deg.D_deg(env['airV'], chemParams['kDeg_air_n'], Z_ij_sub['air'][0]) * Q_air_n
    aer_deg_n = deg.D_deg(env['aerV'], chemParams['kDeg_aer_n'], Z_ij_sub['aer'][0]) * Q_air_n
    aer_deg_i = deg.D_deg(env['aerV'], chemParams['kDeg_aer_i'], Z_ij_sub['aer'][1]) * Q_air_i

    airOut_adv_n = adv.D_advec_air(climate['windspeed_d'][i], env['airA'], env['airH'],
                                   Z_ij_sub['air'][0]) * Q_air_n
    aerOut_adv_n = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'],
                                   env['airH'], env['aerVf'],
                                   Z_ij_sub['aer'][0]) * Q_air_n
    aerOut_adv_i = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'],
                                   env['airH'], env['aerVf'],
                                   Z_ij_sub['aer'][1]) * Q_air_i

    airIn_adv_n = bgConc['gairc_n'] * adv.G_advec_air(climate['windspeed_d'][i],
                                                      env['airA'], env['airH'])
    airIn_adv_i = bgConc['gairc_i'] * adv.G_advec_air(climate['windspeed_d'][i],
                                                      env['airA'], env['airH'])

    # air to other compartments: 1) dry deposition (n), 2) wet deposition (n, i),
    # 3) diffusion (n) between freshwater, seawater, soil bulk
    ## aerosol dry deposition
    k_dep_aer_dry = adv.k_dep_dry(env['aerP'], env['airP'], env['dynViscAir'], env['radiusParticlesAer'])
    aer_dep_dry_n = adv.D_dep_dry(k_dep_aer_dry, env['airA'] * env['aerVf'], Z_ij_sub['aer'][0]) * Q_air_n
    aer_dep_dry_i = adv.D_dep_dry(k_dep_aer_dry, env['airA'] * env['aerVf'], Z_ij_sub['aer'][1]) * Q_air_i
    aer_dep_dry_n_to_fSS = aer_dep_dry_n * (env['fwA'] / env['airA'])
    aer_dep_dry_n_to_sSS = aer_dep_dry_n * (env['swA'] / env['airA'])
    aer_dep_dry_n_to_soil1 = aer_dep_dry_n * (env['soilA1'] / env['airA'])
    aer_dep_dry_n_to_soil2 = aer_dep_dry_n * (env['soilA2'] / env['airA'])
    aer_dep_dry_n_to_soil3 = aer_dep_dry_n * (env['soilA3'] / env['airA'])
    aer_dep_dry_n_to_soil4 = aer_dep_dry_n * (env['soilA4'] / env['airA'])
    aer_dep_dry_i_to_fSS = aer_dep_dry_i * (env['fwA'] / env['airA'])
    aer_dep_dry_i_to_sSS = aer_dep_dry_i * (env['swA'] / env['airA'])
    aer_dep_dry_i_to_soil1 = aer_dep_dry_i * (env['soilA1'] / env['airA'])
    aer_dep_dry_i_to_soil2 = aer_dep_dry_i * (env['soilA2'] / env['airA'])
    aer_dep_dry_i_to_soil3 = aer_dep_dry_i * (env['soilA3'] / env['airA'])
    aer_dep_dry_i_to_soil4 = aer_dep_dry_i * (env['soilA4'] / env['airA'])


    ## aerosol wet deposition

    aer_dep_wet_n = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'],
                                  env['airA'] * env['aerVf'], Z_ij_sub['aer'][0]) * Q_air_n
    aer_dep_wet_i = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'],
                                  env['airA'] * env['aerVf'], Z_ij_sub['aer'][1]) * Q_air_i
    aer_dep_wet_n_to_fSS = aer_dep_wet_n * (env['fwA'] / env['airA'])
    aer_dep_wet_n_to_sSS = aer_dep_wet_n * (env['swA'] / env['airA'])
    aer_dep_wet_n_to_soil1 = aer_dep_wet_n * (env['soilA1'] / env['airA'])
    aer_dep_wet_n_to_soil2 = aer_dep_wet_n * (env['soilA2'] / env['airA'])
    aer_dep_wet_n_to_soil3 = aer_dep_wet_n * (env['soilA3'] / env['airA'])
    aer_dep_wet_n_to_soil4 = aer_dep_wet_n * (env['soilA4'] / env['airA'])

    aer_dep_wet_i_to_fSS = aer_dep_wet_i * (env['fwA'] / env['airA'])
    aer_dep_wet_i_to_sSS = aer_dep_wet_i * (env['swA'] / env['airA'])
    aer_dep_wet_i_to_soil1 = aer_dep_wet_i * (env['soilA1'] / env['airA'])
    aer_dep_wet_i_to_soil2 = aer_dep_wet_i * (env['soilA2'] / env['airA'])
    aer_dep_wet_i_to_soil3 = aer_dep_wet_i * (env['soilA3'] / env['airA'])
    aer_dep_wet_i_to_soil4 = aer_dep_wet_i * (env['soilA4'] / env['airA'])

    # air rain dissolution to water
    air_rain_diss_n = adv.D_rain_diss(climate['precip_m'][i], env['area'], Z=1.0) * Q_air_n
    air_rain_diss_n_to_fw = air_rain_diss_n * (env['freshwA'] / env['area'])
    air_rain_diss_n_to_sw = air_rain_diss_n * (env['seawA'] / env['area'])
    air_rain_diss_n_to_soil1 = air_rain_diss_n * (env['soilA1'] / env['area'])
    air_rain_diss_n_to_soil2 = air_rain_diss_n * (env['soilA2'] / env['area'])
    air_rain_diss_n_to_soil3 = air_rain_diss_n * (env['soilA3'] / env['area'])
    air_rain_diss_n_to_soil4 = air_rain_diss_n * (env['soilA4'] / env['area'])

    # diffusion
    perm_air_to_fw = diff.P_water_air(X_ij['fw'][0])
    perm_air_to_sw = diff.P_water_air(X_ij['sw'][0])
    perm_air_to_soil1 = diff.P_soil_air(X_ij['soil1'][0])
    perm_air_to_soil2 = diff.P_soil_air(X_ij['soil2'][0])
    perm_air_to_soil3 = diff.P_soil_air(X_ij['soil3'][0])
    perm_air_to_soil4 = diff.P_soil_air(X_ij['soil4'][0])

    air_diff_n_to_fw = diff.D_diffu_comp1_comp2(perm_air_to_fw, env['fwA']) * Q_air_n
    air_diff_n_to_sw = diff.D_diffu_comp1_comp2(perm_air_to_sw, env['swA']) * Q_air_n

    # diffusion (n) between soil and air, from soil side, is the combination of soil air and soil water
    air_diff_n_to_soil1 = diff.D_diffu_comp1_comp2(perm_air_to_soil1, env['soilA1']) * Q_air_n
    air_diff_n_to_soil2 = diff.D_diffu_comp1_comp2(perm_air_to_soil2, env['soilA2']) * Q_air_n
    air_diff_n_to_soil3 = diff.D_diffu_comp1_comp2(perm_air_to_soil3, env['soilA3']) * Q_air_n
    air_diff_n_to_soil4 = diff.D_diffu_comp1_comp2(perm_air_to_soil4, env['soilA4']) * Q_air_n


    # freshwater water - fw, Q[1]
    Q_fw_n = Q[1] * Y_ij['fw'][0]
    Q_fw_i = Q[1] * Y_ij['fw'][1]

    # within fw: 1)degradation (n, i), 2) advection in (n, i), 3) advection out (n, i)
    fw_deg_n = deg.D_deg(env['fwV'], chemParams['kDeg_fw_n'], Z_ij_sub['fw'][0]) * Q_fw_n
    fw_deg_i = deg.D_deg(env['fwV'], chemParams['kDeg_fw_i'], Z_ij_sub['fw'][1]) * Q_fw_i
    fSS_deg_n = deg.D_deg(env['fSSV'], chemParams['kDeg_fSS_n'], Z_ij_sub['fSS'][0]) * Q_fw_n
    fSS_deg_i = deg.D_deg(env['fSSV'], chemParams['kDeg_fSS_i'], Z_ij_sub['fSS'][1]) * Q_fw_i

    fwOut_adv_n = adv.D_advec_water(climate['waterflow_d'][i], Z_ij_sub['fw'][0]) * Q_fw_n
    fwOut_adv_i = adv.D_advec_water(climate['waterflow_d'][i], Z_ij_sub['fw'][1]) * Q_fw_i

    fwIn_adv_n = bgConc['gfreshwc_n'] * climate['waterflow_d'][i]
    fwIn_adv_i = bgConc['gfreshwc_i'] * climate['waterflow_d'][i]

    fSSOut_adv_n = adv.D_advec_susSed(climate['waterflow_d'][i], env['fSSVf'], Z_ij_sub['fSS'][0]) * Q_fw_n
    fSSOut_adv_i = adv.D_advec_susSed(climate['waterflow_d'][i], env['fSSVf'], Z_ij_sub['fSS'][1]) * Q_fw_i

    # fw to other compartments: 1) diffusion (air, n), 2) diffusion to fSedW (n, i)
    fw_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_fw, env['fwA']) * Q_fw_n
    perm_fw_to_fSedW = diff.P_water_sed()
    fw_diff_n_to_fSedW = diff.D_diffu_comp1_comp2(perm_fw_to_fSedW, env['fwA']) * Q_fw_n
    fw_diff_i_to_fSedW = diff.D_diffu_comp1_comp2(perm_fw_to_fSedW, env['fwA']) * Q_fw_i


    # fSS to other subcomparts: 1) deposition (n, i) to fSedS
    k_dep_dry_fSS = adv.k_dep_dry(env['freshssP'], env['freshwP'],
                                  env['dynViscFW'], env['radiusParticlesFW'])

    fSS_dep_n = adv.D_dep_dry(k_dep_dry_fSS, env['freshwA'] * env['fSSVf'], Z_ij_sub['fSS'][0]) * Q_fw_n
    fSS_dep_i = adv.D_dep_dry(k_dep_dry_fSS, env['freshwA'] * env['fSSVf'], Z_ij_sub['fSS'][1]) * Q_fw_i

    # freshwater sediment - fwSed, Q[2]
    Q_fwSed_n = Q[2] * Y_ij['fwSed'][0]
    Q_fwSed_i = Q[2] * Y_ij['fwSed'][1]

    # 1) degradation (n, i)
    fSedW_deg_n = deg.D_deg(env['fSedWV'], chemParams['kDeg_fSedW_n'], Z_ij_sub['fSedW'][0]) * Q_fwSed_n
    fSedW_deg_i = deg.D_deg(env['fSedWV'], chemParams['kDeg_fSedW_i'], Z_ij_sub['fSedW'][1]) * Q_fwSed_i
    fSedS_deg_n = deg.D_deg(env['fSedSV'], chemParams['kDeg_fSedS_n'], Z_ij_sub['fSedS'][0]) * Q_fwSed_n
    fSedS_deg_i = deg.D_deg(env['fSedSV'], chemParams['kDeg_fSedS_i'], Z_ij_sub['fSedS'][1]) * Q_fwSed_i


    # fSedW: diffusion (n, i)
    fSedW_diff_n_to_fw = diff.D_diffu_comp1_comp2(perm_fw_to_fSedW, env['fwA']) * Q_fwSed_n
    fSedW_diff_i_to_fw = diff.D_diffu_comp1_comp2(perm_fw_to_fSedW, env['fwA']) * Q_fwSed_i

    # fSedS: burial (n, i)
    fSedS_burial_n = adv.D_burial(env['burialRateFW'], env['fwA'], Z_ij_sub['fSedS'][0]) * Q_fwSed_n
    fSedS_burial_i = adv.D_burial(env['burialRateFW'], env['fwA'], Z_ij_sub['fSedS'][1]) * Q_fwSed_i

    # fSedS: resuspension (n, i)
    fSedS_resusp_n = adv.D_sedResusp(env['resuspensionRateFW'], env['fwA'], Z_ij['fwSed'][0]) * Q_fwSed_n
    fSedS_resusp_i = adv.D_sedResusp(env['resuspensionRateFW'], env['fwA'], Z_ij['fwSed'][1]) * Q_fwSed_i

    # fSedW and fSedS: advection out (n, i)
    fSed_adv_inflow_n = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_n'])
    fSed_adv_inflow_i = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_i'])
    fSed_adv_outflow_n = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], Z_ij['fwSed'][0]) * Q_fwSed_n
    fSed_adv_outflow_i = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], Z_ij['fwSed'][1]) * Q_fwSed_i


    # seawater - sw, Q[3]
    Q_sw_n = Q[3] * Y_ij['sw'][0]
    Q_sw_i = Q[3] * Y_ij['sw'][1]

    # within sw: 1)degradation (n, i), 2) advection in (n, i), 3) advection out (n, i)
    sw_deg_n = deg.D_deg(env['swV'], chemParams['kDeg_sw_n'], Z_ij_sub['sw'][0]) * Q_sw_n
    sw_deg_i = deg.D_deg(env['swV'], chemParams['kDeg_sw_i'], Z_ij_sub['sw'][1]) * Q_sw_i
    sSS_deg_n = deg.D_deg(env['sSSV'], chemParams['kDeg_sSS_n'], Z_ij_sub['sSS'][0]) * Q_sw_n
    sSS_deg_i = deg.D_deg(env['sSSV'], chemParams['kDeg_sSS_i'], Z_ij_sub['sSS'][1]) * Q_sw_i


    # assume waterflow rate is 10 times of freshwater
    swOut_adv_n = adv.D_advec_water(climate['waterflow_d'][i] * 10, Z_ij_sub['sw'][0]) * Q_sw_n
    swOut_adv_i = adv.D_advec_water(climate['waterflow_d'][i] * 10, Z_ij_sub['sw'][1]) * Q_sw_i
    swIn_adv_n = fwOut_adv_n
    swIn_adv_i = fwOut_adv_i

    sSSOut_adv_n = adv.D_advec_susSed(climate['waterflow_d'][i] * 10, env['sSSVf'],
                                      Z_ij_sub['sSS'][0]) * Q_sw_n
    sSSOut_adv_i = adv.D_advec_susSed(climate['waterflow_d'][i] * 10, env['sSSVf'],
                                      Z_ij_sub['sSS'][1]) * Q_sw_i
    sSSIn_adv_n = fSSOut_adv_n
    sSSIn_adv_i = fSSOut_adv_i

    # sw to other subcomparts: 1) diffusion (air, n), 2) diffusion to sSedW (n, i)
    sw_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_sw, env['swA']) * Q_sw_n
    perm_sw_to_sSedW = diff.P_water_sed()
    sw_diff_n_to_sSedW = diff.D_diffu_comp1_comp2(perm_sw_to_sSedW, env['swA']) * Q_sw_n
    sw_diff_i_to_sSedW = diff.D_diffu_comp1_comp2(perm_sw_to_sSedW, env['swA']) * Q_sw_i


    # sSS: 1) deposition (n, i)
    k_dep_dry_sSS = adv.k_dep_dry(env['seassP'], env['seawP'],
                                  env['dynViscSW'], env['radiusParticlesSW'])

    sSS_dep_n = adv.D_dep_dry(k_dep_dry_sSS, env['seawA'] * env['sSSVf'], Z_ij_sub['sSS'][0]) * Q_sw_n
    sSS_dep_i = adv.D_dep_dry(k_dep_dry_sSS, env['seawA'] * env['sSSVf'], Z_ij_sub['sSS'][1]) * Q_sw_i

    sSS_resusp_n = adv.D_aeroResusp(climate['windspeed_s'][i], env['coastalA'],
                                    chemParams['enrichFactor'], env['seawD'], Z_ij_sub['sSS'][0]) * Q_sw_n
    sSS_resusp_i = adv.D_aeroResusp(climate['windspeed_s'][i], env['coastalA'],
                                    chemParams['enrichFactor'], env['seawD'], Z_ij_sub['sSS'][1]) * Q_sw_i


    # seawater sediment - swSed, Q[4]
    Q_swSed_n = Q[4] * Y_ij['swSed'][0]
    Q_swSed_i = Q[4] * Y_ij['swSed'][1]

    # degradation (n, i)
    sSedW_deg_n = deg.D_deg(env['sSedWV'], chemParams['kDeg_sSedW_n'], Z_ij_sub['sSedW'][0]) * Q_swSed_n
    sSedW_deg_i = deg.D_deg(env['sSedWV'], chemParams['kDeg_sSedW_i'], Z_ij_sub['sSedW'][1]) * Q_swSed_i
    sSedS_deg_n = deg.D_deg(env['sSedSV'], chemParams['kDeg_sSedS_n'], Z_ij_sub['sSedS'][0]) * Q_swSed_n
    sSedS_deg_i = deg.D_deg(env['sSedSV'], chemParams['kDeg_sSedS_i'], Z_ij_sub['sSedS'][1]) * Q_swSed_i


    # sSedW: diffusion (n, i) to sw
    sSedW_diff_n_to_sw = diff.D_diffu_comp1_comp2(perm_sw_to_sSedW, env['swA']) * Q_swSed_n
    sSedW_diff_i_to_sw = diff.D_diffu_comp1_comp2(perm_sw_to_sSedW, env['swA']) * Q_swSed_i

    # sSedS: burial (n, i)
    sSedS_burial_n = adv.D_burial(env['burialRateSW'], env['swA'], Z_ij_sub['sSedS'][0]) * Q_swSed_n
    sSedS_burial_i = adv.D_burial(env['burialRateSW'], env['swA'], Z_ij_sub['sSedS'][1]) * Q_swSed_i

    # sSedS to sSS: resuspension (n, i)
    sSedS_resusp_n = adv.D_sedResusp(env['resuspensionRateSW'], env['swA'], Z_ij['swSed'][0]) * Q_swSed_n
    sSedS_resusp_i = adv.D_sedResusp(env['resuspensionRateSW'], env['swA'], Z_ij['swSed'][1]) * Q_swSed_i

    # sSed: advection out
    sSed_adv_outflow_n = adv.D_advec_water(climate['waterflow_d'][i]*10*env['swadvfrac'], Z_ij['swSed'][0]) * Q_swSed_n
    sSed_adv_outflow_i = adv.D_advec_water(climate['waterflow_d'][i]*10*env['swadvfrac'], Z_ij['swSed'][1]) * Q_swSed_i


    # soil 1 - soil1, Q[5]
    Q_soil1_n = Q[5] * Y_ij['soil1'][0]
    Q_soil1_i = Q[5] * Y_ij['soil1'][1]

    # within soil1: 1) degradation (n)
    soilA1_deg_n = deg.D_deg(env['soilAV1'], chemParams['kDeg_soilA1_n'], Z_ij_sub['soilA1'][0]) * Q_soil1_n
    soilW1_deg_n = deg.D_deg(env['soilWV1'], chemParams['kDeg_soilW1_n'], Z_ij_sub['soilW1'][0]) * Q_soil1_n
    soilW1_deg_i = deg.D_deg(env['soilWV1'], chemParams['kDeg_soilW1_i'], Z_ij_sub['soilW1'][1]) * Q_soil1_i
    soilS1_deg_n = deg.D_deg(env['soilSV1'], chemParams['kDeg_soilS1_n'], Z_ij_sub['soilS1'][0]) * Q_soil1_n
    soilS1_deg_i = deg.D_deg(env['soilSV1'], chemParams['kDeg_soilS1_i'], Z_ij_sub['soilS1'][1]) * Q_soil1_i


    # soil1 to other compartments: 1) water runoff (n, i) to fw, 2) infiltration (n, i) to deep soil deepS1,
    # 3) diffusion (n) to air, 4) soil erosion to fSS (n, i)
    soil1_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_soil1, env['soilA1']) * Q_soil1_n
    soilW1_runoff_n = adv.D_runoff(climate['precip_mm'][i], env['CN1'], env['soilA1'], Z_water=1.0) * Q_soil1_n
    soilW1_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN1'], env['soilA1'], Z_water=1.0) * Q_soil1_i

    D_infil_1, k_infil_1 = adv.D_infiltra(climate['precip_mm'][i], env['CN1'], climate['evap_mm'][i], env['FC1'],
                                              env['soilWC1'],env['soilV1'], env['soilA1'], Z_water=1.0)

    soilW1_infil_n = D_infil_1 * Q_soil1_n
    soilW1_infil_i = D_infil_1 * Q_soil1_i



    soilS1_erosion_n = adv.D_erosion(climate['precip_mm'][i], env['slope1'],
                                     env['Kfact1'], env['cropManageFactor1'],
                                     env['supportFactor1'], env['soilA1'],
                                     env['soilP1'], Z_ij_sub['soilS1'][0]) * Q_soil1_n
    soilS1_erosion_i = adv.D_erosion(climate['precip_mm'][i], env['slope1'],
                                     env['Kfact1'], env['cropManageFactor1'],
                                     env['supportFactor1'], env['soilA1'],
                                     env['soilP1'], Z_ij_sub['soilS1'][1]) * Q_soil1_i

    soilS1_windErosion_n = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness1'],
                                             env['Kconstant1'], env['airP'], env['soilA1'], env['A1'],
                                             env['TSV1'], env['TSVmin1'], env['z_wind1'], env['percWind1'],
                                             env['windConstant1'], env['percUncovered1'], env['percSuspended1'],
                                             env['soilP1'], Z_ij_sub['soilS1'][0]) * Q_soil1_n

    soilS1_windErosion_i = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness1'],
                                             env['Kconstant1'], env['airP'], env['soilA1'], env['A1'],
                                             env['TSV1'], env['TSVmin1'], env['z_wind1'], env['percWind1'],
                                             env['windConstant1'], env['percUncovered1'], env['percSuspended1'],
                                             env['soilP1'], Z_ij_sub['soilS1'][1]) * Q_soil1_i


    # deep soil 1 - deepS1, Q[6]
    Q_deepS1_n = Q[6] * Y_ij['deepS1'][0]
    Q_deepS1_i = Q[6] * Y_ij['deepS1'][1]

    # within deepS1: 1) degradation (n, i)
    deepS1_deg_n = deg.D_deg(env['deepSV1'], chemParams['kDeg_deepS1_n'], Z_ij['deepS1'][0]) * Q_deepS1_n
    deepS1_deg_i = deg.D_deg(env['deepSV1'], chemParams['kDeg_deepS1_i'], Z_ij['deepS1'][1]) * Q_deepS1_i

    # deepS1 to other compartments: 1) leaching (n, i) to fw
    deepS1_leach_n = adv.D_leach(k_infil_1, Z_water=1.0) * Q_deepS1_n
    deepS1_leach_i = adv.D_leach(k_infil_1, Z_water=1.0) * Q_deepS1_i


    # soil 2 - soil2, Q[7]
    Q_soil2_n = Q[7] * Y_ij['soil2'][0]
    Q_soil2_i = Q[7] * Y_ij['soil2'][1]

    # within soil2: 1) degradation (n)
    soilA2_deg_n = deg.D_deg(env['soilAV2'], chemParams['kDeg_soilA2_n'], Z_ij_sub['soilA2'][0]) * Q_soil2_n
    soilW2_deg_n = deg.D_deg(env['soilWV2'], chemParams['kDeg_soilW2_n'], Z_ij_sub['soilW2'][0]) * Q_soil2_n
    soilW2_deg_i = deg.D_deg(env['soilWV2'], chemParams['kDeg_soilW2_i'], Z_ij_sub['soilW2'][1]) * Q_soil2_i
    soilS2_deg_n = deg.D_deg(env['soilSV2'], chemParams['kDeg_soilS2_n'], Z_ij_sub['soilS2'][0]) * Q_soil2_n
    soilS2_deg_i = deg.D_deg(env['soilSV2'], chemParams['kDeg_soilS2_i'], Z_ij_sub['soilS2'][1]) * Q_soil2_i


    # soilW2 to other subcomparts: 1) water runoff (n, i) to fw, 2) infiltration (n, i) to deep soil deepS2,
    # 3) diffusion (n) to air, 4) erosion (n, i) to fSS
    soil2_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_soil1, env['soilA2']) * Q_soil2_n

    soilW2_runoff_n = adv.D_runoff(climate['precip_mm'][i], env['CN2'], env['soilA2'], Z_water=1.0) * Q_soil2_n
    soilW2_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN2'], env['soilA2'], Z_water=1.0) * Q_soil2_i

    D_infil_2, k_infil_2 = adv.D_infiltra(climate['precip_mm'][i], env['CN2'], climate['evap_mm'][i], env['FC2'],
                                              env['soilWC2'],env['soilV2'], env['soilA2'], Z_water=1.0)

    soilW2_infil_n = D_infil_2 * Q_soil2_n
    soilW2_infil_i = D_infil_2 * Q_soil2_i


    soilS2_erosion_n = adv.D_erosion(climate['precip_mm'][i], env['slope2'],
                                     env['Kfact2'], env['cropManageFactor2'],
                                     env['supportFactor2'], env['soilA2'],
                                     env['soilP2'], Z_ij_sub['soilS2'][0]) * Q_soil2_n
    soilS2_erosion_i = adv.D_erosion(climate['precip_mm'][i], env['slope2'],
                                     env['Kfact2'], env['cropManageFactor2'],
                                     env['supportFactor2'], env['soilA2'],
                                     env['soilP2'], Z_ij_sub['soilS2'][1]) * Q_soil2_i

    soilS2_windErosion_n = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness2'],
                                             env['Kconstant2'], env['airP'], env['soilA2'], env['A2'],
                                             env['TSV2'], env['TSVmin2'], env['z_wind2'], env['percWind2'],
                                             env['windConstant2'], env['percUncovered2'], env['percSuspended2'],
                                             env['soilP2'], Z_ij_sub['soilS2'][0]) * Q_soil2_n

    soilS2_windErosion_i = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness2'],
                                             env['Kconstant2'], env['airP'], env['soilA2'], env['A2'],
                                             env['TSV2'], env['TSVmin2'], env['z_wind2'], env['percWind2'],
                                             env['windConstant2'], env['percUncovered2'], env['percSuspended2'],
                                             env['soilP2'], Z_ij_sub['soilS2'][1]) * Q_soil2_i


    # deep soil 2 - deepS2, Q[8]
    Q_deepS2_n = Q[8] * Y_ij['deepS2'][0]
    Q_deepS2_i = Q[8] * Y_ij['deepS2'][1]

    # within deepS2: 1) degradation (n, i)
    deepS2_deg_n = deg.D_deg(env['deepSV2'], chemParams['kDeg_deepS2_n'], Z_ij['deepS2'][0]) * Q_deepS2_n
    deepS2_deg_i = deg.D_deg(env['deepSV2'], chemParams['kDeg_deepS2_i'], Z_ij['deepS2'][1]) * Q_deepS2_i
    # deepS2 to other compartments: 1) leaching (n, i) to fw
    deepS2_leach_n = adv.D_leach(k_infil_2, Z_water=1.0) * Q_deepS2_n
    deepS2_leach_i = adv.D_leach(k_infil_2, Z_water=1.0) * Q_deepS2_i


    # soil 3 - soil3, Q[9]
    Q_soil3_n = Q[9] * Y_ij['soil3'][0]
    Q_soil3_i = Q[9] * Y_ij['soil3'][1]

    # within soil3: 1) degradation (n)
    soilA3_deg_n = deg.D_deg(env['soilAV3'], chemParams['kDeg_soilA3_n'], Z_ij_sub['soilA3'][0]) * Q_soil3_n
    soilW3_deg_n = deg.D_deg(env['soilWV3'], chemParams['kDeg_soilW3_n'], Z_ij_sub['soilW3'][0]) * Q_soil3_n
    soilW3_deg_i = deg.D_deg(env['soilWV3'], chemParams['kDeg_soilW3_i'], Z_ij_sub['soilW3'][1]) * Q_soil3_i
    soilS3_deg_n = deg.D_deg(env['soilSV3'], chemParams['kDeg_soilS3_n'], Z_ij_sub['soilS3'][0]) * Q_soil3_n
    soilS3_deg_i = deg.D_deg(env['soilSV3'], chemParams['kDeg_soilS3_i'], Z_ij_sub['soilS3'][1]) * Q_soil3_i


    # soilW3 to other subcomparts: 1) water runoff (n, i) to fw, 2) infiltration (n, i) to deep soil deepS3,
    # 3) diffusion (n) to air, 4) erosion (n, i) to fSS
    # soil3_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_soil3, env['soilA3']) * 1.51e-4 * Q_soil3_n
    soil3_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_soil3, env['soilA3']) * Q_soil3_n

    soilW3_runoff_n = adv.D_runoff(climate['precip_mm'][i], env['CN3'], env['soilA3'], Z_water=1.0) * Q_soil3_n
    soilW3_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN3'], env['soilA3'], Z_water=1.0) * Q_soil3_i

    D_infil_3, k_infil_3 = adv.D_infiltra(climate['precip_mm'][i], env['CN3'], climate['evap_mm'][i], env['FC3'],
                                              env['soilWC3'],env['soilV3'], env['soilA3'], Z_water=1.0)

    soilW3_infil_n = D_infil_3 * Q_soil3_n
    soilW3_infil_i = D_infil_3 * Q_soil3_i


    soilS3_erosion_n = adv.D_erosion(climate['precip_mm'][i], env['slope3'],
                                     env['Kfact3'], env['cropManageFactor3'],
                                     env['supportFactor3'], env['soilA3'],
                                     env['soilP3'], Z_ij_sub['soilS3'][0]) * Q_soil3_n
    soilS3_erosion_i = adv.D_erosion(climate['precip_mm'][i], env['slope3'],
                                     env['Kfact3'], env['cropManageFactor3'],
                                     env['supportFactor3'], env['soilA3'],
                                     env['soilP3'], Z_ij_sub['soilS3'][1]) * Q_soil3_i

    soilS3_windErosion_n = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness3'],
                                             env['Kconstant3'], env['airP'], env['soilA3'],
                                             env['A3'],
                                             env['TSV3'], env['TSVmin3'], env['z_wind3'],
                                             env['percWind3'],
                                             env['windConstant3'], env['percUncovered3'],
                                             env['percSuspended3'],
                                             env['soilP3'], Z_ij_sub['soilS3'][0]) * Q_soil3_n

    soilS3_windErosion_i = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness3'],
                                             env['Kconstant3'], env['airP'], env['soilA3'],
                                             env['A3'],
                                             env['TSV3'], env['TSVmin3'], env['z_wind3'],
                                             env['percWind3'],
                                             env['windConstant3'], env['percUncovered3'],
                                             env['percSuspended3'],
                                             env['soilP3'], Z_ij_sub['soilS3'][1]) * Q_soil3_i


    # deep soil 3 - deepS3, Q[10]
    Q_deepS3_n = Q[10] * Y_ij['deepS3'][0]
    Q_deepS3_i = Q[10] * Y_ij['deepS3'][1]

    # within deepS3: 1) degradation (n, i)
    deepS3_deg_n = deg.D_deg(env['deepSV3'], chemParams['kDeg_deepS3_n'], Z_ij['deepS3'][0]) * Q_deepS3_n
    deepS3_deg_i = deg.D_deg(env['deepSV3'], chemParams['kDeg_deepS3_i'], Z_ij['deepS3'][1]) * Q_deepS3_i
    # deepS3 to other compartment: 1) leaching (n, i) to fw
    deepS3_leach_n = adv.D_leach(k_infil_3, Z_water=1.0) * Q_deepS3_n
    deepS3_leach_i = adv.D_leach(k_infil_3, Z_water=1.0) * Q_deepS3_i


    # soil 4 - soil4, Q[11]
    Q_soil4_n = Q[11] * Y_ij['soil4'][0]
    Q_soil4_i = Q[11] * Y_ij['soil4'][1]

    # within soil4: 1) degradation (n)
    soilA4_deg_n = deg.D_deg(env['soilAV4'], chemParams['kDeg_soilA4_n'], Z_ij_sub['soilA4'][0]) * Q_soil4_n
    soilW4_deg_n = deg.D_deg(env['soilWV4'], chemParams['kDeg_soilW4_n'], Z_ij_sub['soilW4'][0]) * Q_soil4_n
    soilW4_deg_i = deg.D_deg(env['soilWV4'], chemParams['kDeg_soilW4_i'], Z_ij_sub['soilW4'][1]) * Q_soil4_i
    soilS4_deg_n = deg.D_deg(env['soilSV4'], chemParams['kDeg_soilS4_n'], Z_ij_sub['soilS4'][0]) * Q_soil4_n
    soilS4_deg_i = deg.D_deg(env['soilSV4'], chemParams['kDeg_soilS4_i'], Z_ij_sub['soilS4'][1]) * Q_soil4_i

    # soilW4 to other subcomparts: 1) water runoff (n, i) to fw, 2) infiltration (n, i) to deep soil deepS4,
    # 3) diffusion (n) to air, 4) erosion (n, i) to fSS
    soil4_diff_n_to_air = diff.D_diffu_comp1_comp2(perm_air_to_soil4, env['soilA4']) * Q_soil4_n

    soilW4_runoff_n = adv.D_runoff(climate['precip_mm'][i], env['CN4'], env['soilA4'], Z_water=1.0) * Q_soil4_n
    soilW4_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN4'], env['soilA4'], Z_water=1.0) * Q_soil4_i

    D_infil_4, k_infil_4 = adv.D_infiltra(climate['precip_mm'][i], env['CN4'], climate['evap_mm'][i], env['FC4'],
                                              env['soilWC4'],env['soilV4'], env['soilA4'], Z_water=1.0)

    soilW4_infil_n = D_infil_4 * Q_soil4_n
    soilW4_infil_i = D_infil_4 * Q_soil4_i


    soilS4_erosion_n = adv.D_erosion(climate['precip_mm'][i], env['slope4'],
                                     env['Kfact4'], env['cropManageFactor4'],
                                     env['supportFactor4'], env['soilA4'],
                                     env['soilP4'], Z_ij_sub['soilS4'][0]) * Q_soil4_n
    soilS4_erosion_i = adv.D_erosion(climate['precip_mm'][i], env['slope4'],
                                     env['Kfact4'], env['cropManageFactor4'],
                                     env['supportFactor4'], env['soilA4'],
                                     env['soilP4'], Z_ij_sub['soilS4'][1]) * Q_soil4_i

    soilS4_windErosion_n = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness4'],
                                             env['Kconstant4'], env['airP'], env['soilA4'],
                                             env['A4'],
                                             env['TSV4'], env['TSVmin4'], env['z_wind4'],
                                             env['percWind4'],
                                             env['windConstant4'], env['percUncovered4'],
                                             env['percSuspended4'],
                                             env['soilP4'], Z_ij_sub['soilS4'][0]) * Q_soil4_n

    soilS4_windErosion_i = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness4'],
                                             env['Kconstant4'], env['airP'], env['soilA4'],
                                             env['A4'],
                                             env['TSV4'], env['TSVmin4'], env['z_wind4'],
                                             env['percWind4'],
                                             env['windConstant4'], env['percUncovered4'],
                                             env['percSuspended4'],
                                             env['soilP4'], Z_ij_sub['soilS4'][1]) * Q_soil4_i


    # deep soil 4 - deepS4, Q[12]
    Q_deepS4_n = Q[12] * Y_ij['deepS4'][0]
    Q_deepS4_i = Q[12] * Y_ij['deepS4'][1]

    # within deepS4: 1) degradation (n, i)
    deepS4_deg_n = deg.D_deg(env['deepSV4'], chemParams['kDeg_deepS4_n'], Z_ij['deepS4'][0]) * Q_deepS4_n
    deepS4_deg_i = deg.D_deg(env['deepSV4'], chemParams['kDeg_deepS4_i'], Z_ij['deepS4'][1]) * Q_deepS4_i
    # deepS4 to other compartment: 1) leaching (n, i) to fw
    deepS4_leach_n = adv.D_leach(k_infil_4, Z_water=1.0) * Q_deepS4_n
    deepS4_leach_i = adv.D_leach(k_infil_4, Z_water=1.0) * Q_deepS4_i

    ###################################################################
    # processes output to transport rate kg/day
    # N * molar mass = mol/day * kg/mol = kg/day
    ###################################################################

    # 1) degradation process
    deg_air = (air_deg_n + aer_deg_n + aer_deg_i) * chemParams['molar_mass']
    deg_fw = (fw_deg_n + fw_deg_i + fSS_deg_n + fSS_deg_i) * chemParams['molar_mass']
    deg_fwSed = (fSedW_deg_n + fSedW_deg_i + fSedS_deg_n + fSedS_deg_i) * chemParams['molar_mass']
    deg_sw = (sw_deg_n + sw_deg_i + sSS_deg_n + sSS_deg_i) * chemParams['molar_mass']
    deg_swSed = (sSedW_deg_n + sSedW_deg_i + sSedS_deg_n + sSedS_deg_i) * chemParams['molar_mass']
    deg_soil1 = (soilA1_deg_n + soilW1_deg_n + soilW1_deg_i + soilS1_deg_n + soilS1_deg_i) * chemParams['molar_mass']
    deg_deepS1 = (deepS1_deg_n + deepS1_deg_i) * chemParams['molar_mass']
    deg_soil2 = (soilA2_deg_n + soilW2_deg_n + soilW2_deg_i + soilS2_deg_n + soilS2_deg_i) * chemParams['molar_mass']
    deg_deepS2 = (deepS2_deg_n + deepS2_deg_i) * chemParams['molar_mass']
    deg_soil3 = (soilA3_deg_n + soilW3_deg_n + soilW3_deg_i + soilS3_deg_n + soilS3_deg_i) * chemParams['molar_mass']
    deg_deepS3 = (deepS3_deg_n + deepS3_deg_i) * chemParams['molar_mass']
    deg_soil4 = (soilA4_deg_n + soilW4_deg_n + soilW4_deg_i + soilS4_deg_n + soilS4_deg_i) * chemParams['molar_mass']
    deg_deepS4 = (deepS4_deg_n + deepS4_deg_i) * chemParams['molar_mass']

    # 2) advection process
    adv_air_in = (airIn_adv_n + airIn_adv_i) * chemParams['molar_mass']
    adv_air_out = (airOut_adv_n + aerOut_adv_n + aerOut_adv_i) * chemParams['molar_mass']
    adv_fw_in = (fwIn_adv_n + fwIn_adv_i) * chemParams['molar_mass']
    adv_fw_out = (fwOut_adv_n + fwOut_adv_i + fSSOut_adv_n + fSSOut_adv_i) * chemParams['molar_mass']
    adv_fwSed_in = (fSed_adv_inflow_n + fSed_adv_inflow_i) * chemParams['molar_mass']
    adv_fwSed_out = (fSed_adv_outflow_n + fSed_adv_outflow_i) * chemParams['molar_mass']
    adv_sw_in = (swIn_adv_n + swIn_adv_i + sSSIn_adv_n + sSSIn_adv_i) * chemParams['molar_mass']
    adv_sw_out = (swOut_adv_n + swOut_adv_i + sSSOut_adv_n + sSSOut_adv_i) * chemParams['molar_mass']
    adv_swSed_out = (sSed_adv_outflow_n + sSed_adv_outflow_i) * chemParams['molar_mass']

    # 3) deposition process
    dep_dry_air = (aer_dep_dry_n + aer_dep_dry_i) * chemParams['molar_mass']
    dep_dry_air_fw = (aer_dep_dry_n_to_fSS + aer_dep_dry_i_to_fSS) * chemParams['molar_mass']
    dep_dry_air_sw = (aer_dep_dry_n_to_sSS + aer_dep_dry_i_to_sSS) * chemParams['molar_mass']
    dep_dry_air_soil1 = (aer_dep_dry_n_to_soil1 + aer_dep_dry_i_to_soil1) * chemParams['molar_mass']
    dep_dry_air_soil2 = (aer_dep_dry_n_to_soil2 + aer_dep_dry_i_to_soil2) * chemParams['molar_mass']
    dep_dry_air_soil3 = (aer_dep_dry_n_to_soil3 + aer_dep_dry_i_to_soil3) * chemParams['molar_mass']
    dep_dry_air_soil4 = (aer_dep_dry_n_to_soil4 + aer_dep_dry_i_to_soil4) * chemParams['molar_mass']
    
    dep_wet_air = (aer_dep_wet_n + aer_dep_wet_i) * chemParams['molar_mass']
    dep_wet_air_fw = (aer_dep_wet_n_to_fSS + aer_dep_wet_i_to_fSS) * chemParams['molar_mass']
    dep_wet_air_sw = (aer_dep_wet_n_to_sSS + aer_dep_wet_i_to_sSS) * chemParams['molar_mass']
    dep_wet_air_soil1 = (aer_dep_wet_n_to_soil1 + aer_dep_wet_i_to_soil1) * chemParams['molar_mass']
    dep_wet_air_soil2 = (aer_dep_wet_n_to_soil2 + aer_dep_wet_i_to_soil2) * chemParams['molar_mass']
    dep_wet_air_soil3 = (aer_dep_wet_n_to_soil3 + aer_dep_wet_i_to_soil3) * chemParams['molar_mass']
    dep_wet_air_soil4 = (aer_dep_wet_n_to_soil4 + aer_dep_wet_i_to_soil4) * chemParams['molar_mass']

    rain_dis_air = air_rain_diss_n * chemParams['molar_mass']
    rain_dis_air_fw = air_rain_diss_n_to_fw * chemParams['molar_mass']
    rain_dis_air_sw = air_rain_diss_n_to_sw * chemParams['molar_mass']
    rain_dis_air_soil1 = air_rain_diss_n_to_soil1 * chemParams['molar_mass']
    rain_dis_air_soil2 = air_rain_diss_n_to_soil2 * chemParams['molar_mass']
    rain_dis_air_soil3 = air_rain_diss_n_to_soil3 * chemParams['molar_mass']
    rain_dis_air_soil4 = air_rain_diss_n_to_soil4 * chemParams['molar_mass']

    dep_fSS = (fSS_dep_n + fSS_dep_i) * chemParams['molar_mass']
    dep_sSS = (sSS_dep_n + sSS_dep_i) * chemParams['molar_mass']

    # 4) diffusion process
    diff_air_fw = air_diff_n_to_fw * chemParams['molar_mass']
    diff_air_sw = air_diff_n_to_sw * chemParams['molar_mass']
    diff_fw_fSedW = (fw_diff_n_to_fSedW + fw_diff_i_to_fSedW) * chemParams['molar_mass']
    diff_sw_sSedW = (sw_diff_n_to_sSedW + sw_diff_i_to_sSedW) * chemParams['molar_mass']
    diff_air_soil1 = air_diff_n_to_soil1 * chemParams['molar_mass']
    diff_air_soil2 = air_diff_n_to_soil2 * chemParams['molar_mass']
    diff_air_soil3 = air_diff_n_to_soil3 * chemParams['molar_mass']
    diff_air_soil4 = air_diff_n_to_soil4 * chemParams['molar_mass']

    diff_fw_air = fw_diff_n_to_air * chemParams['molar_mass']
    diff_sw_air = sw_diff_n_to_air * chemParams['molar_mass']
    diff_fSedW_fw = (fSedW_diff_n_to_fw + fSedW_diff_i_to_fw) * chemParams['molar_mass']
    diff_sSedW_sw = (sSedW_diff_n_to_sw + sSedW_diff_i_to_sw) * chemParams['molar_mass']
    diff_soil1_air = soil1_diff_n_to_air * chemParams['molar_mass']
    diff_soil2_air = soil2_diff_n_to_air * chemParams['molar_mass']
    diff_soil3_air = soil3_diff_n_to_air * chemParams['molar_mass']
    diff_soil4_air = soil4_diff_n_to_air * chemParams['molar_mass']

    # 5) other process
    burial_fwSed = (fSedS_burial_n + fSedS_burial_i) * chemParams['molar_mass']
    burial_swSed = (sSedS_burial_n + sSedS_burial_i) * chemParams['molar_mass']
    resusp_fwSed = (fSedS_resusp_n + fSedS_resusp_i) * chemParams['molar_mass']
    resusp_swSed = (sSedS_resusp_n + sSedS_resusp_i) * chemParams['molar_mass']
    aero_resusp_sSS = (sSS_resusp_n + sSS_resusp_i) * chemParams['molar_mass']
    runoff_soil1 = (soilW1_runoff_n + soilW1_runoff_i) * chemParams['molar_mass']
    runoff_soil2 = (soilW2_runoff_n + soilW2_runoff_i) * chemParams['molar_mass']
    runoff_soil3 = (soilW3_runoff_n + soilW3_runoff_i) * chemParams['molar_mass']
    runoff_soil4 = (soilW4_runoff_n + soilW4_runoff_i) * chemParams['molar_mass']
    erosion_soil1 = (soilS1_erosion_n + soilS1_erosion_i) * chemParams['molar_mass']
    erosion_soil2 = (soilS2_erosion_n + soilS2_erosion_i) * chemParams['molar_mass']
    erosion_soil3 = (soilS3_erosion_n + soilS3_erosion_i) * chemParams['molar_mass']
    erosion_soil4 = (soilS4_erosion_n + soilS4_erosion_i) * chemParams['molar_mass']
    wind_erosion_soil1 = (soilS1_windErosion_n + soilS1_windErosion_i) * chemParams['molar_mass']
    wind_erosion_soil2 = (soilS2_windErosion_n + soilS2_windErosion_i) * chemParams['molar_mass']
    wind_erosion_soil3 = (soilS3_windErosion_n + soilS3_windErosion_i) * chemParams['molar_mass']
    wind_erosion_soil4 = (soilS4_windErosion_n + soilS4_windErosion_i) * chemParams['molar_mass']
    infiltra_soil1 = (soilW1_infil_n + soilW1_infil_i) * chemParams['molar_mass']
    infiltra_soil2 = (soilW2_infil_n + soilW2_infil_i) * chemParams['molar_mass']
    infiltra_soil3 = (soilW3_infil_n + soilW3_infil_i) * chemParams['molar_mass']
    infiltra_soil4 = (soilW4_infil_n + soilW4_infil_i) * chemParams['molar_mass']
    leach_soil1 = (deepS1_leach_n + deepS1_leach_i) * chemParams['molar_mass']
    leach_soil2 = (deepS2_leach_n + deepS2_leach_i) * chemParams['molar_mass']
    leach_soil3 = (deepS3_leach_n + deepS3_leach_i) * chemParams['molar_mass']
    leach_soil4 = (deepS4_leach_n + deepS4_leach_i) * chemParams['molar_mass']

    processes = [adv_air_in, adv_air_out, adv_fw_in, adv_fw_out, adv_fwSed_in, adv_fwSed_out, adv_sw_in, adv_sw_out,
                 adv_swSed_out, dep_dry_air, dep_dry_air_fw, dep_dry_air_sw, dep_dry_air_soil1, dep_dry_air_soil2,
                 dep_dry_air_soil3, dep_dry_air_soil4, dep_wet_air, dep_wet_air_fw, dep_wet_air_sw, dep_wet_air_soil1,
                 dep_wet_air_soil2, dep_wet_air_soil3, dep_wet_air_soil4, rain_dis_air, rain_dis_air_fw, rain_dis_air_sw,
                 rain_dis_air_soil1, rain_dis_air_soil2, rain_dis_air_soil3, rain_dis_air_soil4, dep_fSS, dep_sSS,
                 diff_air_fw, diff_air_sw, diff_fw_fSedW, diff_sw_sSedW, diff_air_soil1, diff_air_soil2, diff_air_soil3,
                 diff_air_soil4, diff_fw_air, diff_sw_air, diff_fSedW_fw, diff_sSedW_sw, diff_soil1_air, diff_soil2_air,
                 diff_soil3_air, diff_soil4_air, burial_fwSed, burial_swSed, resusp_fwSed, resusp_swSed, aero_resusp_sSS,
                 runoff_soil1, runoff_soil2, runoff_soil3, runoff_soil4, erosion_soil1, erosion_soil2, erosion_soil3,
                 erosion_soil4, wind_erosion_soil1, wind_erosion_soil2, wind_erosion_soil3, wind_erosion_soil4,
                 infiltra_soil1, infiltra_soil2, infiltra_soil3, infiltra_soil4, leach_soil1, leach_soil2, leach_soil3,
                 leach_soil4, deg_air, deg_fw, deg_fwSed, deg_sw, deg_swSed, deg_soil1, deg_deepS1, deg_soil2, deg_deepS2,
                 deg_soil3, deg_deepS3, deg_soil4, deg_deepS4]

    return processes

