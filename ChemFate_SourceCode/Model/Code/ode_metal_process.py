from advective_processes import AdvectiveProcess
from diffusion_process_ion import Diffusion


def metal_process(Q, i, env, chemParams, climate, bgConc, Z_ij, Y_ij):
    # Q is the total aquivalence by each subcompartment, also Q_t, unit: mol/m3
    # Q_n = Q_t * Y_n, and Q_i = Q_t * Y_i,
    # D unit: m3/day, N unit: mol/day, Q unit: mol/m3
    # N = D*Q

    adv = AdvectiveProcess()
    diff = Diffusion(chemParams['molar_mass'], chemParams['Kaw_n'], climate['windspeed_d'][i])

    # Air - air, Q[0]
    Q_air_p = Q[0] * Y_ij['air'][0]
    Q_air_c = Q[0] * Y_ij['air'][1]
    Q_air_i = Q[0] * Y_ij['air'][2]

    # within air: 1) advection in (p, c, i), 3) advection out (p, c, i)
    aerOut_adv_p = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'],
                                   env['airH'], env['aerVf'],
                                   Z_ij['air'][0]) * Q_air_p
    aerOut_adv_c = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'],
                                   env['airH'], env['aerVf'],
                                   Z_ij['air'][1]) * Q_air_c
    aerOut_adv_i = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'],
                                   env['airH'], env['aerVf'],
                                   Z_ij['air'][2]) * Q_air_i

    airIn_adv_p = bgConc['gairc_p'] * adv.G_advec_air(climate['windspeed_d'][i],
                                                      env['airA'], env['airH'])
    airIn_adv_i = bgConc['gairc_i'] * adv.G_advec_air(climate['windspeed_d'][i],
                                                      env['airA'], env['airH'])

    # air to other compartments: 1) dry deposition (p, i), 2) wet deposition (p, i),

    ## aerosol dry deposition
    k_dep_aer_dry = adv.k_dep_dry(env['aerP'], env['airP'], env['dynViscAir'], env['radiusParticlesAer'])
    aer_dep_dry_p = adv.D_dep_dry(k_dep_aer_dry, env['airA'] * env['aerVf'], Z_ij['air'][0]) * Q_air_p
    aer_dep_dry_i = adv.D_dep_dry(k_dep_aer_dry, env['airA'] * env['aerVf'], Z_ij['air'][2]) * Q_air_i
    aer_dep_dry_p_to_fSS = aer_dep_dry_p * (env['fwA'] / env['airA'])
    aer_dep_dry_p_to_sSS = aer_dep_dry_p * (env['swA'] / env['airA'])
    aer_dep_dry_p_to_soil1 = aer_dep_dry_p * (env['soilA1'] / env['airA'])
    aer_dep_dry_p_to_soil2 = aer_dep_dry_p * (env['soilA2'] / env['airA'])
    aer_dep_dry_p_to_soil3 = aer_dep_dry_p * (env['soilA3'] / env['airA'])
    aer_dep_dry_p_to_soil4 = aer_dep_dry_p * (env['soilA4'] / env['airA'])
    aer_dep_dry_i_to_fSS = aer_dep_dry_i * (env['fwA'] / env['airA'])
    aer_dep_dry_i_to_sSS = aer_dep_dry_i * (env['swA'] / env['airA'])
    aer_dep_dry_i_to_soil1 = aer_dep_dry_i * (env['soilA1'] / env['airA'])
    aer_dep_dry_i_to_soil2 = aer_dep_dry_i * (env['soilA2'] / env['airA'])
    aer_dep_dry_i_to_soil3 = aer_dep_dry_i * (env['soilA3'] / env['airA'])
    aer_dep_dry_i_to_soil4 = aer_dep_dry_i * (env['soilA4'] / env['airA'])

    ## aerosol wet deposition
    aer_dep_wet_p = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'],
                                  env['airA'] * env['aerVf'], Z_ij['air'][0]) * Q_air_p
    aer_dep_wet_i = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'],
                                  env['airA'] * env['aerVf'], Z_ij['air'][2]) * Q_air_i

    aer_dep_wet_p_to_fSS = aer_dep_wet_p * (env['fwA'] / env['airA'])
    aer_dep_wet_p_to_sSS = aer_dep_wet_p * (env['swA'] / env['airA'])
    aer_dep_wet_p_to_soil1 = aer_dep_wet_p * (env['soilA1'] / env['airA'])
    aer_dep_wet_p_to_soil2 = aer_dep_wet_p * (env['soilA2'] / env['airA'])
    aer_dep_wet_p_to_soil3 = aer_dep_wet_p * (env['soilA3'] / env['airA'])
    aer_dep_wet_p_to_soil4 = aer_dep_wet_p * (env['soilA4'] / env['airA'])

    aer_dep_wet_i_to_fSS = aer_dep_wet_i * (env['fwA'] / env['airA'])
    aer_dep_wet_i_to_sSS = aer_dep_wet_i * (env['swA'] / env['airA'])
    aer_dep_wet_i_to_soil1 = aer_dep_wet_i * (env['soilA1'] / env['airA'])
    aer_dep_wet_i_to_soil2 = aer_dep_wet_i * (env['soilA2'] / env['airA'])
    aer_dep_wet_i_to_soil3 = aer_dep_wet_i * (env['soilA3'] / env['airA'])
    aer_dep_wet_i_to_soil4 = aer_dep_wet_i * (env['soilA4'] / env['airA'])


    # freshwater water - fw, Q[1]
    Q_fw_p = Q[1] * Y_ij['fw'][0]
    Q_fw_c = Q[1] * Y_ij['fw'][1]
    Q_fw_i = Q[1] * Y_ij['fw'][2]

    # within fw: 1) advection in (n, i), 2) advection out (n, i)
    fwOut_adv_c = adv.D_advec_water(climate['waterflow_d'][i], Z_ij['fw'][1]) * Q_fw_c
    fwOut_adv_i = adv.D_advec_water(climate['waterflow_d'][i], Z_ij['fw'][2]) * Q_fw_i
    fwIn_adv_p = bgConc['gfreshwc_p'] * climate['waterflow_d'][i]
    fwIn_adv_c = bgConc['gfreshwc_c'] * climate['waterflow_d'][i]
    fwIn_adv_i = bgConc['gfreshwc_i'] * climate['waterflow_d'][i]
    fSSOut_adv_p = adv.D_advec_susSed(climate['waterflow_d'][i], env['fSSVf'], Z_ij['fw'][0]) * Q_fw_p

    # fw to other compartments: 1) diffusion to fSedW (c, i)
    fw_diff_c_to_fSedW = diff.D_diffu_comp1_comp2(7.0*10**(-7), env['fwA']*0.6) * Q_fw_c
    fw_diff_i_to_fSedW = diff.D_diffu_comp1_comp2(1.0*10**(-5), env['fwA']*0.6) * Q_fw_i

    # fSS to other subcomparts: 1) deposition (p) to fSedS
    k_dep_dry_fSS = adv.k_dep_dry(env['freshssP'], env['freshwP'],
                                  env['dynViscFW'], env['radiusParticlesFW'])
    fSS_dep_p = adv.D_dep_dry(k_dep_dry_fSS, env['freshwA'], Z_ij['fw'][0]) * Q_fw_p


    # freshwater sediment - fwSed, Q[2]
    Q_fwSed_p = Q[2] * Y_ij['fwSed'][0]
    Q_fwSed_c = Q[2] * Y_ij['fwSed'][1]
    Q_fwSed_i = Q[2] * Y_ij['fwSed'][2]

    # fSedW: diffusion (n, i)
    fSedW_diff_c_to_fw = diff.D_diffu_comp1_comp2(7.0*10**(-7), env['fwA']*0.6) * Q_fwSed_c
    fSedW_diff_i_to_fw = diff.D_diffu_comp1_comp2(1.0*10**(-5), env['fwA']*0.6) * Q_fwSed_i

    # fSedS: burial (p)
    fSedS_burial_p = adv.D_burial(env['burialRateFW'], env['fwA']*0.6, Z_ij['fwSed'][0]) * Q_fwSed_p

    # fSedS: resuspension (p)
    fSedS_resusp_p = adv.D_sedResusp(env['resuspensionRateFW'], env['fwA']*0.6, Z_ij['fwSed'][0]) * Q_fwSed_p

    # fSedW and fSedS: advection out (p, c, i)
    fSed_adv_inflow_p = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_p'])
    fSed_adv_inflow_c = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_c'])
    fSed_adv_inflow_i = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_i'])
    fSed_adv_outflow_p = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], Z_ij['fwSed'][0]) * Q_fwSed_p
    fSed_adv_outflow_c = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], Z_ij['fwSed'][1]) * Q_fwSed_c
    fSed_adv_outflow_i = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], Z_ij['fwSed'][2]) * Q_fwSed_i

    # seawater - sw, Q[3]
    Q_sw_p = Q[3] * Y_ij['sw'][0]
    Q_sw_c = Q[3] * Y_ij['sw'][1]
    Q_sw_i = Q[3] * Y_ij['sw'][2]

    # within sw: 1) advection in (p, c, i), 2) advection out (p, c, i)
    # assume waterflow rate is 90% of freshwater
    swOut_adv_c = adv.D_advec_water(climate['waterflow_d'][i] * 0.9, Z_ij['sw'][0]) * Q_sw_p
    swOut_adv_i = adv.D_advec_water(climate['waterflow_d'][i] * 0.9, Z_ij['sw'][2]) * Q_sw_i
    swIn_adv_c = fwOut_adv_c
    swIn_adv_i = fwOut_adv_i
    sSSOut_adv_p = adv.D_advec_susSed(climate['waterflow_d'][i] * 0.9, env['sSSVf'],
                                      Z_ij['sw'][0]) * Q_sw_p
    sSSIn_adv_p = fSSOut_adv_p

    # sw to other subcomparts: 1) diffusion to sSedW (c, i)
    sw_diff_c_to_sSedW = diff.D_diffu_comp1_comp2(7.0*10**(-7), env['swA']) * Q_sw_c
    sw_diff_i_to_sSedW = diff.D_diffu_comp1_comp2(1.0*10**(-5), env['swA']) * Q_sw_i

    # sSS: 1) deposition (p)
    k_dep_dry_sSS = adv.k_dep_dry(env['seassP'], env['seawP'],
                                  env['dynViscSW'], env['radiusParticlesSW'])
    sSS_dep_p = adv.D_dep_dry(k_dep_dry_sSS, env['seawA'] * env['sSSVf'], Z_ij['sw'][0]) * Q_sw_p

    sSS_resusp_p = adv.D_aeroResusp(climate['windspeed_s'][i], env['coastalA'],
                                    chemParams['enrichFactor'], env['seawD'], Z_ij['sw'][0]) * Q_sw_p


    # seawater sediment - swSed, Q[4]
    Q_swSed_p = Q[4] * Y_ij['swSed'][0]
    Q_swSed_c = Q[4] * Y_ij['swSed'][1]
    Q_swSed_i = Q[4] * Y_ij['swSed'][2]

    # sSedW: diffusion (n, i) to sw
    sSedW_diff_c_to_sw = diff.D_diffu_comp1_comp2(7.0*10**(-7), env['swA']) * Q_swSed_c
    sSedW_diff_i_to_sw = diff.D_diffu_comp1_comp2(1.0*10**(-5), env['swA']) * Q_swSed_i

    # sSedS: burial (p)
    sSedS_burial_p = adv.D_burial(env['burialRateSW'], env['swA'], Z_ij['swSed'][0]) * Q_swSed_p

    # sSedS to sSS: resuspension (p)
    sSedS_resusp_p = adv.D_sedResusp(env['resuspensionRateSW'], env['swA'], Z_ij['swSed'][0]) * Q_swSed_p

    # sSed: advection out
    sSed_adv_outflow_p = adv.D_advec_water(climate['waterflow_d'][i] * env['swadvfrac'], Z_ij['swSed'][0]) * Q_swSed_p
    sSed_adv_outflow_c = adv.D_advec_water(climate['waterflow_d'][i] * env['swadvfrac'], Z_ij['swSed'][1]) * Q_swSed_c
    sSed_adv_outflow_i = adv.D_advec_water(climate['waterflow_d'][i] * env['swadvfrac'], Z_ij['swSed'][2]) * Q_swSed_i


    # soil 1 - soil1, Q[5]
    Q_soil1_p = Q[5] * Y_ij['soil1'][0]
    Q_soil1_c = Q[5] * Y_ij['soil1'][1]
    Q_soil1_i = Q[5] * Y_ij['soil1'][2]

    # soil1 to other compartments: 1) water runoff (c, i) to fw, 2) infiltration (c, i) to deep soil deepS1,
    # 3) soil erosion to fSS (p), 4) wind erosion (p)
    soilW1_runoff_c = adv.D_runoff(climate['precip_mm'][i], env['CN1'], env['soilA1'], Z_ij['soil1'][1]) * Q_soil1_c
    soilW1_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN1'], env['soilA1'], Z_ij['soil1'][2]) * Q_soil1_i
    # precip_mm, CN, evap_mm, FC, soilWC, soilV, soilA, Z_water)
    D_infil_c_1, k_infil_c_1 = adv.D_infiltra(climate['precip_mm'][i], env['CN1'], climate['evap_mm'][i], env['FC1'],
                                              env['soilWC1'], env['soilV1'], env['soilA1'], Z_ij['fw'][1])

    D_infil_i_1, k_infil_i_1 = adv.D_infiltra(climate['precip_mm'][i], env['CN1'], climate['evap_mm'][i], env['FC1'],
                                              env['soilWC1'], env['soilV1'], env['soilA1'], Z_ij['fw'][2])

    soilW1_infil_c = D_infil_c_1 * Q_soil1_c
    soilW1_infil_i = D_infil_i_1 * Q_soil1_i


    soilS1_erosion_p = adv.D_erosion(climate['precip_mm'][i], env['slope1'],
                                     env['Kfact1'], env['cropManageFactor1'],
                                     env['supportFactor1'], env['soilA1'],
                                     env['soilP1'], Z_ij['soil1'][0]) * Q_soil1_p

    soilS1_windErosion_p = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness1'],
                                             env['Kconstant1'], env['airP'], env['soilA1'], env['A1'],
                                             env['TSV1'], env['TSVmin1'], env['z_wind1'], env['percWind1'],
                                             env['windConstant1'], env['percUncovered1'], env['percSuspended1'],
                                             env['soilP1'], Z_ij['soil1'][0]) * Q_soil1_p

    # deep soil 1 - deepS1, Q[6]
    Q_deepS1_p = Q[6] * Y_ij['deepS1'][0]
    Q_deepS1_c = Q[6] * Y_ij['deepS1'][1]
    Q_deepS1_i = Q[6] * Y_ij['deepS1'][2]

    # deepS1 to other compartments: 1) leaching (c, i) to fw
    deepS1_leach_c = adv.D_leach(k_infil_c_1, Z_ij['fw'][1]) * Q_deepS1_c
    deepS1_leach_i = adv.D_leach(k_infil_i_1, Z_ij['fw'][2]) * Q_deepS1_i


    # soil 2 - soil2, Q[7]
    Q_soil2_p = Q[7] * Y_ij['soil2'][0]
    Q_soil2_c = Q[7] * Y_ij['soil2'][1]
    Q_soil2_i = Q[7] * Y_ij['soil2'][2]

    # soilW2 to other subcomparts: 1) water runoff (c, i) to fw, 2) infiltration (c, i) to deep soil deepS2,
    # 3) erosion (p) to fSS, 4) wind erosion (p)
    soilW2_runoff_c = adv.D_runoff(climate['precip_mm'][i], env['CN2'], env['soilA2'], Z_ij['soil2'][1]) * Q_soil2_c
    soilW2_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN2'], env['soilA2'], Z_ij['soil2'][2]) * Q_soil2_i

    D_infil_c_2, k_infil_c_2 = adv.D_infiltra(climate['precip_mm'][i], env['CN2'], climate['evap_mm'][i], env['FC2'],
                                              env['soilWC2'], env['soilV2'], env['soilA2'], Z_ij['fw'][1])

    D_infil_i_2, k_infil_i_2 = adv.D_infiltra(climate['precip_mm'][i], env['CN2'], climate['evap_mm'][i], env['FC2'],
                                              env['soilWC2'], env['soilV2'], env['soilA2'], Z_ij['fw'][2])

    soilW2_infil_c = D_infil_c_2 * Q_soil2_c
    soilW2_infil_i = D_infil_i_2 * Q_soil2_i

    soilS2_erosion_p = adv.D_erosion(climate['precip_mm'][i], env['slope2'],
                                     env['Kfact2'], env['cropManageFactor2'],
                                     env['supportFactor2'], env['soilA2'],
                                     env['soilP2'], Z_ij['soil2'][0]) * Q_soil2_p
    soilS2_windErosion_p = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness2'],
                                             env['Kconstant2'], env['airP'], env['soilA2'], env['A2'],
                                             env['TSV2'], env['TSVmin2'], env['z_wind2'], env['percWind2'],
                                             env['windConstant2'], env['percUncovered2'], env['percSuspended2'],
                                             env['soilP2'], Z_ij['soil2'][0]) * Q_soil2_p


    # deep soil 2 - deepS2, Q[8]
    Q_deepS2_p = Q[8] * Y_ij['deepS2'][0]
    Q_deepS2_c = Q[8] * Y_ij['deepS2'][1]
    Q_deepS2_i = Q[8] * Y_ij['deepS2'][2]

    # deepS2 to other compartments: 1) leaching (c, i) to fw
    deepS2_leach_c = adv.D_leach(k_infil_c_2, Z_ij['fw'][1]) * Q_deepS2_c
    deepS2_leach_i = adv.D_leach(k_infil_i_2, Z_ij['fw'][2]) * Q_deepS2_i

    # soil 3 - soil3, Q[9]
    Q_soil3_p = Q[9] * Y_ij['soil3'][0]
    Q_soil3_c = Q[9] * Y_ij['soil3'][1]
    Q_soil3_i = Q[9] * Y_ij['soil3'][2]

    # soilW3 to other subcomparts: 1) water runoff (c, i) to fw, 2) infiltration (c, i) to deep soil deepS3,
    # 3) erosion (p) to fSS, 4) wind erosion (p)
    soilW3_runoff_c = adv.D_runoff(climate['precip_mm'][i], env['CN3'], env['soilA3'], Z_ij['soil3'][1]) * Q_soil3_c
    soilW3_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN3'], env['soilA3'], Z_ij['soil3'][2]) * Q_soil3_i

    D_infil_c_3, k_infil_c_3 = adv.D_infiltra(climate['precip_mm'][i], env['CN3'], climate['evap_mm'][i], env['FC3'],
                                              env['soilWC3'], env['soilV3'], env['soilA3'], Z_ij['fw'][1])

    D_infil_i_3, k_infil_i_3 = adv.D_infiltra(climate['precip_mm'][i], env['CN3'], climate['evap_mm'][i], env['FC3'],
                                              env['soilWC3'], env['soilV3'], env['soilA3'], Z_ij['fw'][2])

    soilW3_infil_c = D_infil_c_3 * Q_soil3_c
    soilW3_infil_i = D_infil_i_3 * Q_soil3_i

    soilS3_erosion_p = adv.D_erosion(climate['precip_mm'][i], env['slope3'],
                                     env['Kfact3'], env['cropManageFactor3'],
                                     env['supportFactor3'], env['soilA3'],
                                     env['soilP3'], Z_ij['soil3'][0]) * Q_soil3_p
    soilS3_windErosion_p = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness3'],
                                             env['Kconstant3'], env['airP'], env['soilA3'],
                                             env['A3'],
                                             env['TSV3'], env['TSVmin3'], env['z_wind3'],
                                             env['percWind3'],
                                             env['windConstant3'], env['percUncovered3'],
                                             env['percSuspended3'],
                                             env['soilP3'], Z_ij['soil3'][0]) * Q_soil3_p


    # deep soil 3 - deepS3, Q[10]
    Q_deepS3_p = Q[10] * Y_ij['deepS3'][0]
    Q_deepS3_c = Q[10] * Y_ij['deepS3'][1]
    Q_deepS3_i = Q[10] * Y_ij['deepS3'][2]

    # deepS3 to other compartment: 1) leaching (c, i) to fw
    deepS3_leach_c = adv.D_leach(k_infil_c_3, Z_ij['fw'][1]) * Q_deepS3_c
    deepS3_leach_i = adv.D_leach(k_infil_i_3, Z_ij['fw'][2]) * Q_deepS3_i

    # soil 4 - soil4, Q[11]
    Q_soil4_p = Q[11] * Y_ij['soil4'][0]
    Q_soil4_c = Q[11] * Y_ij['soil4'][1]
    Q_soil4_i = Q[11] * Y_ij['soil4'][2]

    # soilW4 to other subcomparts: 1) water runoff (c, i) to fw, 2) infiltration (c, i) to deep soil deepS4,
    # 3) erosion (p) to fSS, 4) wind erosion (p)
    soilW4_runoff_c = adv.D_runoff(climate['precip_mm'][i], env['CN4'], env['soilA4'], Z_ij['soil4'][1]) * Q_soil4_c
    soilW4_runoff_i = adv.D_runoff(climate['precip_mm'][i], env['CN4'], env['soilA4'], Z_ij['soil4'][2]) * Q_soil4_i

    D_infil_c_4, k_infil_c_4 = adv.D_infiltra(climate['precip_mm'][i], env['CN4'], climate['evap_mm'][i], env['FC4'],
                                              env['soilWC4'], env['soilV4'], env['soilA4'], Z_ij['fw'][1])

    D_infil_i_4, k_infil_i_4 = adv.D_infiltra(climate['precip_mm'][i], env['CN4'], climate['evap_mm'][i], env['FC4'],
                                              env['soilWC4'], env['soilV4'], env['soilA4'], Z_ij['fw'][2])

    soilW4_infil_c = D_infil_c_4 * Q_soil4_c
    soilW4_infil_i = D_infil_i_4 * Q_soil4_i

    soilS4_erosion_p = adv.D_erosion(climate['precip_mm'][i], env['slope4'],
                                     env['Kfact4'], env['cropManageFactor4'],
                                     env['supportFactor4'], env['soilA4'],
                                     env['soilP4'], Z_ij['soil4'][0]) * Q_soil4_p
    soilS4_windErosion_p = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i],
                                             env['roughness4'],
                                             env['Kconstant4'], env['airP'], env['soilA4'],
                                             env['A4'],
                                             env['TSV4'], env['TSVmin4'], env['z_wind4'],
                                             env['percWind4'],
                                             env['windConstant4'], env['percUncovered4'],
                                             env['percSuspended4'],
                                             env['soilP4'], Z_ij['soil4'][0]) * Q_soil4_p

    # deep soil 4 - deepS4, Q[12]
    Q_deepS4_p = Q[12] * Y_ij['deepS4'][0]
    Q_deepS4_c = Q[12] * Y_ij['deepS4'][1]
    Q_deepS4_i = Q[12] * Y_ij['deepS4'][2]

    # deepS4 to other compartment: 1) leaching (c, i) to fw
    deepS4_leach_c = adv.D_leach(k_infil_c_4, Z_ij['fw'][1]) * Q_deepS4_c
    deepS4_leach_i = adv.D_leach(k_infil_i_4, Z_ij['fw'][2]) * Q_deepS4_i

    ###################################################################
    # processes output to transport rate kg/day
    # N * molar mass = mol/day * kg/mol = kg/day
    ###################################################################

    # 1) advection process
    adv_air_in = (airIn_adv_p + airIn_adv_i) * chemParams['molar_mass']
    adv_air_out = (aerOut_adv_p + aerOut_adv_c + aerOut_adv_i) * chemParams['molar_mass']
    adv_fw_in = (fwIn_adv_p + fwIn_adv_c + fwIn_adv_i) * chemParams['molar_mass']
    adv_fw_out = (fwOut_adv_c + fwOut_adv_i + fSSOut_adv_p) * chemParams['molar_mass']
    adv_fwSed_in = (fSed_adv_inflow_p + fSed_adv_inflow_c + fSed_adv_inflow_i) * chemParams['molar_mass']
    adv_fwSed_out = (fSed_adv_outflow_p + fSed_adv_outflow_c + fSed_adv_outflow_i) * chemParams['molar_mass']
    adv_sw_in = (swIn_adv_c + swIn_adv_i + sSSIn_adv_p) * chemParams['molar_mass']
    adv_sw_out = (swOut_adv_c + swOut_adv_i + sSSOut_adv_p) * chemParams['molar_mass']
    adv_swSed_out = (sSed_adv_outflow_p + sSed_adv_outflow_c + sSed_adv_outflow_i) * chemParams['molar_mass']

    # 2) deposition process
    dep_dry_air = (aer_dep_dry_p + aer_dep_dry_i) * chemParams['molar_mass']
    dep_dry_air_fw = (aer_dep_dry_p_to_fSS + aer_dep_dry_i_to_fSS) * chemParams['molar_mass']
    dep_dry_air_sw = (aer_dep_dry_p_to_sSS + aer_dep_dry_i_to_sSS) * chemParams['molar_mass']
    dep_dry_air_soil1 = (aer_dep_dry_p_to_soil1 + aer_dep_dry_i_to_soil1) * chemParams['molar_mass']
    dep_dry_air_soil2 = (aer_dep_dry_p_to_soil2 + aer_dep_dry_i_to_soil2) * chemParams['molar_mass']
    dep_dry_air_soil3 = (aer_dep_dry_p_to_soil3 + aer_dep_dry_i_to_soil3) * chemParams['molar_mass']
    dep_dry_air_soil4 = (aer_dep_dry_p_to_soil4 + aer_dep_dry_i_to_soil4) * chemParams['molar_mass']

    dep_wet_air = (aer_dep_wet_p + aer_dep_wet_i) * chemParams['molar_mass']
    dep_wet_air_fw = (aer_dep_wet_p_to_fSS + aer_dep_wet_i_to_fSS) * chemParams['molar_mass']
    dep_wet_air_sw = (aer_dep_wet_p_to_sSS + aer_dep_wet_i_to_sSS) * chemParams['molar_mass']
    dep_wet_air_soil1 = (aer_dep_wet_p_to_soil1 + aer_dep_wet_i_to_soil1) * chemParams['molar_mass']
    dep_wet_air_soil2 = (aer_dep_wet_p_to_soil2 + aer_dep_wet_i_to_soil2) * chemParams['molar_mass']
    dep_wet_air_soil3 = (aer_dep_wet_p_to_soil3 + aer_dep_wet_i_to_soil3) * chemParams['molar_mass']
    dep_wet_air_soil4 = (aer_dep_wet_p_to_soil4 + aer_dep_wet_i_to_soil4) * chemParams['molar_mass']

    dep_fSS = fSS_dep_p * chemParams['molar_mass']
    dep_sSS = sSS_dep_p * chemParams['molar_mass']

    # 3) diffusion process
    diff_fw_fSedW = (fw_diff_c_to_fSedW + fw_diff_i_to_fSedW) * chemParams['molar_mass']
    diff_sw_sSedW = (sw_diff_c_to_sSedW + sw_diff_i_to_sSedW) * chemParams['molar_mass']
    diff_fSedW_fw = (fSedW_diff_c_to_fw + fSedW_diff_i_to_fw) * chemParams['molar_mass']
    diff_sSedW_sw = (sSedW_diff_c_to_sw + sSedW_diff_i_to_sw) * chemParams['molar_mass']

    # 4) other process
    burial_fwSed = fSedS_burial_p * chemParams['molar_mass']
    burial_swSed = sSedS_burial_p * chemParams['molar_mass']
    resusp_fwSed = fSedS_resusp_p * chemParams['molar_mass']
    resusp_swSed = sSedS_resusp_p * chemParams['molar_mass']
    aero_resusp_sSS = sSS_resusp_p * chemParams['molar_mass']
    runoff_soil1 = (soilW1_runoff_c + soilW1_runoff_i) * chemParams['molar_mass']
    runoff_soil2 = (soilW2_runoff_c + soilW2_runoff_i) * chemParams['molar_mass']
    runoff_soil3 = (soilW3_runoff_c + soilW3_runoff_i) * chemParams['molar_mass']
    runoff_soil4 = (soilW4_runoff_c + soilW4_runoff_i) * chemParams['molar_mass']
    erosion_soil1 = soilS1_erosion_p * chemParams['molar_mass']
    erosion_soil2 = soilS2_erosion_p * chemParams['molar_mass']
    erosion_soil3 = soilS3_erosion_p * chemParams['molar_mass']
    erosion_soil4 = soilS4_erosion_p * chemParams['molar_mass']
    wind_erosion_soil1 = soilS1_windErosion_p * chemParams['molar_mass']
    wind_erosion_soil2 = soilS2_windErosion_p * chemParams['molar_mass']
    wind_erosion_soil3 = soilS3_windErosion_p * chemParams['molar_mass']
    wind_erosion_soil4 = soilS4_windErosion_p * chemParams['molar_mass']
    infiltra_soil1 = (soilW1_infil_c + soilW1_infil_i) * chemParams['molar_mass']
    infiltra_soil2 = (soilW2_infil_c + soilW2_infil_i) * chemParams['molar_mass']
    infiltra_soil3 = (soilW3_infil_c + soilW3_infil_i) * chemParams['molar_mass']
    infiltra_soil4 = (soilW4_infil_c + soilW4_infil_i) * chemParams['molar_mass']
    leach_soil1 = (deepS1_leach_c + deepS1_leach_i) * chemParams['molar_mass']
    leach_soil2 = (deepS2_leach_c + deepS2_leach_i) * chemParams['molar_mass']
    leach_soil3 = (deepS3_leach_c + deepS3_leach_i) * chemParams['molar_mass']
    leach_soil4 = (deepS4_leach_c + deepS4_leach_i) * chemParams['molar_mass']

    processes = [adv_air_in, adv_air_out, adv_fw_in, adv_fw_out, adv_fwSed_in, adv_fwSed_out, adv_sw_in, adv_sw_out,
                 adv_swSed_out, dep_dry_air, dep_dry_air_fw, dep_dry_air_sw, dep_dry_air_soil1, dep_dry_air_soil2,
                 dep_dry_air_soil3, dep_dry_air_soil4, dep_wet_air, dep_wet_air_fw, dep_wet_air_sw, dep_wet_air_soil1,
                 dep_wet_air_soil2, dep_wet_air_soil3, dep_wet_air_soil4, dep_fSS, dep_sSS, diff_fw_fSedW, diff_sw_sSedW,
                 diff_fSedW_fw, diff_sSedW_sw, burial_fwSed, burial_swSed, resusp_fwSed, resusp_swSed, aero_resusp_sSS,
                 runoff_soil1, runoff_soil2, runoff_soil3, runoff_soil4, erosion_soil1, erosion_soil2, erosion_soil3,
                 erosion_soil4, wind_erosion_soil1, wind_erosion_soil2, wind_erosion_soil3, wind_erosion_soil4,
                 infiltra_soil1, infiltra_soil2, infiltra_soil3, infiltra_soil4, leach_soil1, leach_soil2, leach_soil3,
                 leach_soil4]

    return processes

