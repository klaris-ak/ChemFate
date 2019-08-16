from __future__ import division
from degradation_process import Degradation
from advective_processes import AdvectiveProcess
from diffusion_process_non_ion import Diffusion, MTC
from Z_non_ion import zValue


def org_process(f, i, env, climate, chemParams, bgConc):
    # differential equation solver for organic chemical fugacity in all compartments
    # t is time, f is the fugacity by compartment and day (y for equations), i is the iteration in
    # the for loop - so the time step, V is the volume vector
    # flux N = D*f, unit: mol/day = mol/(Pa-day) * Pa
    # changes in fugacity in each compartment as a function of time
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


    MTC_process = MTC(chemParams['molar_volume'], climate['temp_K'][i], chemParams['MW'])
    airMD = MTC_process.airMD()
    waterMD = MTC_process.waterMD()
    airWaterMTC = MTC_process.airWaterMTC(airMD)
    airSoilMTC = MTC_process.airSoilMTC(airMD)
    soil1AirMTC = MTC_process.soilAirMTC(airMD, env['soilAC1'], env['soilWC1'])
    soil2AirMTC = MTC_process.soilAirMTC(airMD, env['soilAC2'], env['soilWC2'])
    soil3AirMTC = MTC_process.soilAirMTC(airMD, env['soilAC3'], env['soilWC3'])
    soil4AirMTC = MTC_process.soilAirMTC(airMD, env['soilAC4'], env['soilWC4'])
    soil1WaterMTC = MTC_process.soilWaterMTC(waterMD, env['soilAC1'], env['soilWC1'])
    soil2WaterMTC = MTC_process.soilWaterMTC(waterMD, env['soilAC2'], env['soilWC2'])
    soil3WaterMTC = MTC_process.soilWaterMTC(waterMD, env['soilAC3'], env['soilWC3'])
    soil4WaterMTC = MTC_process.soilWaterMTC(waterMD, env['soilAC4'], env['soilWC4'])
    waterAirMTC = MTC_process.waterAirMTC(waterMD)
    sedmtWaterMTC = MTC_process.sedmtWaterMTC(waterMD, env['fsedpercSolid'])

    deg = Degradation()
    adv = AdvectiveProcess()
    diff = Diffusion()

    # air f[0]
    # air degradation
    dDegAirBulk = deg.D_deg(env['airV'], chemParams['kDeg_air_n'], zAirSub) + deg.D_deg(env['aerV'],
                                                                                        chemParams['kDeg_aer_n'],
                                                                                        zAerSub)
    degradationAir = dDegAirBulk * f[0]  # mol/day

    # air advection (outflow)
    QAdvecAir = adv.G_advec_air(climate['windspeed_d'][i], env['airA'], env['airH'])
    advectionAirOutflow = adv.D_advec_air(climate['windspeed_d'][i], env['airA'], env['airH'], zAirSub) * f[0]
    advectionAerOutflow = adv.D_advec_aer(climate['windspeed_d'][i], env['airA'], env['airH'], env['aerVf'], zAerSub) * \
                          f[0]
    advectionAirBulkOutflow = advectionAirOutflow + advectionAerOutflow
    # air advection (inflow) - air from outside of the compartment
    advectionAirInflow = QAdvecAir * bgConc['gairc_n']  # mol/day

    # aerosol rain dissolution to water
    rainDissolutionAer = adv.D_rain_diss(climate['precip_m'][i], env['area'], zWaterSub) * f[0]
    rainDissolutionAer2FW = rainDissolutionAer * (env['freshwA'] / env['area'])
    rainDissolutionAer2SW = rainDissolutionAer * (env['seawA'] / env['area'])
    rainDissolutionAer2Soil1 = rainDissolutionAer * (env['soilA1'] / env['area'])
    rainDissolutionAer2Soil2 = rainDissolutionAer * (env['soilA2'] / env['area'])
    rainDissolutionAer2Soil3 = rainDissolutionAer * (env['soilA3'] / env['area'])
    rainDissolutionAer2Soil4 = rainDissolutionAer * (env['soilA4'] / env['area'])

    # aerosol dry depostion
    kAerDryDep = adv.k_dep_dry(env['aerP'], env['airP'], env['dynViscAir'], env['radiusParticlesAer'])
    dryDepositionAer = adv.D_dep_dry(kAerDryDep, env['area'], zAerSub) * f[0]
    dryDepositionAer2FW = dryDepositionAer * (env['freshwA'] / env['area'])
    dryDepositionAer2SW = dryDepositionAer * (env['seawA'] / env['area'])
    dryDepositionAer2Soil1 = dryDepositionAer * (env['soilA1'] / env['area'])
    dryDepositionAer2Soil2 = dryDepositionAer * (env['soilA2'] / env['area'])
    dryDepositionAer2Soil3 = dryDepositionAer * (env['soilA3'] / env['area'])
    dryDepositionAer2Soil4 = dryDepositionAer * (env['soilA4'] / env['area'])

    # aerosol wet deposition
    wetDepositionAer = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'], env['area'], zAerSub) * env['aerVf'] * f[0]
    wetDepositionAer2FW = wetDepositionAer * (env['freshwA'] / env['area'])
    wetDepositionAer2SW = wetDepositionAer * (env['seawA'] / env['area'])
    wetDepositionAer2Soil1 = wetDepositionAer * (env['soilA1'] / env['area'])
    wetDepositionAer2Soil2 = wetDepositionAer * (env['soilA2'] / env['area'])
    wetDepositionAer2Soil3 = wetDepositionAer * (env['soilA3'] / env['area'])
    wetDepositionAer2Soil4 = wetDepositionAer * (env['soilA4'] / env['area'])

    # air to water diffusion
    diffusionAir2FW = diff.diffusion_air_water(airWaterMTC, waterAirMTC, env['freshwA'], zAirSub, zWaterSub) * f[0]
    diffusionAir2SW = diff.diffusion_air_water(airWaterMTC, waterAirMTC, env['seawA'], zAirSub, zWaterSub) * f[0]

    # air to soil diffusion
    diffusionAir2S1surf = diff.diffusion_air_soil(airSoilMTC, soil1AirMTC, soil1WaterMTC, env['soilA1'], zAirSub,
                                                  zWaterSub) * f[0]
    diffusionAir2S2surf = diff.diffusion_air_soil(airSoilMTC, soil1AirMTC, soil1WaterMTC, env['soilA2'], zAirSub,
                                                  zWaterSub) * f[0]
    diffusionAir2S3surf = diff.diffusion_air_soil(airSoilMTC, soil1AirMTC, soil1WaterMTC, env['soilA3'], zAirSub,
                                                  zWaterSub) * f[0]
    diffusionAir2S4surf = diff.diffusion_air_soil(airSoilMTC, soil1AirMTC, soil1WaterMTC, env['soilA4'], zAirSub,
                                                  zWaterSub) * f[0]

    # freshwater f[1]
    # freshwater degradation
    dDegFWBulk = deg.D_deg(env['fWaterV'], chemParams['kDeg_fw_n'], zWaterSub) + \
                 deg.D_deg(env['fSSV'], chemParams['kDeg_fSS_n'], zFWSusSedSub)
    degradationFW = dDegFWBulk * f[1]

    # freshwater advection Inflow/outflow
    advectionFWInflow = adv.D_advec_water(climate['waterflow_d'][i], bgConc['gfreshwc_n'])  # mol/day
    advectionFWBulkInflow = advectionFWInflow
    advectionFWOutflow = adv.D_advec_water(climate['waterflow_d'][i], zWaterSub) * f[1]
    advectionFWSusSedOutflow = adv.D_advec_susSed(climate['waterflow_d'][i], env['fSSVf'], zFWSusSedSub) * f[1]
    advectionFWBulkOutflow = advectionFWOutflow + advectionFWSusSedOutflow

    # freshwater to air diffusion
    diffusionFW2Air = diff.diffusion_air_water(airWaterMTC, waterAirMTC, env['freshwA'], zAirSub, zWaterSub) * f[1]
    # freshwater to sediment diffusion
    diffusionFW2FWSed = diff.diffusion_sediment_water(sedmtWaterMTC, env['freshwA'], zWaterSub) * f[1]

    # freshwater sus sed deposition
    kFWSedDep = adv.k_dep_dry(env['freshssP'], env['freshwP'], env['dynViscFW'], env['radiusParticlesFW'])
    sedDepFWSusSed = adv.D_dep_dry(kFWSedDep, env['freshwA'], zFWSusSedSub) * f[1]


    # freshwater sediment f[2]
    # freshwater sediment degradation
    dDegFWSed = deg.D_deg(env['fSedSV'], chemParams['kDeg_fSedS_n'], zFSedSSub) + \
                deg.D_deg(env['fSedWV'], chemParams['kDeg_fSedW_n'], zWaterSub)

    degradationFSed = dDegFWSed * f[2]

    # freshwater sediment advection
    advectionFWSedInflow = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_n'])
    advectionFWSedOutflow = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], zFWSedimentBulk) * f[2]
    # freshwater sediment diffusion
    diffusionFWSed2FW = diff.diffusion_sediment_water(sedmtWaterMTC, env['freshwA'], zWaterSub) * f[2]
    # freshwater sediment solid to freshwater suspended sediment - resuspension
    resuspFWSed = adv.D_sedResusp(env['freshwA'], env['resuspensionRateFW'], zFSedSSub) * f[2]

    # freshwater sediment solid burial
    sedBurialFSed = adv.D_burial(env['freshwA'], env['burialRateFW'], zFSedSSub) * f[2]


    # seawater f[3]
    # seawater degradation
    dDegSWBulk = deg.D_deg(env['sWaterV'], chemParams['kDeg_sw_n'], zWaterSub) + \
                 deg.D_deg(env['sSSV'], chemParams['kDeg_sSS_n'], zSWSusSedSub)
    degradationSW = dDegSWBulk * f[3]

    advectionSWBulkInflow = advectionFWBulkOutflow
    # seawater advection outflow
    advectionSWOutflow = adv.D_advec_water(10 * climate['waterflow_d'][i], zWaterSub) * f[3]
    advectionSWSusSedOutflow = adv.D_advec_susSed(10 * climate['waterflow_d'][i], env['sSSVf'], zSWSusSedSub) * f[3]
    advectionSWBulkOutflow = advectionSWOutflow + advectionSWSusSedOutflow
    # seawater to air diffusion
    diffusionSW2Air = diff.diffusion_air_water(airWaterMTC, waterAirMTC, env['seawA'], zAirSub, zWaterSub) * f[3]
    diffusionSW2SWSed = diff.diffusion_sediment_water(sedmtWaterMTC, env['seawA'], zWaterSub) * f[3]
    # seawater sus sed deposition
    kSWSedDep = adv.k_dep_dry(env['seassP'], env['seawP'], env['dynViscSW'], env['radiusParticlesSW'])
    sedDepSWSusSed = adv.D_dep_dry(kSWSedDep, env['seawA'] * env['sSSVf'], zSWSusSedSub) * f[3]

    # seawater sediment f[4]
    # seawater sediment degradation
    dDegSWSed = deg.D_deg(env['sSedWV'], chemParams['kDeg_sSedW_n'], zWaterSub) + \
                deg.D_deg(env['sSedSV'], chemParams['kDeg_sSedS_n'], zSSedSSub)
    degradationSSed = dDegSWSed * f[4]
    # seawater sediment advection
    advectionSWSedOutflow = adv.D_advec_water(climate['waterflow_d'][i] * env['swadvfrac'], 10 * zSWSedimentBulk) * f[4]
    diffusionSWSed2SW = diff.diffusion_sediment_water(sedmtWaterMTC, env['seawA'], zWaterSub) * f[4]
    # seawater sediment resuspension
    resuspSWSed = adv.D_sedResusp(env['seawA'], env['resuspensionRateSW'], zSSedSSub) * f[4]
    # seawater sediment solid burial
    sedBurialSSed = adv.D_burial(env['seawA'], env['burialRateSW'], zSSedSSub) * f[4]

    # soil 1 surface soil f[5]
    # soil diffusion to air
    diffusionS1surf2Air = diff.diffusion_air_soil(airSoilMTC, soil1AirMTC, soil1WaterMTC, env['soilA1'], zAirSub,
                                                  zWaterSub) * f[5]
    # soil degradation (deg soil air + deg soil water + deg soil solid)
    dDegSoil1Bulk = deg.D_deg(env['soilAV1'], chemParams['kDeg_soilA1_n'], zAirSub) + \
                    deg.D_deg(env['soilWV1'], chemParams['kDeg_soilW1_n'], zWaterSub) + \
                    deg.D_deg(env['soilSV1'], chemParams['kDeg_soilS1_n'], zS1SolidSub)

    degradationSoil1surf = dDegSoil1Bulk * f[5]

    # soil 1 water to freshwater - runoff
    runoffS1Water = adv.D_runoff(climate['precip_mm'][i], env['CN1'], env['soilA1'], zWaterSub) * f[5]

    # soil 1 infiltration from surface soil to deep soil
    D_infil_1, k_infil_1 = adv.D_infiltra(climate['precip_mm'][i], env['CN1'], climate['evap_mm'][i], env['FC1'], env['soilWC1'],
                                          env['soilV1'], env['soilA1'], zWaterSub)
    infiltrationSoil1surf = D_infil_1 * f[5]

    # soil 1 solid to freshwater sus sed - erosion
    erosionS1Solid = adv.D_erosion(climate['precip_mm'][i], env['slope1'], env['Kfact1'], env['cropManageFactor1'],
                                   env['supportFactor1'], env['soilA1'], env['soilP1'], zS1SolidSub) * f[5]

    # soil 1 solid to air - wind erosion
    windErosionS1Solid = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness1'],
                                           env['Kconstant1'], env['airP'], env['soilA1'], env['A1'], env['TSV1'],
                                           env['TSVmin1'], env['z_wind1'], env['percWind1'],
                                           env['windConstant1'], env['percUncovered1'], env['percSuspended1'],
                                           env['soilP1'], zS1SolidSub) * f[5]

    # soil 1 deep soil f[6]
    # soil 1 deep soil to freshwater body
    leachingS1deep = adv.D_leach(k_infil_1, zWaterSub) * f[6]
    degradationSoil1deep = deg.D_deg(env['deepSV1'], chemParams['kDeg_deepS1_n'], zS1DeepSSub) * f[6]

    # soil 2 surface soil f[7]
    # soil diffusion to air
    diffusionS2surf2Air = diff.diffusion_air_soil(airSoilMTC, soil2AirMTC, soil2WaterMTC, env['soilA2'], zAirSub,
                                                  zWaterSub) * f[7]
    # soil degradation (deg soil air + deg soil water + deg soil solid)
    degradationSoil2surf = deg.D_deg(env['soilV2'], chemParams['kDeg_soilS2_n'], zSoil2Bulk) * f[7]
    # soil 2 water to freshwater - runoff
    runoffS2Water = adv.D_runoff(climate['precip_mm'][i], env['CN2'], env['soilA2'], zWaterSub) * f[7]

    # soil 2 infiltration from surface soil to deep soil
    D_infil_2, k_infil_2 = adv.D_infiltra(climate['precip_mm'][i], env['CN2'], climate['evap_mm'][i], env['FC2'], env['soilWC2'],
                                          env['soilV2'], env['soilA2'], zWaterSub)
    infiltrationSoil2surf = D_infil_2 * f[7]

    # soil 2 solid to freshwater sus sed - erosion
    erosionS2Solid = adv.D_erosion(climate['precip_mm'][i], env['slope2'], env['Kfact2'], env['cropManageFactor2'],
                                   env['supportFactor2'], env['soilA2'], env['soilP2'], zS2SolidSub) * f[7]

    # soil 2 solid to air - wind erosion
    windErosionS2Solid = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness2'],
                                           env['Kconstant2'], env['airP'], env['soilA2'], env['A2'], env['TSV2'],
                                           env['TSVmin2'], env['z_wind2'], env['percWind2'],
                                           env['windConstant2'], env['percUncovered2'], env['percSuspended2'],
                                           env['soilP2'], zS2SolidSub) * f[7]

    # soil 2 deep soil f[8]
    leachingS2deep = adv.D_leach(k_infil_2, zWaterSub) * f[8]
    degradationSoil2deep = deg.D_deg(env['deepSV2'], chemParams['kDeg_deepS2_n'], zS2DeepSSub) * f[8]

    # soil 3 surface soil f[9]
    # soil diffusion to air
    diffusionS3surf2Air = diff.diffusion_air_soil(airSoilMTC, soil3AirMTC, soil3WaterMTC, env['soilA3'], zAirSub,
                                                  zWaterSub) * f[9]
    # soil degradation (deg soil air + deg soil water + deg soil solid)
    degradationSoil3surf = deg.D_deg(env['soilV3'], chemParams['kDeg_soilS3_n'], zSoil3Bulk) * f[9]

    # soil water to freshwater - runoff
    runoffS3Water = adv.D_runoff(climate['precip_mm'][i], env['CN3'], env['soilA3'], zWaterSub) * f[9]

    # soil 3 infiltration
    D_infil_3, k_infil_3 = adv.D_infiltra(climate['precip_mm'][i], env['CN3'], climate['evap_mm'][i], env['FC3'], env['soilWC3'],
                                          env['soilV3'], env['soilA3'], zWaterSub)
    infiltrationSoil3surf = D_infil_3 * f[9]

    # soil 3 erosion
    erosionS3Solid = adv.D_erosion(climate['precip_mm'][i], env['slope3'], env['Kfact3'], env['cropManageFactor3'],
                                   env['supportFactor3'], env['soilA3'], env['soilP3'], zS3SolidSub) * f[9]

    # soil 3 solid to air - wind erosion
    windErosionS3Solid = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness3'],
                                           env['Kconstant3'], env['airP'], env['soilA3'], env['A3'], env['TSV3'],
                                           env['TSVmin3'], env['z_wind3'], env['percWind3'],
                                           env['windConstant3'], env['percUncovered3'], env['percSuspended3'],
                                           env['soilP3'], zS3SolidSub) * f[9]


    # soil 3 deep soil f[10]
    leachingS3deep = adv.D_leach(k_infil_3, zWaterSub) * f[10]
    degradationSoil3deep = deg.D_deg(env['deepSV3'], chemParams['kDeg_deepS3_n'], zS3DeepSSub) * f[10]

    # soil 4 surface soil f[11]
    # soil diffusion to air
    diffusionS4surf2Air = diff.diffusion_air_soil(airSoilMTC, soil4AirMTC, soil4WaterMTC, env['soilA4'], zAirSub,
                                                  zWaterSub) * f[11]

    # soil degradation (deg soil air + deg soil water + deg soil solid)
    degradationSoil4surf = deg.D_deg(env['soilV4'], chemParams['kDeg_soilS4_n'], zSoil4Bulk) * f[11]

    # soil water to freshwater - runoff
    runoffS4Water = adv.D_runoff(climate['precip_mm'][i], env['CN4'], env['soilA4'], zWaterSub) * f[11]

    # soil 4 infiltration
    D_infil_4, k_infil_4 = adv.D_infiltra(climate['precip_mm'][i], env['CN4'], climate['evap_mm'][i], env['FC4'], env['soilWC4'],
                                          env['soilV4'], env['soilA4'], zWaterSub)
    infiltrationSoil4surf = D_infil_4 * f[11]

    # soil 4 erosion
    erosionS4Solid = adv.D_erosion(climate['precip_mm'][i], env['slope4'], env['Kfact4'], env['cropManageFactor4'],
                                   env['supportFactor4'], env['soilA4'], env['soilP4'], zS4SolidSub) * f[11]

    # soil 4 solid to air - wind erosion
    windErosionS4Solid = adv.D_windErosion(climate['windspeed_s'][i], climate['precip_mm'][i], env['roughness4'],
                                           env['Kconstant4'], env['airP'], env['soilA4'], env['A4'], env['TSV4'],
                                           env['TSVmin4'], env['z_wind4'], env['percWind4'],
                                           env['windConstant4'], env['percUncovered4'], env['percSuspended4'],
                                           env['soilP4'], zS4SolidSub) * f[11]

    # soil 4 deep soil f[12]
    leachingS4deep = adv.D_leach(k_infil_4, zWaterSub) * f[12]
    degradationSoil4deep = deg.D_deg(env['deepSV4'], chemParams['kDeg_deepS4_n'], zS4DeepSSub) * f[12]

    ###################################################################
    # processes output to transport rate kg/day
    # N * molar mass = mol/day * kg/mol = kg/day
    ###################################################################

    # 1) degradation process
    deg_air = degradationAir * chemParams['molar_mass']
    deg_fw = degradationFW * chemParams['molar_mass']
    deg_fwSed = degradationFSed * chemParams['molar_mass']
    deg_sw = degradationSW * chemParams['molar_mass']
    deg_swSed = degradationSSed * chemParams['molar_mass']
    deg_soil1 = degradationSoil1surf * chemParams['molar_mass']
    deg_deepS1 = degradationSoil1deep * chemParams['molar_mass']
    deg_soil2 = degradationSoil2surf * chemParams['molar_mass']
    deg_deepS2 = degradationSoil2deep * chemParams['molar_mass']
    deg_soil3 = degradationSoil3surf * chemParams['molar_mass']
    deg_deepS3 = degradationSoil3deep * chemParams['molar_mass']
    deg_soil4 = degradationSoil4surf * chemParams['molar_mass']
    deg_deepS4 = degradationSoil4deep * chemParams['molar_mass']

    # 2) advection process
    adv_air_in = advectionAirInflow * chemParams['molar_mass']
    adv_air_out = advectionAirBulkOutflow * chemParams['molar_mass']
    adv_fw_in = advectionFWBulkInflow * chemParams['molar_mass']
    adv_fw_out = advectionFWBulkOutflow * chemParams['molar_mass']
    adv_fwSed_in = advectionFWSedInflow * chemParams['molar_mass']
    adv_fwSed_out = advectionFWSedOutflow * chemParams['molar_mass']
    adv_sw_in = advectionSWBulkInflow * chemParams['molar_mass']
    adv_sw_out = advectionSWBulkOutflow * chemParams['molar_mass']
    adv_swSed_out = advectionSWSedOutflow * chemParams['molar_mass']

    # 3) deposition process
    dep_dry_air = dryDepositionAer * chemParams['molar_mass']
    dep_dry_air_fw = dryDepositionAer2FW * chemParams['molar_mass']
    dep_dry_air_sw = dryDepositionAer2SW * chemParams['molar_mass']
    dep_dry_air_soil1 = dryDepositionAer2Soil1 * chemParams['molar_mass']
    dep_dry_air_soil2 = dryDepositionAer2Soil2 * chemParams['molar_mass']
    dep_dry_air_soil3 = dryDepositionAer2Soil3 * chemParams['molar_mass']
    dep_dry_air_soil4 = dryDepositionAer2Soil4 * chemParams['molar_mass']

    dep_wet_air = wetDepositionAer * chemParams['molar_mass']
    dep_wet_air_fw = wetDepositionAer2FW * chemParams['molar_mass']
    dep_wet_air_sw = wetDepositionAer2SW * chemParams['molar_mass']
    dep_wet_air_soil1 = wetDepositionAer2Soil1 * chemParams['molar_mass']
    dep_wet_air_soil2 = wetDepositionAer2Soil2 * chemParams['molar_mass']
    dep_wet_air_soil3 = wetDepositionAer2Soil3 * chemParams['molar_mass']
    dep_wet_air_soil4 = wetDepositionAer2Soil4 * chemParams['molar_mass']

    rain_dis_air = rainDissolutionAer * chemParams['molar_mass']
    rain_dis_air_fw = rainDissolutionAer2FW * chemParams['molar_mass']
    rain_dis_air_sw = rainDissolutionAer2SW * chemParams['molar_mass']
    rain_dis_air_soil1 = rainDissolutionAer2Soil1 * chemParams['molar_mass']
    rain_dis_air_soil2 = rainDissolutionAer2Soil2 * chemParams['molar_mass']
    rain_dis_air_soil3 = rainDissolutionAer2Soil3 * chemParams['molar_mass']
    rain_dis_air_soil4 = rainDissolutionAer2Soil4 * chemParams['molar_mass']

    dep_fSS = sedDepFWSusSed * chemParams['molar_mass']
    dep_sSS = sedDepSWSusSed * chemParams['molar_mass']

    # 4) diffusion process
    diff_air_fw = diffusionAir2FW * chemParams['molar_mass']
    diff_air_sw = diffusionAir2SW * chemParams['molar_mass']
    diff_fw_fSedW = diffusionFW2FWSed * chemParams['molar_mass']
    diff_sw_sSedW = diffusionSW2SWSed * chemParams['molar_mass']
    diff_air_soil1 = diffusionAir2S1surf * chemParams['molar_mass']
    diff_air_soil2 = diffusionAir2S2surf * chemParams['molar_mass']
    diff_air_soil3 = diffusionAir2S3surf * chemParams['molar_mass']
    diff_air_soil4 = diffusionAir2S4surf * chemParams['molar_mass']

    diff_fw_air = diffusionFW2Air * chemParams['molar_mass']
    diff_sw_air = diffusionSW2Air * chemParams['molar_mass']
    diff_fSedW_fw = diffusionFWSed2FW * chemParams['molar_mass']
    diff_sSedW_sw = diffusionSWSed2SW * chemParams['molar_mass']
    diff_soil1_air = diffusionS1surf2Air * chemParams['molar_mass']
    diff_soil2_air = diffusionS2surf2Air * chemParams['molar_mass']
    diff_soil3_air = diffusionS3surf2Air * chemParams['molar_mass']
    diff_soil4_air = diffusionS4surf2Air * chemParams['molar_mass']

    # 5) other process
    burial_fwSed = sedBurialFSed * chemParams['molar_mass']
    burial_swSed = sedBurialSSed * chemParams['molar_mass']
    resusp_fwSed = resuspFWSed * chemParams['molar_mass']
    resusp_swSed = resuspSWSed * chemParams['molar_mass']
    runoff_soil1 = runoffS1Water * chemParams['molar_mass']
    runoff_soil2 = runoffS2Water * chemParams['molar_mass']
    runoff_soil3 = runoffS3Water * chemParams['molar_mass']
    runoff_soil4 = runoffS4Water * chemParams['molar_mass']
    erosion_soil1 = erosionS1Solid * chemParams['molar_mass']
    erosion_soil2 = erosionS2Solid * chemParams['molar_mass']
    erosion_soil3 = erosionS3Solid * chemParams['molar_mass']
    erosion_soil4 = erosionS4Solid * chemParams['molar_mass']
    wind_erosion_soil1 = windErosionS1Solid * chemParams['molar_mass']
    wind_erosion_soil2 = windErosionS2Solid * chemParams['molar_mass']
    wind_erosion_soil3 = windErosionS3Solid * chemParams['molar_mass']
    wind_erosion_soil4 = windErosionS4Solid * chemParams['molar_mass']
    infiltra_soil1 = infiltrationSoil1surf * chemParams['molar_mass']
    infiltra_soil2 = infiltrationSoil2surf * chemParams['molar_mass']
    infiltra_soil3 = infiltrationSoil3surf * chemParams['molar_mass']
    infiltra_soil4 = infiltrationSoil4surf * chemParams['molar_mass']
    leach_soil1 = leachingS1deep * chemParams['molar_mass']
    leach_soil2 = leachingS2deep * chemParams['molar_mass']
    leach_soil3 = leachingS3deep * chemParams['molar_mass']
    leach_soil4 = leachingS4deep * chemParams['molar_mass']

    processes = [adv_air_in, adv_air_out, adv_fw_in, adv_fw_out, adv_fwSed_in, adv_fwSed_out, adv_sw_in, adv_sw_out,
                 adv_swSed_out, dep_dry_air, dep_dry_air_fw, dep_dry_air_sw, dep_dry_air_soil1, dep_dry_air_soil2,
                 dep_dry_air_soil3, dep_dry_air_soil4, dep_wet_air, dep_wet_air_fw, dep_wet_air_sw, dep_wet_air_soil1,
                 dep_wet_air_soil2, dep_wet_air_soil3, dep_wet_air_soil4, rain_dis_air, rain_dis_air_fw,
                 rain_dis_air_sw, rain_dis_air_soil1, rain_dis_air_soil2, rain_dis_air_soil3, rain_dis_air_soil4, dep_fSS, dep_sSS,
                 diff_air_fw, diff_air_sw, diff_fw_fSedW, diff_sw_sSedW, diff_air_soil1, diff_air_soil2, diff_air_soil3,
                 diff_air_soil4, diff_fw_air, diff_sw_air, diff_fSedW_fw, diff_sSedW_sw, diff_soil1_air, diff_soil2_air,
                 diff_soil3_air, diff_soil4_air, burial_fwSed, burial_swSed, resusp_fwSed, resusp_swSed,
                 runoff_soil1, runoff_soil2, runoff_soil3, runoff_soil4, erosion_soil1, erosion_soil2, erosion_soil3,
                 erosion_soil4, wind_erosion_soil1, wind_erosion_soil2, wind_erosion_soil3, wind_erosion_soil4,
                 infiltra_soil1, infiltra_soil2, infiltra_soil3, infiltra_soil4, leach_soil1, leach_soil2, leach_soil3,
                 leach_soil4, deg_air, deg_fw, deg_fwSed, deg_sw, deg_swSed, deg_soil1, deg_deepS1, deg_soil2,
                 deg_deepS2, deg_soil3, deg_deepS3, deg_soil4, deg_deepS4]

    return processes

