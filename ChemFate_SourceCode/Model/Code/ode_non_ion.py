from __future__ import division
from degradation_process import Degradation
from advective_processes import AdvectiveProcess
from diffusion_process_non_ion import Diffusion, MTC
from Z_non_ion import zValue


def org_ode(t, f, i, presence, env, climate, chemParams, release, bgConc):
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
    dryDepositionAer = adv.D_dep_dry(kAerDryDep, env['area']*env['aerVf'], zAerSub) * f[0]
    dryDepositionAer2FW = dryDepositionAer * (env['freshwA'] / env['area'])
    dryDepositionAer2SW = dryDepositionAer * (env['seawA'] / env['area'])
    dryDepositionAer2Soil1 = dryDepositionAer * (env['soilA1'] / env['area'])
    dryDepositionAer2Soil2 = dryDepositionAer * (env['soilA2'] / env['area'])
    dryDepositionAer2Soil3 = dryDepositionAer * (env['soilA3'] / env['area'])
    dryDepositionAer2Soil4 = dryDepositionAer * (env['soilA4'] / env['area'])

    # aerosol wet deposition
    wetDepositionAer = adv.D_dep_wet(climate['precip_m'][i], env['scavenging'], env['area']*env['aerVf'], zAerSub) * env['aerVf'] * f[0]
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
    sedDepFWSusSed = adv.D_dep_dry(kFWSedDep, env['freshwA']* env['fSSVf'], zFWSusSedSub) * f[1]


    # freshwater sediment f[2]
    # freshwater sediment degradation
    dDegFWSed = deg.D_deg(env['fSedSV'], chemParams['kDeg_fSedS_n'], zFSedSSub) + \
                deg.D_deg(env['fSedWV'], chemParams['kDeg_fSedW_n'], zWaterSub)

    degradationFSed = dDegFWSed * f[2]

    # freshwater sediment advection
    advectionFWSedInflow = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], bgConc['gfSedc_n'])
    advectionFSedOutflow = adv.D_advec_water(climate['waterflow_d'][i] * env['fwadvfrac'], zFWSedimentBulk) * f[2]
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
    advectionSSedOutflow = adv.D_advec_water(climate['waterflow_d'][i] * env['swadvfrac'], 10 * zSWSedimentBulk) * f[4]
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

    ###################################################################################################################
    # the code below will combine the processes together 1) entering the compartment, 2) loss from the
    # compartment, and 3) transfer between compartments
    # if one compartment doesn't exist, re-route scenarios are considered
    ###################################################################################################################

    # air f[0]
    AirEnter = release['air'][i] + advectionAirInflow  # mol/day

    diffFW2Air = diffusionFW2Air
    diffSW2Air = diffusionSW2Air
    diffS1surf2Air = diffusionS1surf2Air
    diffS2surf2Air = diffusionS2surf2Air
    diffS3surf2Air = diffusionS3surf2Air
    diffS4surf2Air = diffusionS4surf2Air

    # without air compartment, all of the air related processes would be 0
    # the transfer process from other compartments still exist, but the transfer to this compartment is 0
    if presence['air'] == 0:
        advectionAirInflow, degradationAir, advectionAirBulkOutflow, rainDissolutionAer, \
        dryDepositionAer, wetDepositionAer, diffusionAir2FW, diffusionAir2SW, diffusionAir2S1surf, diffusionAir2S2surf, \
        diffusionAir2S3surf, diffusionAir2S4surf, diffFW2Air, diffSW2Air, diffS1surf2Air, diffS2surf2Air, \
        diffS3surf2Air, diffS4surf2Air, windErosionS1Solid, windErosionS2Solid, windErosionS3Solid, windErosionS4Solid \
            = [0] * 22

    AirLoss = -(degradationAir + advectionAirBulkOutflow + rainDissolutionAer + dryDepositionAer +
                wetDepositionAer + diffusionAir2FW + diffusionAir2SW + diffusionAir2S1surf + diffusionAir2S2surf +
                diffusionAir2S3surf + diffusionAir2S4surf)

    AirTransfer = diffFW2Air + diffSW2Air + diffS1surf2Air + diffS2surf2Air + diffS3surf2Air + diffS4surf2Air + \
                  windErosionS1Solid + windErosionS2Solid + windErosionS3Solid + windErosionS4Solid

    # freshwater f[1]
    FWEnter = release['fw'][i] + release['fSS'][i] + advectionFWBulkInflow

    rainDisAir2FW = rainDissolutionAer2FW
    dryDepAer2FW = dryDepositionAer2FW
    wetDepAer2FW = wetDepositionAer2FW
    runoffS1 = runoffS1Water
    runoffS2 = runoffS2Water
    runoffS3 = runoffS3Water
    runoffS4 = runoffS4Water
    erosionS1 = erosionS1Solid
    erosionS2 = erosionS2Solid
    erosionS3 = erosionS3Solid
    erosionS4 = erosionS4Solid
    diffAir2FW = diffusionAir2FW

    # without freshwater compartment, freshwater related processes would be 0
    if presence['fw'] == 0:
        advectionFWBulkInflow, diffusionFW2Air, diffusionFW2FWSed, degradationFW, sedDepFWSusSed, advectionFWBulkOutflow, \
        rainDisAir2FW, dryDepAer2FW, wetDepAer2FW, runoffS1, runoffS2, runoffS3, runoffS4, erosionS1, erosionS2, erosionS3, \
        erosionS4, diffAir2FW, diffusionFWSed2FW, resuspFWSed = [0] * 20

    FWLoss = - (diffusionFW2Air + diffusionFW2FWSed + degradationFW + sedDepFWSusSed + advectionFWBulkOutflow)
    FWTransfer = rainDisAir2FW + dryDepAer2FW + wetDepAer2FW \
                 + runoffS1 + runoffS2 + runoffS3 + runoffS4 \
                 + erosionS1 + erosionS2 + erosionS3 + erosionS4 \
                 + diffAir2FW + diffusionFWSed2FW + resuspFWSed


    # freshwater sediment f[2]
    FSedEnter = release['fwSed'][i] + advectionFWSedInflow
    diffFW2FWSed = diffusionFW2FWSed
    sedDepFWSusSed2FSed = sedDepFWSusSed
    # no freshwater, no freshwater sediment, freshwater sediment related processes would be 0
    if presence['fw'] == 0:
        advectionFWSedInflow, resuspFWSed, degradationFSed, diffusionFWSed2FW, advectionFSedOutflow, sedBurialFSed, \
        diffFW2FWSed, sedDepFWSusSed2FSed = [0] * 8

    FSedLoss = - (resuspFWSed + degradationFSed + diffusionFWSed2FW + advectionFSedOutflow + sedBurialFSed)
    FSedTransfer = diffFW2FWSed + sedDepFWSusSed2FSed

    # seawater f[3]
    SWEnter = release['sw'][i] + release['sSS'][i]
    rainDisAir2SW = rainDissolutionAer2SW
    dryDepAer2SW = dryDepositionAer2SW
    wetDepAer2SW = wetDepositionAer2SW
    advectFWOutflow = advectionFWOutflow
    diffAir2SW = diffusionAir2SW

    runoffS1SW, runoffS2SW, runoffS3SW, runoffS4SW, erosionS1SW, erosionS2SW, erosionS3SW, erosionS4SW = [0] * 8

    if presence['fw'] == 0:
        runoffS1SW = runoffS1Water
        runoffS2SW = runoffS2Water
        runoffS3SW = runoffS3Water
        runoffS4SW = runoffS4Water
        erosionS1SW = erosionS1Solid
        erosionS2SW = erosionS2Solid
        erosionS3SW = erosionS3Solid
        erosionS4SW = erosionS4Solid

    if presence['sw'] == 0:
        diffusionSW2Air, diffusionSW2SWSed, degradationSW, sedDepSWSusSed, advectionSWBulkOutflow, rainDisAir2SW, \
        dryDepAer2SW, wetDepAer2SW, advectFWOutflow, diffAir2SW, diffusionSWSed2SW, resuspSWSed = [0] * 12

    SWLoss = - (diffusionSW2Air + diffusionSW2SWSed + degradationSW + sedDepSWSusSed + advectionSWBulkOutflow)
    SWTransfer = runoffS1SW + runoffS2SW + runoffS3SW + runoffS4SW + erosionS1SW + erosionS2SW + erosionS3SW + \
                 erosionS4SW + rainDisAir2SW + dryDepAer2SW + wetDepAer2SW + advectFWOutflow \
                 + diffAir2SW + diffusionSWSed2SW + resuspSWSed

    # seawater sediment f[4]
    SSedEnter = release['swSed'][i]
    if presence['sw'] == 0:
        resuspSWSed, degradationSSed, diffusionSWSed2SW, advectionSSedOutflow, diffusionSW2SWSed, sedDepSWSusSed = [0] * 6

    SSedLoss = - (resuspSWSed + degradationSSed + diffusionSWSed2SW + advectionSSedOutflow + sedBurialSSed)
    SSedTransfer = diffusionSW2SWSed + sedDepSWSusSed

    # soil 1 surface soil f[5]
    S1surfEnter = release['soil1'][i]
    diffAir2S1surf = diffusionAir2S1surf
    rainDisAir2Soil1 = rainDissolutionAer2Soil1
    dryDepAer2Soil1 = dryDepositionAer2Soil1
    wetDepAer2Soil1 = wetDepositionAer2Soil1

    if presence['soil1'] == 0:
        diffusionS1surf2Air, degradationSoil1surf, runoffS1Water, erosionS1Solid, infiltrationSoil1surf, \
        diffAir2S1surf, rainDisAir2Soil1, dryDepAer2Soil1, wetDepAer2Soil1, windErosionS1Solid = [0] * 10

    S1surfLoss = - (diffusionS1surf2Air + degradationSoil1surf + runoffS1Water + erosionS1Solid +
                    infiltrationSoil1surf + windErosionS1Solid)
    S1surfTransfer = diffAir2S1surf + rainDisAir2Soil1 + dryDepAer2Soil1 + wetDepAer2Soil1

    # soil 1 deep soil f[6]
    S1DeepSEnter = release['dsoil1'][i]
    if presence['soil1'] == 0:
        degradationSoil1deep, leachingS1deep, infiltrationSoil1surf = 0
    S1DeepSLoss = - degradationSoil1deep - leachingS1deep
    S1DeepSTransfer = infiltrationSoil1surf

    # soil 2 surface soil f[7]
    S2surfEnter = release['soil2'][i]
    diffAir2S2surf = diffusionAir2S2surf
    rainDisAir2Soil2 = rainDissolutionAer2Soil2
    dryDepAer2Soil2 = dryDepositionAer2Soil2
    wetDepAer2Soil2 = wetDepositionAer2Soil2
    if presence['soil2'] == 0:
        diffusionS2surf2Air, degradationSoil2surf, runoffS2Water, erosionS2Solid, infiltrationSoil2surf, \
        diffAir2S2surf, rainDisAir2Soil2, dryDepAer2Soil2, wetDepAer2Soil2, windErosionS2Solid = [0] * 10
    S2surfLoss = - (diffusionS2surf2Air + degradationSoil2surf + runoffS2Water + erosionS2Solid +
                    infiltrationSoil2surf + windErosionS2Solid)
    S2surfTransfer = diffAir2S2surf + rainDisAir2Soil2 + dryDepAer2Soil2 + wetDepAer2Soil2

    # soil 2 deep soil f[8]
    S2DeepSEnter = release['dsoil2'][i]
    if presence['soil2'] == 0:
        degradationSoil2deep, leachingS2deep, infiltrationSoil2surf = [0] * 3
    S2DeepSLoss = - degradationSoil2deep - leachingS2deep
    S2DeepSTransfer = infiltrationSoil2surf

    # soil 3 surface soil f[9]
    S3surfEnter = release['soil3'][i]
    diffAir2S3surf = diffusionAir2S3surf
    rainDisAir2Soil3 = rainDissolutionAer2Soil3
    dryDepAer2Soil3 = dryDepositionAer2Soil3
    wetDepAer2Soil3 = wetDepositionAer2Soil3
    if presence['soil3'] == 0:
        diffusionS3surf2Air, degradationSoil3surf, runoffS3Water, erosionS3Solid, infiltrationSoil3surf, \
        diffAir2S3surf, rainDisAir2Soil3, dryDepAer2Soil3, wetDepAer2Soil3, windErosionS3Solid = [0] * 10
    S3surfLoss = - (diffusionS3surf2Air + degradationSoil3surf + runoffS3Water + erosionS3Solid +
                    infiltrationSoil3surf + windErosionS3Solid)
    S3surfTransfer = diffAir2S3surf + rainDisAir2Soil3 + dryDepAer2Soil3 + wetDepAer2Soil3

    # soil 3 deep soil f[10]
    S3DeepSEnter = release['dsoil3'][i]
    if presence['soil3'] == 0:
        degradationSoil3deep, leachingS3deep, infiltrationSoil3surf = [0]*3

    S3DeepSLoss = - degradationSoil3deep - leachingS3deep
    S3DeepSTransfer = infiltrationSoil3surf

    # soil 4 surface soil f[11]
    S4surfEnter = release['soil4'][i]
    diffAir2S4surf = diffusionAir2S4surf
    rainDisAir2Soil4 = rainDissolutionAer2Soil4
    dryDepAer2Soil4 = dryDepositionAer2Soil4
    wetDepAer2Soil4 = wetDepositionAer2Soil4
    if presence['soil4'] == 0:
        diffusionS4surf2Air, degradationSoil4surf, runoffS4Water, erosionS4Solid, infiltrationSoil4surf, \
        diffAir2S4surf, rainDisAir2Soil4, dryDepAer2Soil4, wetDepAer2Soil4, windErosionS4Solid = [0] * 10
    S4surfLoss = - (diffusionS4surf2Air + degradationSoil4surf + runoffS4Water + erosionS4Solid +
                    infiltrationSoil4surf + windErosionS4Solid)
    S4surfTransfer = diffAir2S4surf + rainDisAir2Soil4 + dryDepAer2Soil4 + wetDepAer2Soil4

    # soil 4 deep soil f[12]
    S4DeepSEnter = release['dsoil4'][i]
    if presence['soil4'] == 0:
        degradationSoil4deep, leachingS4deep, infiltrationSoil4surf = [0] * 3
    S4DeepSLoss = - degradationSoil4deep - leachingS4deep
    S4DeepSTransfer = infiltrationSoil4surf

    # unit: mol/day / (m3 * mol/m3-Pa) = Pa/day
    if (env['areaV'] * zAirBulk) != 0:
        air = (AirEnter + AirLoss + AirTransfer) / (env['areaV'] * zAirBulk)
    else:
        air = 0

    if (env['fwV'] * zFWBulk) != 0:
        fw = (FWEnter + FWLoss + FWTransfer) / (env['fwV'] * zFWBulk)
    else:
        fw = 0

    if (env['sedFWV'] * zFWSedimentBulk) != 0:
        fwSed = (FSedEnter + FSedLoss + FSedTransfer) / (env['sedFWV'] * zFWSedimentBulk)
    else:
        fwSed = 0

    if (env['swV'] * zSWBulk) != 0:
        sw = (SWEnter + SWLoss + SWTransfer) / (env['swV'] * zSWBulk)
    else:
        sw = 0

    if (env['sedSWV'] * zSWSedimentBulk) != 0:
        swSed = (SSedEnter + SSedLoss + SSedTransfer) / (env['sedSWV'] * zSWSedimentBulk)
    else:
        swSed = 0

    if (env['soilV1'] * zSoil1Bulk) != 0:
        s1surf = (S1surfEnter + S1surfLoss + S1surfTransfer) / (env['soilV1'] * zSoil1Bulk)
    else:
        s1surf = 0

    if (env['deepSV1'] * zS1DeepSSub) != 0:
        s1deep = (S1DeepSEnter + S1DeepSLoss + S1DeepSTransfer) / (env['deepSV1'] * zS1DeepSSub)
    else:
        s1deep = 0

    if (env['soilV2'] * zSoil2Bulk) != 0:
        s2surf = (S2surfEnter + S2surfLoss + S2surfTransfer) / (env['soilV2'] * zSoil2Bulk)
    else:
        s2surf = 0

    if (env['deepSV2'] * zS2DeepSSub) != 0:
        s2deep = (S2DeepSEnter + S2DeepSLoss + S2DeepSTransfer) / (env['deepSV2'] * zS2DeepSSub)
    else:
        s2deep = 0

    if (env['soilV3'] * zSoil3Bulk) != 0:
        s3surf = (S3surfEnter + S3surfLoss + S3surfTransfer) / (env['soilV3'] * zSoil3Bulk)
    else:
        s3surf = 0

    if (env['deepSV3'] * zS3DeepSSub) != 0:
        s3deep = (S3DeepSEnter + S3DeepSLoss + S3DeepSTransfer) / (env['deepSV3'] * zS3DeepSSub)
    else:
        s3deep = 0

    if (env['soilV4'] * zSoil4Bulk) != 0:
        s4surf = (S4surfEnter + S4surfLoss + S4surfTransfer) / (env['soilV4'] * zSoil4Bulk)
    else:
        s4surf = 0

    if (env['deepSV4'] * zS4DeepSSub) != 0:
        s4deep = (S4DeepSEnter + S4DeepSLoss + S4DeepSTransfer) / (env['deepSV4'] * zS4DeepSSub)
    else:
        s4deep = 0

    dYdt = [air, fw, fwSed, sw, swSed, s1surf, s2deep, s2surf, s2deep, s3surf, s3deep, s4surf, s4deep]

    return dYdt
