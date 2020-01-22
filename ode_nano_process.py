import numpy as np

from advective_processes_nano import dryDepAir
from advective_processes_nano import wetDepAir
from advective_processes_nano import odeHetagg
from advective_processes_nano import airAdvection
from advective_processes_nano import dryDepAer
from advective_processes_nano import wetDep
from advective_processes_nano import waterAdv
from advective_processes_nano import odeDissolution
from advective_processes_nano import sedDeposition
from advective_processes_nano import resuspensionSed
from advective_processes_nano import burial
from advective_processes_nano import aerosolResuspension
from advective_processes_nano import windErosion
from advective_processes_nano import erosion
from advective_processes_nano import soilwaterPartition
from advective_processes_nano import runoff
from advective_processes_nano import vertFlow
from advective_processes_nano import horiFlow

#################################################################
#
#   CLiCC Nanomaterials F&T Model Developed by Dr. Garner
#	Model original devloped in MATLAB
#   Date: December 26th, 2016
#   Converted to Python by Jill Farnsworth
#
#################################################################

def nano_process(f,i,V,presence,env,climate,ENM,DIS):
	# %   Differential equation solver for ENM mass in all compartments
	# %   t is time, f is the mass by compartment and day, i is the iteration in
	# %   the for loop - so the time step, V is the volume vector
	# Note: i in Python should be one less than i in Matlab, since Python is zero indexed

	# %% Process calculations

	# % Air
	if presence['air']==1:
		# % dry deposition from air
		dryDepositionAir = dryDepAir(ENM['density'],env['airP'],env['dynViscAir'],(ENM['radiusENMagg']*(10**-9)),env['airH'])*f[0]
		# % wet deposition from air
		wetDepositionAir = wetDepAir(climate['precip'][i],env['scavengingENM'],env['area'],V[0])*f[0]
		# % heteroaggregation in air (pseudofirst order rate constant)
		heteroaggregationAirAer = odeHetagg(ENM['khetA'],env['aerC'])*f[0]
		# % advection in air
		advectionAir = airAdvection(climate['windspeed'][i],env['area'],env['airH'],V[0],np.true_divide(f[0],V[0]),ENM['density'])*f[0]
	else:
		dryDepositionAir = 0
		wetDepositionAir = 0
		heteroaggregationAirAer = 0
		advectionAir = 0

	# % Aerosols
	if presence['aer']==1:
		# % dry deposition of aerosols
		dryDepositionAer = dryDepAer(env['aerP'],env['airP'],env['dynViscAir'],env['radiusParticlesAer'],env['airH'])*f[1]
		# % wet deposition of aerosols
		wetDepositionAer = wetDep(climate['precip'][i],V[0],env['scavenging'],env['area'])*f[1]
		# % advection of aerosols
		advectionAer = airAdvection(climate['windspeed'][i],env['area'],env['airH'],V[0],env['aerC'],env['aerP'])*f[1]
	else:
		dryDepositionAer = 0
		wetDepositionAer = 0
		advectionAer = 0

	# % Freshwater
	if presence['fw']==1:
		# % sedimentation of free ENMs
		sedimentationFW = np.true_divide(ENM['ksedFW'],env['freshwD'])*f[2]
		# % heteroaggregation in freshwater with suspended sediment
		heteroaggregationFW = odeHetagg(ENM['khetFW'],env['freshssC'])*f[2]
		# % advective flow
		advectionFW = waterAdv(climate['flow'][i],V[2],np.true_divide(f[2],V[2]),ENM['density'])*f[2]
		# % Dissolution in freshwater
		dissolutionFW = odeDissolution(DIS['percfitaFW'],DIS['percfitbFW'],f[2],f[16],ENM['kdisFW'],V[2])
	else:
		sedimentationFW = 0
		heteroaggregationFW = 0
		advectionFW = 0
		dissolutionFW = 0

	# % Freshwater Suspended Sediment
	if presence['fSS']==1:
		# % deposition of freshwater suspended sediment
		sedimentationFWSS = sedDeposition(env['freshssP'],env['freshwP'],env['dynViscFW'],env['radiusParticlesFW'],env['freshwD'])*f[3]
		# % advection of freshwater suspended sediment
		advectionFWSS = waterAdv(climate['flow'][i],V[2],env['freshssC'],env['freshssP'])*f[3]
	else:
		sedimentationFWSS = 0
		advectionFWSS = 0

	# % Freshwater Sediment
	if presence['fSed']==1:
		# % resuspension from freshwater sediment
		resuspensionFWSed = resuspensionSed(env['sedFWA'],env['resuspensionRateFW'],V[4])*f[4]
		# % burial in freshwater sediment
		burialFWSed = burial(env['sedFWA'],env['burialRateFW'],V[4])*f[4]
		# % advection of sediment
		advectionFWSed = waterAdv(climate['flow'][i],V[4],f[4]/V[4],env['sedFWP'])*env['fwadvfrac']*f[4]
		# % dissolution in sediment
		dissolutionFWSed = odeDissolution(DIS['percfitaFW'],DIS['percfitbFW'],f[4],f[17],ENM['kdisFWsed'],V[4])
	else:
		resuspensionFWSed = 0
		burialFWSed = 0
		advectionFWSed = 0
		dissolutionFWSed = 0

	# % Seawater
	if presence['sw']==1:
		# % sedimentation of free ENMs in seawater
		sedimentationSW = (np.true_divide(ENM['ksedSW'],env['seawD']))*f[5]
		# % heteroaggregation with suspended sediment in seawater
		heteroaggregationSW = odeHetagg(ENM['khetSW'],env['seassC'])*f[5]
		# % aerosolization of particles from seawater to aerosols or air
		aerosolizationSW = aerosolResuspension(climate['windspeed'][i],env['coastalA'],ENM['enrichFactor'],env['seawD'],env['seawV'])*f[5]
		# % advection out of the coastal marine to larger system
		advectionSW = waterAdv(climate['flow'][i],V[5],np.true_divide(f[5],V[5]),ENM['density'])*f[5]
		# % dissolution in marine
		dissolutionSW = odeDissolution(DIS['percfitaSW'],DIS['percfitbSW'],f[5],f[18],ENM['kdisSW'],V[5])
	else:
		sedimentationSW = 0
		heteroaggregationSW = 0
		aerosolizationSW = 0
		advectionSW = 0
		dissolutionSW = 0

	# % Seawater suspended sediment
	if presence['sSS']==1:
		# % deposition of suspended sediment
		sedimentationSWSS = sedDeposition(env['seassP'],env['seawP'],env['dynViscSW'],env['radiusParticlesSW'],env['seawD'])*f[6]
		# % advective flow
		advectionSWSS = waterAdv(climate['flow'][i],V[5],env['seassC'],env['seassP'])*f[6]
	else:
		sedimentationSWSS = 0
		advectionSWSS = 0

	# % Seawater sediment
	if presence['sSed']==1:
		# % resuspension from marine sediment
		resuspensionSWSed = resuspensionSed(env['sedSWA'],env['resuspensionRateSW'],V[7])*f[7]
		# % burial in marine sediment
		burialSWSed = burial(env['sedSWA'],env['burialRateSW'],V[7])*f[7]
		# % advection of sediment
		advectionSWSed = waterAdv(climate['flow'][i],V[7],f[7]/V[7],env['sedSWP'])*env['swadvfrac']*f[7]
		# % dissolution in sediment
		dissolutionSWSed = odeDissolution(DIS['percfitaSW'],DIS['percfitbSW'],f[7],f[19],ENM['kdisSWsed'],V[7])
	else:
		resuspensionSWSed = 0
		burialSWSed = 0
		advectionSWSed = 0
		dissolutionSWSed = 0

	# % Soil 1
	if presence['soil1']==1:
		# % wind erosion
		windErosionSoil1 = windErosion(climate['windspeed'][i],climate['precip'][i], \
									   env['roughness1'],env['Kconstant1'],env['airP'],env['soilA1'],\
									   env['A1'],env['TSV1'],env['TSVmin1'],env['z_wind1'],env['percWind1'],\
									   env['windConstant1'],env['percUncovered1'],env['percSuspended1'],env['soilP1'],V[8])*f[8]
		# % solid soil erosion
		solidErosionSoil1 = erosion(climate['precip'][i],env['Kfact1'],env['lenslope1'],\
									env['cropManageFactor1'],env['supportFactor1'],env['soilA1'],env['soilP1'])*np.true_divide(f[8],V[8])
		# % loss by partitioning to soil water
		# [0] calls first value in function
		soil2soilwater1 = soilwaterPartition(f[8],f[9],ENM['elutionS1'],1)[0]
	else:
		windErosionSoil1 = 0
		solidErosionSoil1 = 0
		soil2soilwater1 = 0

	# % Soil Water 1
	if presence['soilW1']==1:
		# % runoff
		runoffSoil1 = runoff(climate['precip'][i],env['CN1'],env['soilA1'], V[8])*f[9]
		# % infiltra
		k_infil1, infiltraSoil1 = vertFlow(climate['precip'][i], env['CN1'], climate['evap'][i], env['FC1'], env['soilWC1'],
								 env['soilV1'], env['soilA1'])
		infiltraSoil1 = infiltraSoil1*f[9]
		# % loss to partitioning to soil solids
		# [1] calls second value in function
		soilwater2soil1 = soilwaterPartition(f[8],f[9],ENM['elutionS1'],2)[1]
		# % dissolution
		dissolutionSoil1 = odeDissolution(DIS['percfitaGW1'],DIS['percfitbGW1'],f[9],f[20],ENM['kdisS1'],V[9])
	else:
		runoffSoil1 = 0
		infiltraSoil1 = 0
		soilwater2soil1 = 0
		dissolutionSoil1 = 0

	# % Soil 2
	if presence['soil2']==1:
		# % wind erosion
		windErosionSoil2 = windErosion(climate['windspeed'][i],climate['precip'][i], env['roughness2'],\
									   env['Kconstant2'],env['airP'],env['soilA2'],env['A2'],env['TSV2'],\
									   env['TSVmin2'],env['z_wind2'],env['percWind2'],env['windConstant2'],\
									   env['percUncovered2'],env['percSuspended2'],env['soilP2'],V[10])*f[10]
		# % solid soil erosion
		solidErosionSoil2 = erosion(climate['precip'][i],env['Kfact2'],env['lenslope2'],\
									env['cropManageFactor2'],env['supportFactor2'],env['soilA2'],env['soilP2'])*np.true_divide(f[10],V[10])
		# % loss by partitioning to soil water
		soil2soilwater2 = soilwaterPartition(f[10],f[11],ENM['elutionS2'],1)[0]
	else:
		windErosionSoil2 = 0
		solidErosionSoil2 = 0
		soil2soilwater2 = 0

	# % Soil Water 2
	if presence['soilW2']==1:
		# % runoff
		runoffSoil2 = runoff(climate['precip'][i],env['CN2'],env['soilA2'],V[10])*f[11]
		# % infiltra
		k_infil2, infiltraSoil2 = vertFlow(climate['precip'][i], env['CN2'], climate['evap'][i], env['FC2'], env['soilWC2'],
								 env['soilV2'], env['soilA2'])
		infiltraSoil2 = infiltraSoil2*f[11]
		# % loss to partitioning to soil solids
		soilwater2soil2 = soilwaterPartition(f[10],f[11],ENM['elutionS2'],2)[1]
		# % dissolution
		dissolutionSoil2 = odeDissolution(DIS['percfitaGW2'],DIS['percfitbGW2'],f[11],f[21],ENM['kdisS2'],V[11])
	else:
		runoffSoil2 = 0
		infiltraSoil2 = 0
		soilwater2soil2 = 0
		dissolutionSoil2 = 0

	# % Soil 3
	if presence['soil3']==1:
		# % wind erosion
		windErosionSoil3 = windErosion(climate['windspeed'][i],climate['precip'][i], env['roughness3'],\
									   env['Kconstant3'],env['airP'],env['soilA3'],env['A3'],env['TSV3'],\
									   env['TSVmin3'],env['z_wind3'],env['percWind3'],env['windConstant3'],\
									   env['percUncovered3'],env['percSuspended3'],env['soilP3'],V[12])*f[12]
		# % solid soil erosion
		solidErosionSoil3 = erosion(climate['precip'][i],env['Kfact3'],env['lenslope3'],\
									env['cropManageFactor3'],env['supportFactor3'],env['soilA3'],env['soilP3'])*np.true_divide(f[12],V[12])
		# % loss by partitioning to soil water
		soil2soilwater3 = soilwaterPartition(f[12],f[13],ENM['elutionS3'],1)[0]
	else:
		windErosionSoil3 = 0
		solidErosionSoil3 = 0
		soil2soilwater3 = 0

	# % Soil Water 3
	if presence['soilW3']==1:
		# % runoff
		runoffSoil3 = runoff(climate['precip'][i],env['CN3'],env['soilA3'],V[12])*f[13]
		# % infiltra
		k_infil3, infiltraSoil3 = vertFlow(climate['precip'][i], env['CN3'], climate['evap'][i], env['FC3'], env['soilWC3'],
								 env['soilV3'], env['soilA3'])
		infiltraSoil3 = infiltraSoil3*f[13]
		# % loss to partitioning to soil solids
		soilwater2soil3 = soilwaterPartition(f[12],f[13],ENM['elutionS3'],2)[1]
		# % dissolution
		dissolutionSoil3 = odeDissolution(DIS['percfitaGW3'],DIS['percfitbGW3'],f[13],f[22],ENM['kdisS3'],V[13])
	else:
		runoffSoil3 = 0
		infiltraSoil3 = 0
		soilwater2soil3 = 0
		dissolutionSoil3 = 0

	# % Soil 4
	if presence['soil4']==1:
		# % wind erosion
		windErosionSoil4 = windErosion(climate['windspeed'][i],climate['precip'][i], env['roughness4'],env['Kconstant4'],\
									   env['airP'],env['soilA4'],env['A4'],env['TSV4'],env['TSVmin4'],env['z_wind4'],\
									   env['percWind4'],env['windConstant4'],env['percUncovered4'],env['percSuspended4'],\
									   env['soilP4'],V[14])*f[14]
		# % solid soil erosion
		solidErosionSoil4 = erosion(climate['precip'][i],env['Kfact4'],env['lenslope4'],\
									env['cropManageFactor4'],env['supportFactor4'],env['soilA4'],env['soilP4'])*np.true_divide(f[14],V[14])
		# % loss by partitioning to soil water
		soil2soilwater4 = soilwaterPartition(f[14],f[15],ENM['elutionS4'],1)[0]
	else:
		windErosionSoil4 = 0
		solidErosionSoil4 = 0
		soil2soilwater4 = 0

	# % Soil Water 4
	if presence['soilW4']==1:
		# % runoff
		runoffSoil4 = runoff(climate['precip'][i],env['CN4'],env['soilA4'],V[14])*f[15]
		# % infiltra
		k_infil4, infiltraSoil4 = vertFlow(climate['precip'][i], env['CN4'], climate['evap'][i], env['FC4'], env['soilWC4'],
								 env['soilV4'], env['soilA4'])
		infiltraSoil4 = infiltraSoil4*f[15]
		# % loss to partitioning to soil solids
		soilwater2soil4 = soilwaterPartition(f[14],f[15],ENM['elutionS4'],2)[1]
		# % dissolution
		dissolutionSoil4 = odeDissolution(DIS['percfitaGW4'],DIS['percfitbGW4'],f[15],f[23],ENM['kdisS4'],V[15])
	else:
		runoffSoil4 = 0
		infiltraSoil4 = 0
		soilwater2soil4 = 0
		dissolutionSoil4 = 0
	
	# % Dissolved FW mass transport
	if presence['fw'] == 1:
		advectionFWDis = waterAdv(climate['flow'][i],V[2],f[16]/V[2],env['freshwP'])*f[16]
	else:
		advectionFWDis = 0
	
	# % Dissolved FW sediment mass transport
	if presence['fw'] == 1:
		advectionFWSedDis = waterAdv(climate['flow'][i],V[4],f[17]/V[4],env['sedFWP'])*env['fwadvfrac']*f[17]
	else:
		advectionFWSedDis = 0
		
	# % Dissolved SW mass transport
	if presence['sw'] == 1:
		advectionSWDis = waterAdv(climate['flow'][i],V[5],f[18]/V[5],env['seawP'])*f[18]
	else:
		advectionSWDis = 0
	
	# % Dissolved SW sediment mass transport
	if presence['sw'] == 1:
		advectionSWSedDis = waterAdv(climate['flow'][i],V[7],f[19]/V[7],env['sedSWP'])*env['fwadvfrac']*f[19]
	else:
		advectionSWSedDis = 0
	
	# % Dissolved Soil Water 1 runoff
	if presence['soil1'] == 1:
		runoffSoilDis1 = runoff(climate['precip'][i],env['CN1'],env['soilA1'],V[8])*f[20]
	else:
		runoffSoilDis1 = 0
	
	# % Dissolved Soil Water 2 runoff
	if presence['soil2'] == 1:
		runoffSoilDis2 = runoff(climate['precip'][i],env['CN2'],env['soilA2'],V[10])*f[21]
	else:
		runoffSoilDis2 = 0
	
	# % Dissolved Soil Water 3 runoff
	if presence['soil3'] == 1:
		runoffSoilDis3 = runoff(climate['precip'][i],env['CN3'],env['soilA3'],V[12])*f[22]
	else:
		runoffSoilDis3 = 0
	
	# % Dissolved Soil Water 4 runoff
	if presence['soil4'] == 1:
		runoffSoilDis4 = runoff(climate['precip'][i],env['CN4'],env['soilA4'],V[14])*f[23]
	else:
		runoffSoilDis4 = 0

	# deep soil 1
	if presence['soil1'] == 1:
		leachSoil1 = horiFlow(k_infil1, V[24]) * f[24]
	else:
		leachSoil1 = 0

	# deep soil 2
	if presence['soil2'] == 1:
		leachSoil2 = horiFlow(k_infil2, V[25]) * f[25]
	else:
		leachSoil2 = 0

	# deep soil 3
	if presence['soil3'] == 1:
		leachSoil3 = horiFlow(k_infil3, V[26]) * f[26]
	else:
		leachSoil3 = 0

	# deep soil 4
	if presence['soil4'] == 1:
		leachSoil4 = horiFlow(k_infil4, V[27]) * f[27]
	else:
		leachSoil4 = 0

	###################################################################
	# processes output to transport rate kg/day
	# it's already in ng/day
	# divide by (10 ** 9) to get kg/day
	###################################################################

	# 1) advection process
	adv_air_in = (advectionAir + advectionAer)/(10 ** 9)
	adv_fw_out = (advectionFW + advectionFWSS)/(10 ** 9)
	adv_fwSed_out = advectionFWSed/(10 ** 9)
	adv_sw_out = (advectionSW + advectionSWSS)/(10 ** 9)
	adv_swSed_out = advectionSWSed/(10 ** 9)

	adv_fw_dissolved = advectionFWDis/(10 ** 9)
	adv_fwSed_dissolved = advectionFWSedDis/(10 ** 9)
	adv_sw_dissolved = advectionSWDis/(10 ** 9)
	adv_swSed_dissolved = advectionSWSedDis/(10 ** 9)

	# 2) deposition process
	dep_dry_air = (dryDepositionAir + dryDepositionAer)/(10 ** 9)
	dep_dry_air_fw = dep_dry_air * np.true_divide(env['freshwA'],env['area'])/(10 ** 9)
	dep_dry_air_sw = dep_dry_air * np.true_divide(env['seawA'],env['area'])/(10 ** 9)
	dep_dry_air_soil1 = dep_dry_air * np.true_divide(env['soilA1'],env['area'])/(10 ** 9)
	dep_dry_air_soil2 = dep_dry_air * np.true_divide(env['soilA2'],env['area'])/(10 ** 9)
	dep_dry_air_soil3 = dep_dry_air * np.true_divide(env['soilA3'],env['area'])/(10 ** 9)
	dep_dry_air_soil4 = dep_dry_air * np.true_divide(env['soilA4'],env['area'])/(10 ** 9)

	dep_wet_air = (wetDepositionAir + wetDepositionAer)/(10 ** 9)
	dep_wet_air_fw = dep_wet_air * np.true_divide(env['freshwA'],env['area'])/(10 ** 9)
	dep_wet_air_sw = dep_wet_air * np.true_divide(env['seawA'],env['area'])/(10 ** 9)
	dep_wet_air_soil1 = dep_wet_air * np.true_divide(env['soilA1'],env['area'])/(10 ** 9)
	dep_wet_air_soil2 = dep_wet_air * np.true_divide(env['soilA2'],env['area'])/(10 ** 9)
	dep_wet_air_soil3 = dep_wet_air * np.true_divide(env['soilA3'],env['area'])/(10 ** 9)
	dep_wet_air_soil4 = dep_wet_air * np.true_divide(env['soilA4'],env['area'])/(10 ** 9)

	dep_fSS = (sedimentationFW + sedimentationFWSS)/(10 ** 9)
	dep_sSS = (sedimentationSW + sedimentationSWSS)/(10 ** 9)

	# 3) heteroaggregation process
	heteroagg_air = heteroaggregationAirAer/(10 ** 9)
	heteroagg_fw = heteroaggregationFW/(10 ** 9)
	heteroagg_sw = heteroaggregationSW/(10 ** 9)

	# 4) dissolution process
	# dissolutionFW may return an array [] or just 0
	try:
		dissolution_fw = dissolutionFW[0]/(10 ** 9)
	except:
		dissolution_fw = dissolutionFW/(10 ** 9)
	try:
		dissolution_fwSed = dissolutionFWSed[0]/(10 ** 9)
	except:
		dissolution_fwSed = dissolutionFWSed/(10 ** 9)
	try:
		dissolution_sw = dissolutionSW[0]/(10 ** 9)
	except:
		dissolution_sw = dissolutionSW/(10 ** 9)
	try:
		dissolution_swSed = dissolutionSWSed[0]/(10 ** 9)
	except:
		dissolution_swSed = dissolutionSWSed/(10 ** 9)
	try:
		dissolution_soil1 = dissolutionSoil1[0]/(10 ** 9)
	except:
		dissolution_soil1 = dissolutionSoil1/(10 ** 9)
	try:
		dissolution_soil2 = dissolutionSoil2[0]/(10 ** 9)
	except:
		dissolution_soil2 = dissolutionSoil2/(10 ** 9)
	try:
		dissolution_soil3 = dissolutionSoil3[0]/(10 ** 9)
	except:
		dissolution_soil3 = dissolutionSoil3/(10 ** 9)
	try:
		dissolution_soil4= dissolutionSoil4[0]/(10 ** 9)
	except:
		dissolution_soil4 = dissolutionSoil4/(10 ** 9)

	# 5) partitioning process
	partition_soil2soilw1 = soil2soilwater1/(10 ** 9)
	partition_soil2soilw2 = soil2soilwater2/(10 ** 9)
	partition_soil2soilw3 = soil2soilwater3/(10 ** 9)
	partition_soil2soilw4 = soil2soilwater4/(10 ** 9)
	partition_soilw2soil1 = soilwater2soil1/(10 ** 9)
	partition_soilw2soil2 = soilwater2soil2/(10 ** 9)
	partition_soilw2soil3 = soilwater2soil3/(10 ** 9)
	partition_soilw2soil4 = soilwater2soil4/(10 ** 9)

	# 6) runoff process
	runoff_soil1 = runoffSoil1/(10 ** 9)
	runoff_soil2 = runoffSoil2/(10 ** 9)
	runoff_soil3 = runoffSoil3/(10 ** 9)
	runoff_soil4 = runoffSoil4/(10 ** 9)
	runoff_soil1_dissolved = runoffSoilDis1/(10 ** 9)
	runoff_soil2_dissolved = runoffSoilDis2/(10 ** 9)
	runoff_soil3_dissolved = runoffSoilDis3/(10 ** 9)
	runoff_soil4_dissolved = runoffSoilDis4/(10 ** 9)

	# 7) erosion process
	erosion_soil1 = solidErosionSoil1/(10 ** 9)
	erosion_soil2 = solidErosionSoil2/(10 ** 9)
	erosion_soil3 = solidErosionSoil3/(10 ** 9)
	erosion_soil4 = solidErosionSoil4/(10 ** 9)
	wind_erosion_soil1 = windErosionSoil1/(10 ** 9)
	wind_erosion_soil2 = windErosionSoil2/(10 ** 9)
	wind_erosion_soil3 = windErosionSoil3/(10 ** 9)
	wind_erosion_soil4 = windErosionSoil4/(10 ** 9)

	# 8) infiltra process
	infiltra_soil1 = infiltraSoil1/(10 ** 9)
	infiltra_soil2 = infiltraSoil2/(10 ** 9)
	infiltra_soil3 = infiltraSoil3/(10 ** 9)
	infiltra_soil4 = infiltraSoil4/(10 ** 9)

	leach_soil1 = leachSoil1/(10 ** 9)
	leach_soil2 = leachSoil2/(10 ** 9)
	leach_soil3 = leachSoil3/(10 ** 9)
	leach_soil4 = leachSoil4/(10 ** 9)

	# 9) other process
	burial_fwSed = burialFWSed/(10 ** 9)
	burial_swSed = burialSWSed/(10 ** 9)
	resusp_fwSed = resuspensionFWSed/(10 ** 9)
	resusp_swSed = resuspensionSWSed/(10 ** 9)
	aero_resusp_sSS = aerosolizationSW/(10 ** 9)

	processes = [adv_air_in, adv_fw_out, adv_fwSed_out, adv_sw_out, adv_swSed_out, adv_fw_dissolved, adv_fwSed_dissolved,
				 adv_sw_dissolved, adv_swSed_dissolved, dep_dry_air, dep_dry_air_fw, dep_dry_air_sw, dep_dry_air_soil1,
				 dep_dry_air_soil2, dep_dry_air_soil3, dep_dry_air_soil4, dep_wet_air, dep_wet_air_fw, dep_wet_air_sw,
				 dep_wet_air_soil1, dep_wet_air_soil2, dep_wet_air_soil3, dep_wet_air_soil4, dep_fSS, dep_sSS,
				 heteroagg_air, heteroagg_fw, heteroagg_sw, dissolution_fw, dissolution_fwSed, dissolution_sw,
				 dissolution_swSed, dissolution_soil1, dissolution_soil2, dissolution_soil3, dissolution_soil4,
				 partition_soil2soilw1, partition_soil2soilw2, partition_soil2soilw3, partition_soil2soilw4,
				 partition_soilw2soil1, partition_soilw2soil2, partition_soilw2soil3, partition_soilw2soil4,
				 runoff_soil1, runoff_soil2, runoff_soil3, runoff_soil4, runoff_soil1_dissolved, runoff_soil2_dissolved,
				 runoff_soil3_dissolved, runoff_soil4_dissolved, erosion_soil1, erosion_soil2, erosion_soil3, erosion_soil4,
				 wind_erosion_soil1, wind_erosion_soil2, wind_erosion_soil3, wind_erosion_soil4, infiltra_soil1, infiltra_soil2,
				 infiltra_soil3, infiltra_soil4, leach_soil1, leach_soil2, leach_soil3, leach_soil4,
				 burial_fwSed, burial_swSed, resusp_fwSed, resusp_swSed, aero_resusp_sSS]

	return processes








