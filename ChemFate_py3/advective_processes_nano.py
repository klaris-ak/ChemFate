import numpy as np
import datetime
import math


#################################################################
#
#   CLiCC Nanomaterials F&T Model Developed by Dr. Garner
#	Model original devloped in MATLAB
#   Date: December 22nd, 2016
#   Converted to Python by Jill Farnsworth
#
#################################################################


time1 = datetime.datetime.now()
#print time1, "1"

def aerosolResuspension(windspeed,coastalA,enrichFactor,depth,vol):
    # %   Water to aerosol particle transport
    # %   Enrichment factor values (EF) estimated from Eisenreich 1980, Sievering et al.,
    # %   1980 and Eisenreich 1982
    # %   Transport equation adapted from Eisenreich 1982
    # %   Output is in m3/day

    # % assumes bubbles can only be formed in the top foot (1 m of the marine
    # % environment) 1 m/(depth m)

    # % bubbleRate=bubbleFormation(windspeed); % m/d
    # Twaer=enrichFactor*bubbleFormation(windspeed)*coastalA*(1/depth)/vol #% m/d * m2 * m/m / m3
    factor1 = np.true_divide(1, depth)
    factor2 = np.true_divide(factor1, vol)
    Twaer = enrichFactor * bubbleFormation(windspeed) * coastalA * factor2

    return Twaer


def airAdvection(wind,airA,height,airV,partC,partP):
    # Wind Speed given in m/s, converted to m/d

    # airFlow has the units of m^3/day, m/day * m * m
    # Maadv=(np.multiply(wind,86400)*math.sqrt(airA)*height)*((partC/partP)/airV)
    factor1 = np.true_divide(partC,partP)
    factor2 = np.true_divide(factor1,airV)
    Maadv = (np.multiply(wind, 86400) * math.sqrt(airA) * height) * factor2
    return Maadv


def dryDepAer(partP,airP,dynVisc,radiusParticles,height):
    # Stoke's Law to calculate deposition of aerosols and the associated ENMs
    # gravitation acceleration in m/s^2, hardcoded since this doesn't really vary
    g = 9.8
    # average aerosols or ENMs particle radius in meters
    # dynamic viscosity of air in kg/ms
    # density of aeroosl (or ENMs) and of air
    # Vset=(2*(partP-airP)*g*(radiusParticles**2))/(9*dynVisc)
    Vset = np.true_divide(2 * (partP - airP) * g * (radiusParticles ** 2), 9 * dynVisc)
    # second to day conversion (m/s) to (m/day)
    Vset=Vset*86400
    # rate per day (m/d) * (1/m) * kg aer/m3 air * m3 aer/kg aer * V air
    # Kdep=Vset*(1/(height*0.1))
    Kdep = Vset * np.true_divide(1, height * 0.1)
    return Kdep


def dryDepAir(partP,airP,dynVisc,radiusParticles,height):
    # Stoke's Law to calculate deposition of aerosols and the associated ENMs
    # Gravitation acceleration in m/s^2 (hardcoded since this doens't really vary)
    g = 9.8
    # Average aerosols or ENMs particle radius in meters
    # dynamic viscosity of air in kg/m*s
    # density of aerosols (or ENMs) and of air
    # Vset=(2*(partP-airP)*9.8*(radiusParticles**2))/(9*dynVisc)
    Vset = np.true_divide(2 * (partP - airP) * 9.8 * (radiusParticles ** 2), 9 * dynVisc)
    # second to day conversion (m/s) to (m/day)
    Vset=Vset*86400
    # rate per day (m/d) * (1/m)
    # Kdep=Vset*(1/(height*0.2))
    Kdep = Vset * np.true_divide(1, height * 0.2)
    return Kdep


def wetDepAir(precip,scavengingENM,area,vol):
    # Output is in m^3/day
    # A single rain drop sweeps through 200000 times the golume of area
    # Mackay 2000, page 58
    # Precip is the precip vector based on your climate's daily precipitation in mm/day
    # airH is height
    # ConcA/p is the volF of the nanoparticles in air
    # precip = precip*(1/1000); mm/d to m/d
    wetDepA=np.true_divide(precip,1000)*np.true_divide(area,vol)*scavengingENM
    return wetDepA


def wetDep(precip,airV,scavenging,area):
    # A single rain drop sweeps through 200000 times the volume of area
    # Mackay 2000, page 58
    # Precip is the precip vector based on your climate's daily precipitation, data given in mm/d
    # Area is either the area of the water or the soil
    # aerVf is the volume fraction of aerosols in air, given in the env data

    wetDep=np.true_divide(precip,1000)*np.true_divide(area,airV)*scavenging
    return wetDep


def waterAdv(flow,wV,partC,partP):
    # % water flow is given in m3/s
    # factor1 = np.true_divide(partC,partP)
    # factor2 = np.true_divide(factor1,wV)
    # Mwadv=np.multiply(flow,86400)*factor2
    Mwadv = np.multiply(flow, 86400)/wV
    return Mwadv


# def vertFlow(area,leachingR,vol):
# 	# % Mass transfer of soil water to deep soil transfer
# 	Qvert=np.true_divide(leachingR*area*0.05*24,vol) #% m3/m2-hr * m2 * hr / m3= 1/day
# 	return Qvert


def vertFlow(precip_mm, CN, evap_mm, FC, soilWC, soilV, soilA):
    # % Mass transfer of soil water to deep soil transfer
    # FC = field capacity(m^3/m^3): https://stormwater.pca.state.mn.us/index.php?title=Soil_water_storage_properties
    runoff_mm = 0

    if (precip_mm * 0.0393701) > (0.2 * CN):
        runoff_mm = (((precip_mm * 0.0393701) - 0.2 * CN) ** 2) / (((precip_mm * 0.0393701) + 0.8 * CN))

    # unit: mm/day
    infil_mm = precip_mm - runoff_mm - evap_mm

    if infil_mm <= 0:
        k_infil, Qvert = 0, 0
    else:
        # unit: unitless * m3 + mm/day * m/1000mm * m2 = m3
        soil_water = soilWC * soilV + infil_mm * 0.001 * soilA

        # unit: unitless * m3 = m3
        if soil_water >= FC * soilV:
            # unit: m3/day
            k_infil = soil_water - FC * soilV
        else:
            k_infil = 0

        # unit: m3/day / m3 = 1/day
        Qvert = k_infil / soilV
    return k_infil, Qvert

def horiFlow(k_infil, soilV):
    # unit: m3/day / m3 = 1/day
    Qhori = k_infil / soilV
    return Qhori


def soilwaterPartition(MassS,MassSW,elution,direction):
    # %   If the ratio of mass in sed to mass is sed water is less than 30%,
    # %   then transfer occurs to sed water to balance that ratio
    # %   If the ratio is greater than 30%, then the transfer occurs towards sed
    # %   to balance out the ratio
    # %   dir is the direction of the transport reported, 1(solid to water or 2 (water to solid) (of which one
    # %   will always be 0)
    Mpart1=0
    Mpart2=0
    if (MassS+MassSW) != 0: #prevent divide by zero
        if np.true_divide((MassSW),(MassS+MassSW)) < (1-elution):
            Mpart1=((MassS+MassSW)*(1-elution))-(MassSW)
            Mpart2=0

        if direction==2:
            if np.true_divide((MassSW),(MassS+MassSW))>(1-elution) and (MassSW)>(((MassS+MassSW)*elution)-MassS):
                Mpart2=((MassS+MassSW)*elution)-MassS
                Mpart1=0
            else:
                Mpart2=MassSW
                Mpart1=0

    return [Mpart1, Mpart2]


def sedDeposition(ssP,wP,dynVisc,radiusParticles,depth):
    # %  Stokes law to calculate sedimentation of suspended sediment
    # % g=9.8; % gravitational acceleration in m/s2
    # % average suspended sediment particle radius in m
    # % dynamic viscosity of water in kg/ms
    # Vset=(2*(ssP-wP)*9.8*(radiusParticles**2))/(9*dynVisc)
    Vset = np.true_divide(2*(ssP-wP)*9.8*(radiusParticles**2),9*dynVisc)
    Vset=Vset*86400 #% second to day conversion (m/s) to (m/day)
    Kdep=Vset*np.true_divide(1,depth) #% rate per day (m/d) * (1/m)
    return Kdep


def resuspensionSed(sedA,resuspensionRate,sedV):
    # %   Resuspension rate from sediment to water m3/m2/day
    Msedw=np.true_divide(sedA*resuspensionRate*24,sedV)  #% m3*m2-d * m2 / m3
    return Msedw


def burial(sedA,burialRate,sedV):
    # %   Burial is given a m3/m2-hr
    # %   accounts for only the fraction adsorbed to the sediments

    # burial=burialRate*sedA*24/sedV #% m3/m2-hr * m2 / m3
    burial = np.true_divide(burialRate * sedA * 24, sedV)

    return burial


def runoff(precipI,CN,area, soilV):
    # %   precip is your precipitation dataset
    # %   CN is the NRCS curve number, ranges from 30 -> 100 with 30 being low
    # %   runoff and 100 being high runoff and high imprevious surface area
    # %   area is the area of the compartment that the water falls onto
    # % calculations need to be done in inches,
    # % precip=precip*(1/25.4)*(1/24)
    # % conversion factor from mm/d to inches
    Q=0
    if (precipI*(1/25.4))>(0.2*CN):
        Q=(((precipI*(1/25.4))-0.2*CN)**2)/(((precipI*(1/25.4))+0.8*CN))

    # %http://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1082989.pdf
    # % 1 in/hr = 1 cfs/acre
    # % ft3/s-acre * M3/ft3 * s/d * acre/m2 * m2 /m3 = 1/day
    Qw=(Q)*(1/35.312)*(86400)*(1/4046.86)*area/soilV
    return Qw


def erosion(precipI,Kfact,lenslope,cropManageFactor,supportFactor,area,soilP):
    # % Rainfall erosivity factor
    # % http://www.sciencedirect.com.proxy.library.ucsb.edu:2048/science/article/pii/S004896971500011X

    er=0.29*(1-0.72*math.exp(-0.05*precipI)) # % precip in mm/d
    # % this implicitly contains day because we don't use 30 minute storm
    # % intensity, we use 1 day storm total
    EI=er*precipI*(precipI)*np.true_divide(1,24) # % (MJ mm ha-1 hr-1 d-1)
    #% er MJ/mm-ha-day * V in mm * intensity in mm/d (assumed same as 30-min
    #% intensity) * 1 day/24 hr

    R_unit=EI*(1/2.471)*(1/25.4)*np.true_divide(7375.62,2000) # % converts R to SI units
    #% MJ*mm/ha*hr*day * ha/2.471 acre * in/25.4 mm * 7375.62 hundred ft-lbs/MJ *
    #% tons/2000lbs = hundred ft-ton-in/acre-hr-day

    # % soil erodibility factor, taken from statsgo dataset
    # %K=Kfact; % % tons hr/100s ft ton inches

    # % slope lenght gradient factor
    # % length will always be 1000 m since we are looking at large distances/
    # % long term
    # % slope can come from the soil data
    # % factor calculated in excel then read in here
    # % do we have to account for actual length relative to the 1000 m???
    # %LS=lenslope;

    # % crop vegetation management factor
    # % 0.005 for continuous forests to 0.5 for crops, maybe depends on land
    # % cover
    # %C=cropManageFactor;

    # % support practice factor
    # % values range from 0 to 1.  1 for straight row farming.  The lower the
    # % value the better technology implemented to prevent erosion
    # %P=supportFactor;

    # % RUSLE
    # % units are tons/acre/year
    # % convert to kg/m2/year to kg/year over total area
    # converted to ug/yr
    # converted to m3/yr using density
    soilLoss=np.true_divide(((R_unit*Kfact*lenslope*cropManageFactor*supportFactor)*((907.185*(10**9))/4046.86)*area),soilP);

    # % fractional amount of soil lost by day (precipI - daily precip (mm/d) relative to
    # % total precip P (mm/d))
    # soilLoss=(Atotal*(precipI/(np.sum(precip)/years)))/soilP  # % kg/d to m3/d
    #if i==0:
        # precip i cancels out and this is separate because can't divide by 0
     #   soilLoss=np.true_divide(Atotal/365.24, soilP)
    #elif i<365:
    # if time is less than 1 year, first year annual average or partial year average
    #	factor0 = np.true_divide(365.24,i)
    #	factor1 = np.sum(precip[0:i+1])*factor0
    #	factor2 = np.true_divide(precipI, factor1)
    #	soilLoss = np.true_divide(Atotal*factor2, soilP)
    #elif i+183<=time:
    # else take the precip of the year around that day in time
    #	factor2 = np.true_divide(precipI, np.sum(precip[i-183:i+184]))
    #	soilLoss = np.true_divide(Atotal*factor2, soilP)
    #else:
    # up until the final year then shift the average back to still include 1 full year
    #	t2oneyear = 365-(time-(i-183))
    #	factor2 = np.true_divide(precipI, np.sum(precip[i-183-t2oneyear:time]))
    #	soilLoss = np.true_divide(Atotal*factor2, soilP)

    return soilLoss


def windErosion(wind,precipI,roughness,Kconstant,airP,soilA,A,TSV,TSVmin,z_wind,percWind,windConstant,percUncovered,percSuspended,soilP,soilV):
    # %   TSV,TSVmin,z_wind,percWind,windConstant,percUncovered,percSuspended,soilP,soilV)
    # %   annualFlux per land use type is given in kg/m2 year
    # %   We assume that the starting soil moisture is currently 0, though we may
    # %   want to change this.
    # %   We assume that if there was no rain that day erosion can occur.
    # %   We use the saltation equation, given in Kelly et al. 2004, and the
    # %   vertical flux conversion to estimate total transport of soil between
    # %   soil and aerosols
    # %   Wind profile and shear rates are taken from Gillette 1978
    # %   z_wind is height at which wind measurements were taken (in US, typically 1.5 m)
    # %   roughness value needs to be in m
    # %   K constant to convert from horizontal to vertical in 1/m
    # %   airP density of air, kg/m3
    # % Threshold velocity for erosion (m/s)
    # % minimum wind speed needed to cause erosion
    # % if there was precipitation recently, then the minimum threshold is
    # % increased to 30 m/s windspeed, cause not much wind erosion when soil
    # % is saturated

    if precipI>10:
        TSV=TSVmin #% m/s

    # % Uz=wind; % m/s
    # % k=0.41; % von Karmen constant, unitless
    # % u* is the wind shear velocity in m/s for each day based on the
    # % existing wind speed
    ustar=(wind*0.41)/(math.log(np.true_divide(z_wind,roughness))) #% m/s / m/m = m/s

    # % g=9.81; % acceleration due to gravity in m/s2
    # % calculate horizontal Qtot and vertical Fa fluxes
    Fa=0
    if ustar>TSV:
        # % A * (airP/9.81)*(ustar*((ustar^2)-(TSV^2)))) if Qtot
        Fa=(A*(airP/9.81)*(ustar*((ustar**2)-(TSV**2))))*Kconstant
        # % K is the conversion constant between horizontal and vertical fluxes in 1/m
        # % kg/m-s * m-1 = kg/m2 s

    # % conversion to proper units
    # % if wind constantly throughout the day Fa*86400
    # % but it doesn't so Fa*86400*0.3 because wind dies down at night (~10
    # % hours) and is not consistant throughout the day (10% of the time)
    Fa=Fa*86400*percWind*windConstant #% convert to (kg/m2 day)
    erosionAer= np.true_divide(Fa*soilA*percUncovered*percSuspended,soilP*soilV) #% 1/day

    return erosionAer


def odeHetagg(k,ssC):
    # %   Assuming heteroaggregation is a pseudo first order rate constant where
    # %   rate = k'[C]^1 because k'=k*[ssC] the equation to solver for the
    # %   concentration at time (t) is simply C=k*initC*ssC
    # %   We do not need to use a proper differential equation solver because you
    # %   just need the solution to update with each time step.
    # %   Since we are solving for concentration in kg/m3, we need the units of k
    # %   to be in terms of 1/day (technically m3/kg per day) because we are solving over the change
    # %   of a day
    Mwss = k*ssC #% m3/kg-d * kg/m3
    return Mwss


def odeDissolution(a, b, MassENM, MassDiss, k, V):
    # % a is the predicted exponential equation a value
    # % b is the predicted exponential equation b value
    # % MassDiss is the dissolved fraction kg
    # % MassENM is the mass of the ENM
    # % k is the dissolution rate constant
    # % V is the volume of the compartment
    # % this assumes first order rate constant for dissolution
    # % Exponential fit equation for dissolution data
    eqPerc = 99.9999
    if np.true_divide((MassENM + MassDiss),
                      V) > 10 ** 2:  # % typically equilibrium dissolution is 100% at anything less than 10^2 pg/m3
        # % one place to improve efficiency would sbe to calculate out this
        # % concentration for each metal and pH
        eqPerc = a * math.exp(b * (np.true_divide(MassENM + MassDiss, V)))
        if eqPerc >= 100:
            eqPerc = 99.9999
        elif eqPerc < 0.0000001:
            eqPerc = 0

    Mdisw = k * MassENM  # %  1/d * kg
    if k == 0 or MassENM == 0:
        Mdisw = 0  # prevent divide by zero
    elif (np.true_divide(MassDiss, (MassENM + MassDiss)) * 100) >= eqPerc:
        Mdisw = 0
    # (Mdisw+MassDiss)  is new dissolved mass kg + kg
    # (MassENM+MassDiss) total mass kg + kg
    elif ((Mdisw + MassDiss) / (MassENM + MassDiss)) * 100 >= eqPerc:
        Mdisw = np.true_divide((eqPerc - 0.0000001), 100) * (MassENM + MassDiss)  # % kg
        if Mdisw < 0:
            Mdisw = 0

    if Mdisw > MassDiss:
        lossENM = Mdisw - MassDiss
    else:
        lossENM = 0

    # elif (np.true_divide((np.true_divide(MassDiss,V)),np.true_divide((MassENM+MassDiss),V)))*100>=eqPerc:  #% the additional dissolved concentration cannot exceed the equilibrium concentration
    #	Mdisw=np.true_divide((eqPerc-0.00001),100)*MassENM*V #% kg
    # if eqPerc is very low, Mdisw could go negative, so reset to zero
    #	if Mdisw<0:
    #	    Mdisw=0
    #	#% if the ratio of the dissolved concentration to the total is greater than the equilibrium dissolve concentrations then  no more can dissolve
    # if (np.true_divide(MassDiss,(MassENM+MassDiss))*100)>=eqPerc:
    #	Mdisw=0
    return [Mdisw, lossENM]


def bubbleFormation(windspeed):
    # %   Droplet volume flux (bubbleRate) taken from Kerman 1986, dependent on
    # %   windspeed (1*10^-9 cm/s at windspeeds of 6 m/s, and (9*10^-4 cm/s at windspeeds
    # %   greater than 12 m/s)

    if windspeed < 6:  # % m/s
        bubbleRate = 0
    elif windspeed > 12:  # % m/s
        # bubbleRate=(9*10**-4)*86400/100 #% convert from cm/s to m/day,
        bubbleRate = (9 * 10 ** -4) * np.true_divide(86400, 100)
    # %rate from Wu et al. 1984
    else:
        # bubbleRate=(1*10**-9)*86400/100 #% convert to m/day,
        bubbleRate = (1 * 10 ** -9) * np.true_divide(86400, 100)
    # %rate from blanchard 1963

    return bubbleRate


def lsFactor(slope):
    if slope<=0.2:
        lenslope=0.06
    elif slope<=0.5:
        lenslope=0.1
    elif slope<=1:
        lenslope=0.2
    elif slope<=2:
        lenslope=0.47
    elif slope<=3:
        lenslope=0.8
    elif slope<=4:
        lenslope=1.19
    elif slope<=5:
        lenslope=1.63
    elif slope<=6:
        lenslope=2.11
    elif slope<=8:
        lenslope=3.15
    elif slope<=10:
        lenslope=4.56
    elif slope<=12:
        lenslope=6.28
    elif slope<=14:
        lenslope=8.11
    elif slope<=16:
        lenslope=10.02
    elif slope<=20:
        lenslope=13.99
    elif slope<=25:
        lenslope=19.13
    elif slope<=30:
        lenslope=24.31
    elif slope<=40:
        lenslope=34.48
    elif slope<=50:
        lenslope=44.02
    elif slope<=60:
        lenslope=52.7
    return lenslope
