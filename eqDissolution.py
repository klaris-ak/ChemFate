import datetime
import scipy.io as sio

#################################################################
#
#   CLiCC Nanomaterials F&T Model Developed by Dr. Garner
#	Model original devloped in MATLAB
#   Date: December 22nd, 2016
#   Converted to Python by Jill Farnsworth
#
#################################################################


def eqDissolution(ENM,FWpH,SWpH,GWpH1,GWpH2,GWpH3,GWpH4,presence):
	# %   Inputs include which ENM you are modeling, the water pH for each of the
	# %   compartments that dissolution can occur, and the presence indicator for
	# %   each compartment

	# %   This calculates the equilibrium dissolution values based on the pH
	# %   accross a range of ENM concentrations

	# substitution for following lines of Matlab:
	# 	FWlow=find(vec<=FWpH,1,'last');
	# 	FWhigh=find(vec>=FWpH,1,'first');
	# 	FWlower=vec(FWlow);
	# 	FWhigher=vec(FWhigh);
	def sub(a, func, pos):
		res = [i for (i, val) in enumerate(a) if func(val)]
		if res == []:
			return res
		else:
			if pos == 'first':
				ind = res[0]
				return a[ind]
			elif pos == 'last':
				ind = res[-1]
				return a[ind]
			else:
				print("invalid position\n")

	# load correct mat file
	if ENM == 'Ag':
		contents = sio.loadmat('./Ag_eq_dis.mat')
	elif ENM == 'CuO':
		contents = sio.loadmat('./CuO_eq_dis.mat')
	elif ENM == 'nZVI':
		contents = sio.loadmat('./nZVI_eq_dis.mat')
	elif ENM == 'TiO2':
		contents = sio.loadmat('./TiO2_eq_dis.mat')
	elif ENM == 'ZnO':
		contents = sio.loadmat('./ZnO_eq_dis.mat')
	elif ENM == 'CeO2':
		contents = sio.loadmat('./CeO2_eq_dis.mat')
	elif ENM == 'SiO2':
		contents = sio.loadmat('./SiO2_eq_dis.mat')

	#%% Freshwater

	if presence['fw']==0: #% if there is no freshwater compartment
		percfitaFW=0 
		percfitbFW=0
	else:
		# % Finding the closest FW pH
		vec=[7.5,8,8.33,8.5,9]
		FWlower = sub(vec, lambda x: x<=FWpH, 'last')
		FWhigher = sub(vec, lambda x: x>=FWpH, 'first')
		if FWlower == []:
			FWpH=FWhigher
		elif FWhigher == []:
			FWpH=FWlower
		elif abs(FWlower-FWpH)<abs(FWhigher-FWpH):
			FWpH=FWlower
		else:
			FWpH=FWhigher

	if ENM == 'Ag':
		if FWpH==7.5:
			percfitaFW=contents['Ag75FWpercfita'][0][0]
			percfitbFW=contents['Ag75FWpercfitb'][0][0]
		elif FWpH==8:
			percfitaFW=contents['Ag8FWpercfita'][0][0]
			percfitbFW=contents['Ag8FWpercfitb'][0][0]
		elif FWpH==8.33:
			percfitaFW=contents['Ag833FWpercfita'][0][0]
			percfitbFW=contents['Ag833FWpercfitb'][0][0]
		elif FWpH==8.5:
			percfitaFW=contents['Ag85FWpercfita'][0][0]
			percfitbFW=contents['Ag85FWpercfitb'][0][0]
		elif FWpH==9:
			percfitaFW=contents['Ag9FWpercfita'][0][0]
			percfitbFW=contents['Ag9FWpercfitb'][0][0]
	elif ENM == 'CuO':
		if FWpH==7.5:
			percfitaFW=contents['CuO75FWpercfita'][0][0]
			percfitbFW=contents['CuO75FWpercfitb'][0][0]
		elif FWpH==8:
			percfitaFW=contents['CuO8FWpercfita'][0][0]
			percfitbFW=contents['CuO8FWpercfitb'][0][0]
		elif FWpH==8.33:
			percfitaFW=contents['CuO833FWpercfita'][0][0]
			percfitbFW=contents['CuO833FWpercfitb'][0][0]
		elif FWpH==8.5:
			percfitaFW=contents['CuO85FWpercfita'][0][0]
			percfitbFW=contents['CuO85FWpercfitb'][0][0]
		elif FWpH==9:
			percfitaFW=contents['CuO9FWpercfita'][0][0]
			percfitbFW=contents['CuO9FWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if FWpH==7.5:
			percfitaFW=contents['nZVI75FWpercfita'][0][0]
			percfitbFW=contents['nZVI75FWpercfitb'][0][0]
		elif FWpH==8:
			percfitaFW=contents['nZVI8FWpercfita'][0][0]
			percfitbFW=contents['nZVI8FWpercfitb'][0][0]
		elif FWpH==8.33:
			percfitaFW=contents['nZVI833FWpercfita'][0][0]
			percfitbFW=contents['nZVI833FWpercfitb'][0][0]
		elif FWpH==8.5:
			percfitaFW=contents['nZVI85FWpercfita'][0][0]
			percfitbFW=contents['nZVI85FWpercfitb'][0][0]
		elif FWpH==9:
			percfitaFW=contents['nZVI9FWpercfita'][0][0]
			percfitbFW=contents['nZVI9FWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if FWpH==7.5:
			percfitaFW=contents['ZnO75FWpercfita'][0][0]
			percfitbFW=contents['ZnO75FWpercfitb'][0][0]
		elif FWpH==8:
			percfitaFW=contents['ZnO8FWpercfita'][0][0]
			percfitbFW=contents['ZnO8FWpercfitb'][0][0]
		elif FWpH==8.33:
			percfitaFW=contents['ZnO833FWpercfita'][0][0]
			percfitbFW=contents['ZnO833FWpercfitb'][0][0]
		elif FWpH==8.5:
			percfitaFW=contents['ZnO85FWpercfita'][0][0]
			percfitbFW=contents['ZnO85FWpercfitb'][0][0]
		elif FWpH==9:
			percfitaFW=contents['ZnO9FWpercfita'][0][0]
			percfitbFW=contents['ZnO9FWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaFW=contents['TiO2FWpercfita'][0][0]
		percfitbFW=contents['TiO2FWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaFW=contents['CeO2FWpercfita'][0][0]
		percfitbFW=contents['CeO2FWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaFW=contents['SiO2FWpercfita'][0][0]
		percfitbFW=contents['SiO2FWpercfitb'][0][0]


	#%% Marine

	if presence['sw']==0: #% if there is no freshwater compartment
		percfitaSW=0 
		percfitbSW=0
	else:
		# % Finding the closest SW pH
		vec=[7,7.5,8.05,8.5,9]
		SWlower = sub(vec, lambda x: x<=SWpH, 'last')
		SWhigher = sub(vec, lambda x: x>=SWpH, 'first')
		if SWlower == []:
			SWpH=SWhigher
		elif SWhigher == []:
			SWpH=SWlower
		elif abs(SWlower-SWpH)<abs(SWhigher-SWpH):
			SWpH=SWlower
		else:
			SWpH=SWhigher

	if ENM == 'Ag':
		if SWpH==7:
			percfitaSW=contents['Ag7SWpercfita'][0][0]
			percfitbSW=contents['Ag7SWpercfitb'][0][0]
		elif SWpH==7.5:
			percfitaSW=contents['Ag75SWpercfita'][0][0]
			percfitbSW=contents['Ag75SWpercfitb'][0][0]
		elif SWpH==8.05:
			percfitaSW=contents['Ag805SWpercfita'][0][0]
			percfitbSW=contents['Ag805SWpercfitb'][0][0]
		elif SWpH==8.5:
			percfitaSW=contents['Ag85SWpercfita'][0][0]
			percfitbSW=contents['Ag85SWpercfitb'][0][0]
		elif SWpH==9:
			percfitaSW=contents['Ag9SWpercfita'][0][0]
			percfitbSW=contents['Ag9SWpercfitb'][0][0]
	elif ENM == 'CuO':
		if SWpH==7:
			percfitaSW=contents['CuO7SWpercfita'][0][0]
			percfitbSW=contents['CuO7SWpercfitb'][0][0]
		elif SWpH==7.5:
			percfitaSW=contents['CuO75SWpercfita'][0][0]
			percfitbSW=contents['CuO75SWpercfitb'][0][0]
		elif SWpH==8.05:
			percfitaSW=contents['CuO805SWpercfita'][0][0]
			percfitbSW=contents['CuO805SWpercfitb'][0][0]
		elif SWpH==8.5:
			percfitaSW=contents['CuO85SWpercfita'][0][0]
			percfitbSW=contents['CuO85SWpercfitb'][0][0]
		elif SWpH==9:
			percfitaSW=contents['CuO9SWpercfita'][0][0]
			percfitbSW=contents['CuO9SWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if SWpH==7:
			percfitaSW=contents['nZVI7SWpercfita'][0][0]
			percfitbSW=contents['nZVI7SWpercfitb'][0][0]
		elif SWpH==7.5:
			percfitaSW=contents['nZVI75SWpercfita'][0][0]
			percfitbSW=contents['nZVI75SWpercfitb'][0][0]
		elif SWpH==8.05:
			percfitaSW=contents['nZVI805SWpercfita'][0][0]
			percfitbSW=contents['nZVI805SWpercfitb'][0][0]
		elif SWpH==8.5:
			percfitaSW=contents['nZVI85SWpercfita'][0][0]
			percfitbSW=contents['nZVI85SWpercfitb'][0][0]
		elif SWpH==9:
			percfitaSW=contents['nZVI9SWpercfita'][0][0]
			percfitbSW=contents['nZVI9SWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if SWpH==7:
			percfitaSW=contents['ZnO7SWpercfita'][0][0]
			percfitbSW=contents['ZnO7SWpercfitb'][0][0]
		elif SWpH==7.5:
			percfitaSW=contents['ZnO75SWpercfita'][0][0]
			percfitbSW=contents['ZnO75SWpercfitb'][0][0]
		elif SWpH==8.05:
			percfitaSW=contents['ZnO805SWpercfita'][0][0]
			percfitbSW=contents['ZnO805SWpercfitb'][0][0]
		elif SWpH==8.5:
			percfitaSW=contents['ZnO85SWpercfita'][0][0]
			percfitbSW=contents['ZnO85SWpercfitb'][0][0]
		elif SWpH==9:
			percfitaSW=contents['ZnO9SWpercfita'][0][0]
			percfitbSW=contents['ZnO9SWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaSW=contents['TiO2SWpercfita'][0][0]
		percfitbSW=contents['TiO2SWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaSW=contents['CeO2SWpercfita'][0][0]
		percfitbSW=contents['CeO2SWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaSW=contents['SiO2SWpercfita'][0][0]
		percfitbSW=contents['SiO2SWpercfitb'][0][0]


	#%% Groundwater 1

	if presence['soilW1']==0: #% if there is no freshwater compartment
		percfitaGW1=0 
		percfitbGW1=0
	else:
		# % Finding the closest GW pH
		vec=[6,6.5,7,7.5,8,8.5]
		GW1lower = sub(vec, lambda x: x<=GWpH1, 'last')
		GW1higher = sub(vec, lambda x: x>=GWpH1, 'first')
		if GW1lower == []:
			GWpH1=GW1higher
		elif GW1higher == []:
			GWpH1=GW1lower
		elif abs(GW1lower-GWpH1)<abs(GW1higher-GWpH1):
			GWpH1=GW1lower
		else:
			GWpH1=GW1higher

	if ENM == 'Ag':
		if GWpH1==6:
			percfitaGW1=contents['Ag6GWpercfita'][0][0]
			percfitbGW1=contents['Ag6GWpercfitb'][0][0]
		elif GWpH1==6.5:
			percfitaGW1=contents['Ag65GWpercfita'][0][0]
			percfitbGW1=contents['Ag65GWpercfitb'][0][0]
		elif GWpH1==7:
			percfitaGW1=contents['Ag7GWpercfita'][0][0]
			percfitbGW1=contents['Ag7GWpercfitb'][0][0]
		elif GWpH1==7.5:
			percfitaGW1=contents['Ag75GWpercfita'][0][0]
			percfitbGW1=contents['Ag75GWpercfitb'][0][0]
		elif GWpH1==8:
			percfitaGW1=contents['Ag8GWpercfita'][0][0]
			percfitbGW1=contents['Ag8GWpercfitb'][0][0]
	elif ENM == 'CuO':
		if GWpH1==6:
			percfitaGW1=contents['CuO6GWpercfita'][0][0]
			percfitbGW1=contents['CuO6GWpercfitb'][0][0]
		elif GWpH1==6.5:
			percfitaGW1=contents['CuO65GWpercfita'][0][0]
			percfitbGW1=contents['CuO65GWpercfitb'][0][0]
		elif GWpH1==7:
			percfitaGW1=contents['CuO7GWpercfita'][0][0]
			percfitbGW1=contents['CuO7GWpercfitb'][0][0]
		elif GWpH1==7.5:
			percfitaGW1=contents['CuO75GWpercfita'][0][0]
			percfitbGW1=contents['CuO75GWpercfitb'][0][0]
		elif GWpH1==8:
			percfitaGW1=contents['CuO8GWpercfita'][0][0]
			percfitbGW1=contents['CuO8GWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if GWpH1==6:
			percfitaGW1=contents['nZVI6GWpercfita'][0][0]
			percfitbGW1=contents['nZVI6GWpercfitb'][0][0]
		elif GWpH1==6.5:
			percfitaGW1=contents['nZVI65GWpercfita'][0][0]
			percfitbGW1=contents['nZVI65GWpercfitb'][0][0]
		elif GWpH1==7:
			percfitaGW1=contents['nZVI7GWpercfita'][0][0]
			percfitbGW1=contents['nZVI7GWpercfitb'][0][0]
		elif GWpH1==7.5:
			percfitaGW1=contents['nZVI75GWpercfita'][0][0]
			percfitbGW1=contents['nZVI75GWpercfitb'][0][0]
		elif GWpH1==8:
			percfitaGW1=contents['nZVI8GWpercfita'][0][0]
			percfitbGW1=contents['nZVI8GWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if GWpH1==6:
			percfitaGW1=contents['ZnO6GWpercfita'][0][0]
			percfitbGW1=contents['ZnO6GWpercfitb'][0][0]
		elif GWpH1==6.5:
			percfitaGW1=contents['ZnO65GWpercfita'][0][0]
			percfitbGW1=contents['ZnO65GWpercfitb'][0][0]
		elif GWpH1==7:
			percfitaGW1=contents['ZnO7GWpercfita'][0][0]
			percfitbGW1=contents['ZnO7GWpercfitb'][0][0]
		elif GWpH1==7.5:
			percfitaGW1=contents['ZnO75GWpercfita'][0][0]
			percfitbGW1=contents['ZnO75GWpercfitb'][0][0]
		elif GWpH1==8:
			percfitaGW1=contents['ZnO8GWpercfita'][0][0]
			percfitbGW1=contents['ZnO8GWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaGW1=contents['TiO2GWpercfita'][0][0]
		percfitbGW1=contents['TiO2GWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaGW1=contents['CeO2GWpercfita'][0][0]
		percfitbGW1=contents['CeO2GWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaGW1=contents['SiO2GWpercfita'][0][0]
		percfitbGW1=contents['SiO2GWpercfitb'][0][0]


	#%% Groundwater 2

	if presence['soilW2']==0: #% if there is no freshwater compartment
		percfitaGW2=0 
		percfitbGW2=0
	else:
		# % Finding the closest GW pH
		vec=[6,6.5,7,7.5,8,8.5]
		GW2lower = sub(vec, lambda x: x<=GWpH2, 'last')
		GW2higher = sub(vec, lambda x: x>=GWpH2, 'first')
		if GW2lower == []:
			GWpH2=GW2higher
		elif GW2higher == []:
			GWpH2=GW2lower
		elif abs(GW2lower-GWpH2)<abs(GW2higher-GWpH2):
			GWpH2=GW2lower
		else:
			GWpH2=GW2higher

	if ENM == 'Ag':
		if GWpH2==6:
			percfitaGW2=contents['Ag6GWpercfita'][0][0]
			percfitbGW2=contents['Ag6GWpercfitb'][0][0]
		elif GWpH2==6.5:
			percfitaGW2=contents['Ag65GWpercfita'][0][0]
			percfitbGW2=contents['Ag65GWpercfitb'][0][0]
		elif GWpH2==7:
			percfitaGW2=contents['Ag7GWpercfita'][0][0]
			percfitbGW2=contents['Ag7GWpercfitb'][0][0]
		elif GWpH2==7.5:
			percfitaGW2=contents['Ag75GWpercfita'][0][0]
			percfitbGW2=contents['Ag75GWpercfitb'][0][0]
		elif GWpH2==8:
			percfitaGW2=contents['Ag8GWpercfita'][0][0]
			percfitbGW2=contents['Ag8GWpercfitb'][0][0]
	elif ENM == 'CuO':
		if GWpH2==6:
			percfitaGW2=contents['CuO6GWpercfita'][0][0]
			percfitbGW2=contents['CuO6GWpercfitb'][0][0]
		elif GWpH2==6.5:
			percfitaGW2=contents['CuO65GWpercfita'][0][0]
			percfitbGW2=contents['CuO65GWpercfitb'][0][0]
		elif GWpH2==7:
			percfitaGW2=contents['CuO7GWpercfita'][0][0]
			percfitbGW2=contents['CuO7GWpercfitb'][0][0]
		elif GWpH2==7.5:
			percfitaGW2=contents['CuO75GWpercfita'][0][0]
			percfitbGW2=contents['CuO75GWpercfitb'][0][0]
		elif GWpH2==8:
			percfitaGW2=contents['CuO8GWpercfita'][0][0]
			percfitbGW2=contents['CuO8GWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if GWpH2==6:
			percfitaGW2=contents['nZVI6GWpercfita'][0][0]
			percfitbGW2=contents['nZVI6GWpercfitb'][0][0]
		elif GWpH2==6.5:
			percfitaGW2=contents['nZVI65GWpercfita'][0][0]
			percfitbGW2=contents['nZVI65GWpercfitb'][0][0]
		elif GWpH2==7:
			percfitaGW2=contents['nZVI7GWpercfita'][0][0]
			percfitbGW2=contents['nZVI7GWpercfitb'][0][0]
		elif GWpH2==7.5:
			percfitaGW2=contents['nZVI75GWpercfita'][0][0]
			percfitbGW2=contents['nZVI75GWpercfitb'][0][0]
		elif GWpH2==8:
			percfitaGW2=contents['nZVI8GWpercfita'][0][0]
			percfitbGW2=contents['nZVI8GWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if GWpH2==6:
			percfitaGW2=contents['ZnO6GWpercfita'][0][0]
			percfitbGW2=contents['ZnO6GWpercfitb'][0][0]
		elif GWpH2==6.5:
			percfitaGW2=contents['ZnO65GWpercfita'][0][0]
			percfitbGW2=contents['ZnO65GWpercfitb'][0][0]
		elif GWpH2==7:
			percfitaGW2=contents['ZnO7GWpercfita'][0][0]
			percfitbGW2=contents['ZnO7GWpercfitb'][0][0]
		elif GWpH2==7.5:
			percfitaGW2=contents['ZnO75GWpercfita'][0][0]
			percfitbGW2=contents['ZnO75GWpercfitb'][0][0]
		elif GWpH2==8:
			percfitaGW2=contents['ZnO8GWpercfita'][0][0]
			percfitbGW2=contents['ZnO8GWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaGW2=contents['TiO2GWpercfita'][0][0]
		percfitbGW2=contents['TiO2GWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaGW2=contents['CeO2GWpercfita'][0][0]
		percfitbGW2=contents['CeO2GWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaGW2=contents['SiO2GWpercfita'][0][0]
		percfitbGW2=contents['SiO2GWpercfitb'][0][0]


	#%% Groundwater 3

	if presence['soilW3']==0: #% if there is no freshwater compartment
		percfitaGW3=0 
		percfitbGW3=0
	else:
		# % Finding the closest GW pH
		vec=[6,6.5,7,7.5,8,8.5]
		GW3lower = sub(vec, lambda x: x<=GWpH3, 'last')
		GW3higher = sub(vec, lambda x: x>=GWpH3, 'first')
		if GW3lower == []:
			GWpH3=GW3higher
		elif GW3higher == []:
			GWpH3=GW3lower
		elif abs(GW3lower-GWpH3)<abs(GW3higher-GWpH3):
			GWpH3=GW3lower
		else:
			GWpH3=GW3higher

	if ENM == 'Ag':
		if GWpH3==6:
			percfitaGW3=contents['Ag6GWpercfita'][0][0]
			percfitbGW3=contents['Ag6GWpercfitb'][0][0]
		elif GWpH3==6.5:
			percfitaGW3=contents['Ag65GWpercfita'][0][0]
			percfitbGW3=contents['Ag65GWpercfitb'][0][0]
		elif GWpH3==7:
			percfitaGW3=contents['Ag7GWpercfita'][0][0]
			percfitbGW3=contents['Ag7GWpercfitb'][0][0]
		elif GWpH3==7.5:
			percfitaGW3=contents['Ag75GWpercfita'][0][0]
			percfitbGW3=contents['Ag75GWpercfitb'][0][0]
		elif GWpH3==8:
			percfitaGW3=contents['Ag8GWpercfita'][0][0]
			percfitbGW3=contents['Ag8GWpercfitb'][0][0]
	elif ENM == 'CuO':
		if GWpH3==6:
			percfitaGW3=contents['CuO6GWpercfita'][0][0]
			percfitbGW3=contents['CuO6GWpercfitb'][0][0]
		elif GWpH3==6.5:
			percfitaGW3=contents['CuO65GWpercfita'][0][0]
			percfitbGW3=contents['CuO65GWpercfitb'][0][0]
		elif GWpH3==7:
			percfitaGW3=contents['CuO7GWpercfita'][0][0]
			percfitbGW3=contents['CuO7GWpercfitb'][0][0]
		elif GWpH3==7.5:
			percfitaGW3=contents['CuO75GWpercfita'][0][0]
			percfitbGW3=contents['CuO75GWpercfitb'][0][0]
		elif GWpH3==8:
			percfitaGW3=contents['CuO8GWpercfita'][0][0]
			percfitbGW3=contents['CuO8GWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if GWpH3==6:
			percfitaGW3=contents['nZVI6GWpercfita'][0][0]
			percfitbGW3=contents['nZVI6GWpercfitb'][0][0]
		elif GWpH3==6.5:
			percfitaGW3=contents['nZVI65GWpercfita'][0][0]
			percfitbGW3=contents['nZVI65GWpercfitb'][0][0]
		elif GWpH3==7:
			percfitaGW3=contents['nZVI7GWpercfita'][0][0]
			percfitbGW3=contents['nZVI7GWpercfitb'][0][0]
		elif GWpH3==7.5:
			percfitaGW3=contents['nZVI75GWpercfita'][0][0]
			percfitbGW3=contents['nZVI75GWpercfitb'][0][0]
		elif GWpH3==8:
			percfitaGW3=contents['nZVI8GWpercfita'][0][0]
			percfitbGW3=contents['nZVI8GWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if GWpH3==6:
			percfitaGW3=contents['ZnO6GWpercfita'][0][0]
			percfitbGW3=contents['ZnO6GWpercfitb'][0][0]
		elif GWpH3==6.5:
			percfitaGW3=contents['ZnO65GWpercfita'][0][0]
			percfitbGW3=contents['ZnO65GWpercfitb'][0][0]
		elif GWpH3==7:
			percfitaGW3=contents['ZnO7GWpercfita'][0][0]
			percfitbGW3=contents['ZnO7GWpercfitb'][0][0]
		elif GWpH3==7.5:
			percfitaGW3=contents['ZnO75GWpercfita'][0][0]
			percfitbGW3=contents['ZnO75GWpercfitb'][0][0]
		elif GWpH3==8:
			percfitaGW3=contents['ZnO8GWpercfita'][0][0]
			percfitbGW3=contents['ZnO8GWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaGW3=contents['TiO2GWpercfita'][0][0]
		percfitbGW3=contents['TiO2GWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaGW3=contents['CeO2GWpercfita'][0][0]
		percfitbGW3=contents['CeO2GWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaGW3=contents['SiO2GWpercfita'][0][0]
		percfitbGW3=contents['SiO2GWpercfitb'][0][0]


	#%% Groundwater 4

	if presence['soilW4']==0: #% if there is no freshwater compartment
		percfitaGW4=0 
		percfitbGW4=0
	else:
		# % Finding the closest GW pH
		vec=[6,6.5,7,7.5,8,8.5]
		GW4lower = sub(vec, lambda x: x<=GWpH4, 'last')
		GW4higher = sub(vec, lambda x: x>=GWpH4, 'first')
		if GW4lower == []:
			GWpH4=GW4higher
		elif GW4higher == []:
			GWpH4=GW4lower
		elif abs(GW4lower-GWpH4)<abs(GW4higher-GWpH4):
			GWpH4=GW4lower
		else:
			GWpH4=GW4higher

	if ENM == 'Ag':
		if GWpH4==6:
			percfitaGW4=contents['Ag6GWpercfita'][0][0]
			percfitbGW4=contents['Ag6GWpercfitb'][0][0]
		elif GWpH4==6.5:
			percfitaGW4=contents['Ag65GWpercfita'][0][0]
			percfitbGW4=contents['Ag65GWpercfitb'][0][0]
		elif GWpH4==7:
			percfitaGW4=contents['Ag7GWpercfita'][0][0]
			percfitbGW4=contents['Ag7GWpercfitb'][0][0]
		elif GWpH4==7.5:
			percfitaGW4=contents['Ag75GWpercfita'][0][0]
			percfitbGW4=contents['Ag75GWpercfitb'][0][0]
		elif GWpH4==8:
			percfitaGW4=contents['Ag8GWpercfita'][0][0]
			percfitbGW4=contents['Ag8GWpercfitb'][0][0]
	elif ENM == 'CuO':
		if GWpH4==6:
			percfitaGW4=contents['CuO6GWpercfita'][0][0]
			percfitbGW4=contents['CuO6GWpercfitb'][0][0]
		elif GWpH4==6.5:
			percfitaGW4=contents['CuO65GWpercfita'][0][0]
			percfitbGW4=contents['CuO65GWpercfitb'][0][0]
		elif GWpH4==7:
			percfitaGW4=contents['CuO7GWpercfita'][0][0]
			percfitbGW4=contents['CuO7GWpercfitb'][0][0]
		elif GWpH4==7.5:
			percfitaGW4=contents['CuO75GWpercfita'][0][0]
			percfitbGW4=contents['CuO75GWpercfitb'][0][0]
		elif GWpH4==8:
			percfitaGW4=contents['CuO8GWpercfita'][0][0]
			percfitbGW4=contents['CuO8GWpercfitb'][0][0]
	elif ENM == 'nZVI':
		if GWpH4==6:
			percfitaGW4=contents['nZVI6GWpercfita'][0][0]
			percfitbGW4=contents['nZVI6GWpercfitb'][0][0]
		elif GWpH4==6.5:
			percfitaGW4=contents['nZVI65GWpercfita'][0][0]
			percfitbGW4=contents['nZVI65GWpercfitb'][0][0]
		elif GWpH4==7:
			percfitaGW4=contents['nZVI7GWpercfita'][0][0]
			percfitbGW4=contents['nZVI7GWpercfitb'][0][0]
		elif GWpH4==7.5:
			percfitaGW4=contents['nZVI75GWpercfita'][0][0]
			percfitbGW4=contents['nZVI75GWpercfitb'][0][0]
		elif GWpH4==8:
			percfitaGW4=contents['nZVI8GWpercfita'][0][0]
			percfitbGW4=contents['nZVI8GWpercfitb'][0][0]
	elif ENM == 'ZnO':
		if GWpH4==6:
			percfitaGW4=contents['ZnO6GWpercfita'][0][0]
			percfitbGW4=contents['ZnO6GWpercfitb'][0][0]
		elif GWpH4==6.5:
			percfitaGW4=contents['ZnO65GWpercfita'][0][0]
			percfitbGW4=contents['ZnO65GWpercfitb'][0][0]
		elif GWpH4==7:
			percfitaGW4=contents['ZnO7GWpercfita'][0][0]
			percfitbGW4=contents['ZnO7GWpercfitb'][0][0]
		elif GWpH4==7.5:
			percfitaGW4=contents['ZnO75GWpercfita'][0][0]
			percfitbGW4=contents['ZnO75GWpercfitb'][0][0]
		elif GWpH4==8:
			percfitaGW4=contents['ZnO8GWpercfita'][0][0]
			percfitbGW4=contents['ZnO8GWpercfitb'][0][0]
	elif ENM== 'TiO2':
		percfitaGW4=contents['TiO2GWpercfita'][0][0]
		percfitbGW4=contents['TiO2GWpercfitb'][0][0]
	elif ENM== 'CeO2':
		percfitaGW4=contents['CeO2GWpercfita'][0][0]
		percfitbGW4=contents['CeO2GWpercfitb'][0][0]
	elif ENM== 'SiO2':
		percfitaGW4=contents['SiO2GWpercfita'][0][0]
		percfitbGW4=contents['SiO2GWpercfitb'][0][0]

	DIS = {}
	DIS['percfitaFW']=percfitaFW
	DIS['percfitbFW']=percfitbFW
	DIS['percfitaSW']=percfitaSW
	DIS['percfitbSW']=percfitbSW
	DIS['percfitaGW1']=percfitaGW1
	DIS['percfitbGW1']=percfitbGW1
	DIS['percfitaGW2']=percfitaGW2
	DIS['percfitbGW2']=percfitbGW2
	DIS['percfitaGW3']=percfitaGW3
	DIS['percfitbGW3']=percfitbGW3
	DIS['percfitaGW4']=percfitaGW4
	DIS['percfitbGW4']=percfitbGW4
	
	return DIS
