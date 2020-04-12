
import json
from collections import OrderedDict

from species_fraction_ion import SpeciesFraction
from partition_coef_ion import PartitionCoefficient


with open('./IonizableChem_helper.json') as f:
    data = json.load(f)

foc_map = data['foc_map']
solidP_map = data['solidP_map']
CEC_map = data['CEC_map']
pH_map = data['pH_map']
metal_map = data['total_metal_map']

compart_list = ['air', 'fw', 'fwSed', 'sw', 'swSed', 'soil1', 'deepS1', 'soil2', 'deepS2',
                'soil3', 'deepS3', 'soil4', 'deepS4']
soil_compart_list = ['soil1', 'soil2', 'soil3', 'soil4']
deepSoil_compart_list = ['deepS1', 'deepS2', 'deepS3', 'deepS4']

class Y_Value:

    def __init__(self, chem_type, chemParams, env_dict):
        self.chem_type = chem_type
        self.chemParams = chemParams
        self.env_dict = env_dict


    ''' Species fraction '''
    def X_ij(self):
        solid_compart_list = ['fwSed', 'swSed', 'soil1', 'soil2', 'soil3', 'soil4']
        # for ionizable organic chems, the list contains two types of species: neutral and ionic
        # for metals, the list contains three types of species: 1) attached to particles;
        # 2) attached to dissolved organic matters; 3) dissolved ions
        X_dict = OrderedDict()
        for compart in compart_list:
            X_dict[compart] = []

        # compart_txt_list = ['air']
        # check if user input species fraction value
        partition = PartitionCoefficient(self.chemParams['type'])
        if self.chem_type != 'Metal':
            speciation = SpeciesFraction(self.chemParams['pKa'], self.chemParams['type'])
            for compart in compart_list:
                Fr_txt_n = 'Fr_' + compart + '_n'
                Fr_txt_i = 'Fr_' + compart + '_i'

                # if self.chemParams[Fr_txt_n] == '--':
                if isinstance(self.chemParams[Fr_txt_n], str):
                    # if compart not in solid_compart_list:
                    #     Fr_n, Fr_i = speciation.species_fraction_water(self.env_dict[pH_map[compart]])
                    # else:
                    #     Fr_n, Fr_i = speciation.species_fraction_solid(self.env_dict[pH_map[compart]])
                    Fr_n, Fr_i = speciation.species_fraction_water(self.env_dict[pH_map[compart]])
                    X_dict[compart].append(Fr_n)
                    X_dict[compart].append(Fr_i)
                else:
                    X_dict[compart].append(self.chemParams[Fr_txt_n])
                    X_dict[compart].append(self.chemParams[Fr_txt_i])

        if self.chem_type == 'Metal':
            for compart in compart_list:
                Fr_txt_p = 'Fr_' + compart + '_p'
                Fr_txt_c = 'Fr_' + compart + '_c'
                Fr_txt_i = 'Fr_' + compart + '_i'

                if compart not in deepSoil_compart_list and isinstance(self.chemParams[Fr_txt_p], float):
                    X_dict[compart].append(self.chemParams[Fr_txt_p] / 100.0)
                    X_dict[compart].append(self.chemParams[Fr_txt_c] / 100.0)
                    X_dict[compart].append(self.chemParams[Fr_txt_i] / 100.0)

                if compart in soil_compart_list and not isinstance(self.chemParams[Fr_txt_p], float):
                    Kd_j = partition.Kd_soil_metal_j_regression(self.chemParams['chem_formula'],
                                                                self.env_dict[pH_map[compart]],
                                                                self.env_dict[metal_map[compart]],
                                                                self.env_dict[foc_map[compart]])

                    #Fr_txt_c need to be provided in the input file.
                    Fr_txt_p = Kd_j/(Kd_j + 1.0)
                    b = self.chemParams[Fr_txt_c]/self.chemParams[Fr_txt_i]
                    Fr_txt_c = b/((Kd_j+1.0)*(b+1.0))
                    Fr_txt_i = 1.0/((Kd_j+1.0)*(b+1.0))
                    X_dict[compart].append(Fr_txt_p)
                    X_dict[compart].append(Fr_txt_c)
                    X_dict[compart].append(Fr_txt_i)

                # assume deep soil species fraction distributions are the same as surface soil
                if compart in soil_compart_list:
                    index = soil_compart_list.index(compart)
                    X_dict[deepSoil_compart_list[index]] = X_dict[compart]

        return X_dict


    ''' Partition coefficient '''
    def Kd_ij(self):
        X_dict = self.X_ij()
        partition = PartitionCoefficient(self.chemParams['type'])
        # each list contains two values, one for neutral and the other one for ionic
        Kd_dict = OrderedDict()
        for compart in compart_list:
            Kd_dict[compart] = [0, 0]

        # calculate Kd from partition coefficient script
        if self.chem_type != 'Metal':
            if isinstance(self.chemParams['Kp_n'], float):
                Kd_dict['air'][0] = self.chemParams['Kp_n']
            # elif self.chemParams['Koa_n'] != '--':
            elif isinstance(self.chemParams['Koa_n'], float):
                Kd_dict['air'][0] = partition.Kd_aero_n(self.chemParams['Koa_n'], self.chemParams['Kaw_n'],
                                                        self.env_dict[foc_map['air']], self.env_dict[solidP_map['air']])
            # for rest of the compartment neural species Kd can use a loop
            Kd_n_compart = ['fw', 'fwSed', 'sw', 'swSed', 'soil1', 'deepS1', 'soil2', 'deepS2', 'soil3',
                            'deepS3', 'soil4', 'deepS4']
            for compart in Kd_n_compart:
                Kd_dict[compart][0] = partition.Kd_n(self.chemParams['Koc_n'],
                                                     self.env_dict[foc_map[compart]],
                                                     self.env_dict[solidP_map[compart]])
            # Kd [1] ionic
            Kd_i_compart = ['air', 'fw', 'fwSed', 'sw', 'swSed', 'soil1', 'deepS1', 'soil2', 'deepS2', 'soil3',
                            'deepS3', 'soil4', 'deepS4']
            Kd_i_water = ['air', 'fw', 'fwSed', 'sw', 'swSed']
            X_dict_2 = X_dict.copy()
            # Kd_i value are optional for user to input, so check first if exist use, if not, calculate it
            for compart in Kd_i_compart:
                key_txt = 'Kd_' + compart + '_i'

                if isinstance(self.chemParams[key_txt], float):
                    Kd_dict[compart][1] = self.chemParams[key_txt]
                else:
                    # calculate the value from equations
                    # check the compartment water or soil, aerosl contains water too
                    if X_dict_2[compart][1] == 0:
                        Kd_j = 0
                    else:
                        if compart in Kd_i_water:
                            Koc_j = partition.Koc_j(self.chemParams['Kow_n'], self.chemParams['pKa'],
                                                    X_dict_2[compart][0], X_dict_2[compart][1])

                            Kd_j = partition.Kd_j(Koc_j, self.env_dict[foc_map[compart]],
                                                  self.env_dict[solidP_map[compart]])

                        else:
                            if 'base' in self.chemParams['type'] or 'metal' in self.chemParams['type']:
                                Koc_acid = None
                            else:
                                # check if Koc in the Koc_organic_acid_db
                                Koc_acid = self.chemParams['Koc_acid'] # it might still be None is not in the db

                            Kd_j = partition.Kd_soil_organic_j(self.chemParams['smiles'],
                                                               self.env_dict[foc_map[compart]],
                                                               self.env_dict[CEC_map[compart]],
                                                               Koc_acid, self.chemParams['Kow_n'],
                                                               self.chemParams['pKa'],
                                                               X_dict_2[compart][0], X_dict_2[compart][1],
                                                               self.env_dict[solidP_map[compart]])
                    Kd_dict[compart][1] = Kd_j
        else:
            # use the species fraction values for metal to calculate Kd's
            # Kd_dict[compart][0] is KdPD, Kd_dict[compart][1] is KdCD
            for compart in compart_list:
                Kd_dict[compart][0] = X_dict[compart][0] / X_dict[compart][2]
                Kd_dict[compart][1] = X_dict[compart][1] / X_dict[compart][2]

        return Kd_dict


    ''' Aquivalence capacity '''
    def Z_ij(self):
        """
        Z values are species-specific and can be calculated using species-specific partition coefficients;

        if the chemical is metal
            - particulate form, colloidal form, dissolved form
            - Z values: [Zp, Zc, Zw=1]

        """
        Kd_dict = self.Kd_ij()
        # Z_air = Z_water * Kaw_n, for neutral phase, ionic species = 0
        # Z_air doesn't exist. Z_water = 1
        # Zsolid = Zwater * Kd = Kd
        Z_air = self.chemParams['Kaw_n']
        Z_dict = OrderedDict()
        Z_dict_sub = {}

        if self.chem_type != 'Metal':
            Z_dict_sub = {'air': [Z_air, 0.0, 1.0], 'aer': [0.0, 0.0, 1.0], 'fw': [1.0, 1.0, 1.0], 'fSS': [0.0, 0.0, 1.0],
                          'fSedW': [1.0, 1.0, 1.0], 'fSedS': [0.0, 0.0, 1.0], 'sw': [1.0, 1.0, 1.0], 'sSS': [0.0, 0.0, 1.0],
                          'sSedW': [1.0, 1.0, 1.0], 'sSedS': [0.0, 0.0, 1.0], 'soilA1': [Z_air, 0.0, 1.0],
                          'soilW1': [1.0, 1.0, 1.0], 'soilS1': [0.0, 0.0, 1.0], 'deepS1': [0.0, 0.0, 1.0],
                          'soilA2': [Z_air, 0.0, 1.0], 'soilW2': [1.0, 1.0, 1.0], 'soilS2': [0.0, 0.0, 1.0],
                          'deepS2': [0.0, 0.0, 1.0], 'soilA3': [Z_air, 0.0, 1.0], 'soilW3': [1.0, 1.0, 1.0],
                          'soilS3': [0.0, 0.0, 1.0], 'deepS3': [0.0, 0.0, 1.0], 'soilA4': [Z_air, 0.0, 1.0],
                          'soilW4': [1.0, 1.0, 1.0], 'soilS4': [0.0, 0.0, 1.0], 'deepS4': [0.0, 0.0, 1.0]}
            Z_dict_sub = OrderedDict(Z_dict_sub)
            Z_dict_sub['aer'][0:2] = Kd_dict['air']
            Z_dict_sub['fSS'][0:2] = Kd_dict['fw']
            Z_dict_sub['fSedS'][0:2] = Kd_dict['fwSed']
            Z_dict_sub['sSS'][0:2] = Kd_dict['sw']
            Z_dict_sub['sSedS'][0:2] = Kd_dict['swSed']
            Z_dict_sub['soilS1'][0:2] = Kd_dict['soil1']
            Z_dict_sub['deepS1'][0:2] = Kd_dict['deepS1']
            Z_dict_sub['soilS2'][0:2] = Kd_dict['soil2']
            Z_dict_sub['deepS2'][0:2] = Kd_dict['deepS2']
            Z_dict_sub['soilS3'][0:2] = Kd_dict['soil3']
            Z_dict_sub['deepS3'][0:2] = Kd_dict['deepS3']
            Z_dict_sub['soilS4'][0:2] = Kd_dict['soil4']
            Z_dict_sub['deepS4'][0:2] = Kd_dict['deepS4']


            for compart in compart_list:
                Z_dict[compart] = [0.0, 0.0]
            for i in range(0, 2):
                Z_dict['air'][i] = self.env_dict['airVf'] * Z_dict_sub['air'][i] + \
                                   self.env_dict['aerVf'] * Z_dict_sub['aer'][i]
                Z_dict['fw'][i] = self.env_dict['fwVf'] * Z_dict_sub['fw'][i] + \
                                  self.env_dict['fSSVf'] * Z_dict_sub['fSS'][i]
                Z_dict['fwSed'][i] = (1 - self.env_dict['fsedpercSolid']) * Z_dict_sub['fSedW'][i] + \
                                     self.env_dict['fsedpercSolid'] * Z_dict_sub['fSedS'][i]
                Z_dict['sw'][i] = self.env_dict['swVf'] * Z_dict_sub['sw'][i] + self.env_dict['sSSVf'] * \
                                  Z_dict_sub['sSS'][i]
                Z_dict['swSed'][i] = (1 - self.env_dict['ssedpercSolid']) * Z_dict_sub['sSedW'][i] + \
                                     self.env_dict['ssedpercSolid'] * Z_dict_sub['sSedS'][i]
                Z_dict['soil1'][i] = self.env_dict['soilAC1'] * Z_dict_sub['soilA1'][i] + \
                                     self.env_dict['soilWC1'] * Z_dict_sub['soilW1'][i] + \
                                     self.env_dict['soilSC1'] * Z_dict_sub['soilS1'][i]
                Z_dict['deepS1'][i] = Z_dict_sub['deepS1'][i]
                Z_dict['soil2'][i] = self.env_dict['soilAC2'] * Z_dict_sub['soilA2'][i] + \
                                     self.env_dict['soilWC2'] * Z_dict_sub['soilW2'][i] + \
                                     self.env_dict['soilSC2'] * Z_dict_sub['soilS2'][i]
                Z_dict['deepS2'][i] = Z_dict_sub['deepS2'][i]
                Z_dict['soil3'][i] = self.env_dict['soilAC3'] * Z_dict_sub['soilA3'][i] + \
                                     self.env_dict['soilWC3'] * Z_dict_sub['soilW3'][i] + \
                                     self.env_dict['soilSC3'] * Z_dict_sub['soilS3'][i]
                Z_dict['deepS3'][i] = Z_dict_sub['deepS3'][i]
                Z_dict['soil4'][i] = self.env_dict['soilAC4'] * Z_dict_sub['soilA4'][i] + \
                                     self.env_dict['soilWC4'] * Z_dict_sub['soilW4'][i] + \
                                     self.env_dict['soilSC4'] * Z_dict_sub['soilS4'][i]
                Z_dict['deepS4'][i] = Z_dict_sub['deepS4'][i]
        elif self.chem_type == 'Metal':
            for compart in compart_list:
                # the dissolved for is purely in the water, Zw = 1.0 by definition
                Z_dict[compart] = [0.0, 0.0, 1.0]
                Z_dict[compart][0:2] = Kd_dict[compart]
        return Z_dict, Z_dict_sub


    ''' Aquivalence fraction '''
    def Y_ij(self):
        X_dict = self.X_ij()
        Z_dict, Z_dict_sub = self.Z_ij()
        Y_dict = OrderedDict()

        for compart in compart_list:
            if self.chem_type != 'Metal':
                # pre-define the aquivaence fraction, using air phase, neutral fraction is 1
                Y_dict[compart] = [1.0, 0.0]
                
                # Yij = (Xij/Zij)/sum(Xij/Zij), summation of neural and ionic species
                # avoid bottom value is 0
                if Z_dict[compart][0] == 0:
                    Z_dict[compart][0] = 10**(-20)
                if Z_dict[compart][1] == 0:
                    Z_dict[compart][1] = 10**(-20)
                sum_val = X_dict[compart][0]/Z_dict[compart][0] + X_dict[compart][1]/Z_dict[compart][1]
                Y_dict[compart][0] = (X_dict[compart][0] / Z_dict[compart][0]) / sum_val
                Y_dict[compart][1] = (X_dict[compart][1] / Z_dict[compart][1]) / sum_val

            elif self.chem_type == 'Metal':
                Y_dict[compart] = [0.0, 0.0, 0.0]
                # Yij = (Xij/Zij)/sum(Xij/Zij), summation of neural and ionic species
                # avoid bottom value is 0
                if Z_dict[compart][0] == 0:
                    Z_dict[compart][0] = 10 ** (-20)
                if Z_dict[compart][1] == 0:
                    Z_dict[compart][1] = 10 ** (-20)
                if Z_dict[compart][2] == 0:
                    Z_dict[compart][2] = 10 ** (-20)

                sum_val = X_dict[compart][0]/Z_dict[compart][0] + X_dict[compart][1]/Z_dict[compart][1] + X_dict[compart][2]/Z_dict[compart][2]

                Y_dict[compart][0] = (X_dict[compart][0] / Z_dict[compart][0]) / sum_val
                Y_dict[compart][1] = (X_dict[compart][1] / Z_dict[compart][1]) / sum_val
                Y_dict[compart][2] = (X_dict[compart][2] / Z_dict[compart][2]) / sum_val

        return Y_dict


    ''' Z value by each compartment with all species '''

    def Z_i(self):
        Z_ij_dict, Z_ij_dict_sub = self.Z_ij()
        Y_ij_dict = self.Y_ij()

        Z_i_dict = OrderedDict()
        for compart in compart_list:
            if self.chem_type != 'Metal':
                Z_i_dict[compart] = Y_ij_dict[compart][0] * Z_ij_dict[compart][0] + \
                                    Y_ij_dict[compart][1] * Z_ij_dict[compart][1]
            else:
                Z_i_dict[compart] = Y_ij_dict[compart][0] * Z_ij_dict[compart][0] + \
                                    Y_ij_dict[compart][1] * Z_ij_dict[compart][1] + \
                                    Y_ij_dict[compart][2] * Z_ij_dict[compart][2]

        return Z_i_dict



