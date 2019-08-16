
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors


class PartitionCoefficient:

    """
    Kd is species-specific, even compartment specific too
    Kaw - Kaw_n can be generated from EPI Suite
    Kow - Kow_n can be generated from EPI Suite
    """

    def __init__(self, chem_type):
        '''
        :param chem_type: organic acid or organic base or metal
        '''
        self.chem_type = chem_type


    def Kd_aero_n(self, Koa_n, Kaw_n, foc, solidP):
        '''
        :param Koa:
        :param foc:
        :param solidP:
        :return:

        Source: Franco, 2009, equation was developed by Harner and Bidleman and corrected as suggested by Gotz et al.
        '''
        Kd_val = 0.54 * Koa_n * Kaw_n * foc * solidP/1000.0
        return Kd_val


    def Kd_n(self, Koc_n, foc, solidP):
        '''
        It can be applied to sediment and soil
        :param Koc_n: neutral Koc
        :param foc: organic carbon concentent
        :param solidP: solid density (kg-soil/m3-soil)
        :return:
        '''
        Kd_n = Koc_n * foc * solidP/1000.0
        return Kd_n


    def Kd_j(self, Koc_j, foc, solidP):
        '''
        :param foc: organic carbon content of the dry matter (kg/kg)
        :param Koc: solid-water partition coefficient (L/kg)
        :param density_solid: kg/L, kg/m3

        :return: dimensionless solid-water sorption coefficient (L/L)
        '''

        Kd_val = Koc_j * foc * solidP / 1000.0

        return Kd_val


    def Kd_soil_organic_j(self, smiles, foc, CEC_soil, Koc_acid, Kow_n, pKa, Fr_n, Fr_i, solidP):
        """
        This is used for Kd values in soil.
        If the user doesn't input the values, we will calcuate them. Otherwise, use the user provided values.
        Kd_soil_n is handled in the load data function through Koc values
        :param smiles:
        :param foc: organic carbon fraction, value is between 0 to 1
        :param CEC_soil: in unit mol/kg
        :param Koc_acid: from a search function, if Koc_acid not empty, use Koc_acid, otherwise use the function results
        :return:
        """
        Kd_j = 0.0
        if self.chem_type == 'organic base':
            Kd_j = self.Kd_organic_base(smiles, CEC_soil, foc)
        elif self.chem_type == 'organic acid':
            if Koc_acid is None:
                Koc_j_val = self.Koc_j(Kow_n, pKa, Fr_n, Fr_i)
            else:
                Koc_j_val = Koc_acid
            Kd_j = self.Kd_j(foc, Koc_j_val, solidP)
        return Kd_j


    def Kd_soil_metal_j_regression(self, chem_formula, pH, total_metal, SOC):
        '''
        :param metal_str:
        :return:

        Source: Solid-solution partitioning of metals in contaminated soils: dependence on pH, total metal burden,
            and organic matter, 2000, ES&T
            https://pubs.acs.org/doi/abs/10.1021/es9907764
        '''
        chem_formula = chem_formula.lower()

        # check if SOC and total_metal are provided
        if isinstance(total_metal, float) and isinstance(SOC, float):
            # multiply 2.0 because organic matter is roughly twice of organic carbon
            if 'cd' in chem_formula:
                Kd_j = 10 ** (0.48 * pH + 0.82 * np.log10(SOC*2*100) - 0.65)
            elif 'cu' in chem_formula:
                Kd_j = 10 ** (0.21 * pH + 0.51 * np.log10(SOC*2*100) + 1.75)
            elif 'ni' in chem_formula:
                Kd_j = 10 ** (1.02 * pH + 0.80 * np.log10(SOC*2*100) - 4.16)
            elif 'pb' in chem_formula:
                Kd_j = 10 ** (0.37 * pH + 0.44 * np.log10(total_metal) + 1.19)
            elif 'zn' in chem_formula:
                Kd_j = 10 ** (0.60 * pH + 0.21 * np.log10(total_metal) - 1.34)
            else:
                # if the metal doesn't belong to any of the five (Cd, Cu, Ni, Pb, Zn)
                # use the regression model of Cd to represent it
                Kd_j = 10 ** (0.48 * pH + 0.82 * np.log10(SOC*2) - 0.65)
        else:
            # use simplier models with lower R2
            if 'cd' in chem_formula:
                Kd_j = 10 ** (0.49 * pH - 0.60)
            elif 'cu' in chem_formula:
                Kd_j = 10 ** (0.27 * pH + 1.49)
            elif 'ni' in chem_formula:
                Kd_j = 10 ** (0.72 * pH - 1.75)
            elif 'pb' in chem_formula:
                Kd_j = 10 ** (0.49 * pH + 1.37)
            elif 'zn' in chem_formula:
                Kd_j = 10 ** (0.62 * pH - 0.97)
            else:
                # if the metal doesn't belong to any of the five (Cd, Cu, Ni, Pb, Zn)
                # use the regression model of Cd to represent it
                Kd_j = 10 ** (0.49 * pH - 0.60)
        return Kd_j



    def Kd_organic_base(self, smiles, CEC_soil, foc):
        '''
        Source: Development and evaluation of a new sorption model for organic cations in soil:
        contributions from organic matter and clay minerals
            https://pubs.acs.org/doi/10.1021/es4031886
        '''

        RDKit_dic = self.cal_descriptor(smiles)
        Vx = RDKit_dic['VMcGowan'] / 100.0
        NAi = RDKit_dic['NAi']
        # reference sorption coefficient for the clay fraction in soils
        K_cec_clay = 10.0 ** (1.22 * Vx - 0.22 * NAi + 1.09)
        # ion-exchange-based sorption coefficient to Pahokee peat
        D_oc_ie = 10.0 ** (1.53 * Vx + 0.32 * NAi - 0.27)
        Kd_j = K_cec_clay * (CEC_soil - 3.4 * foc) + foc * D_oc_ie

        return Kd_j


    # Load RDKit data to the exposure table in JSON
    def cal_descriptor(self, smiles):
        # call rdkit mordred
        # McGowan's Volume and number of hydrogens bound by the charged nitrogen
        RDKIT_KEYS = ["VMcGowan"]
        mols = [Chem.MolFromSmiles(smi) for smi in [smiles]]
        calc = Calculator(descriptors, ignore_3D=True)
        df = calc.pandas(mols, nproc=1)
        mol = Chem.MolFromSmiles(smiles)
        NH0 = Chem.Fragments.fr_NH0(mol)  # number of Tertiary amines
        NH1 = Chem.Fragments.fr_NH1(mol)  # number of Secondary amines
        NH2 = Chem.Fragments.fr_NH2(mol)  # number of Primary amines

        NAi = 0
        if NH2 != 0:
            NAi = 3
        elif NH1 != 0:
            NAi = 2
        elif NH0 != 0:
            NAi = 1
        RDKit_dic = {}
        RDKit_dic['NAi'] = NAi
        for key in RDKIT_KEYS:
            RDKit_dic[key] = float(df[key][0])
        return RDKit_dic


    def Koc_j(self, Kow_n, pKa, Fr_n, Fr_i):
        '''
        for chem_type is organic base or organic acid, this equation is used for organic acid if not exp data avalable
        :param Kow_n: octanol-water partition coefficient of neutral species
        :param Fr_n: fraction of neutral species in water phase
        :param Fr_i: fraction of ionic species in water phase
        :param pKa: dissociation acid constant

        :return:

        # the equation is applicable for pKa (0,12), pKb (2,12)
        # logKow_n (-2.18, 8.5) for acides, (-1.66, 7.03) for bases

        Source: A multimedia activity model for ionizable compounds - validation study with 2,4-D, aniline and trimethoprim
            Antonio Franco and Stefan Trapp, 2010
            https://www.ncbi.nlm.nih.gov/pubmed/20821507

        '''

        Kow_i = 10.0 ** (np.log10(Kow_n) - 3.5)
        Kow_apparent = Fr_n * Kow_n + Fr_i * Kow_i

        Koc_val = 0.0
        if self.chem_type == 'organic acid':
            Koc_val = 10.0 ** (0.11 * np.log10(Kow_n) + 1.54)

        elif self.chem_type == 'organic base':
            f = Kow_apparent / (Kow_apparent + 1)
            Koc_val = (pKa ** 0.65) * (f ** 0.14)

        return Koc_val
