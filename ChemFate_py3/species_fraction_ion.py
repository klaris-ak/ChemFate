
class SpeciesFraction:

    def __init__(self, pKa, chem_type):
        '''
        :param pKa: the negative lag of the acid dissociation constant

        species fraction equations are used for organic acid or base
        for metals, species fraction need to be generated from WHAM or MINTEQ software
        '''
        self.pKa = pKa
        self.chem_type = chem_type
        pass


    def species_fraction_water(self, pH):
        '''
        :param species_type: neutral, cation, and anion
        :param pH: the pH in the compartment
        :return:

        Source: A multimedia activity model for ionizable compounds:
            validation study with 2,4-dichlorophenoxyacetic acid, aniline, and trimethoprim.
            Franco, 2009, SI
            (https://www.ncbi.nlm.nih.gov/pubmed/20821507)

        applicapable in: aerosol water, freshwater, seawater, sediment water, soil pore water

        '''

        Fr_n = 0.0 # netural fraction (Qn)
        Fr_i = 0.0 # anionic fraction (Q-, ionic acid)

        # Henderson-Hasselbalch equation (Henderson, 1908)

        if self.chem_type == 'organic acid': # cation
            Fr_n = 1.0/(1.0 + (10.0 ** (pH - self.pKa)))
            Fr_i = 1.0 - Fr_n
        elif self.chem_type == 'organic base': # anion
            Fr_n = 1.0/(1.0 + (10.0 ** (self.pKa - pH)))
            Fr_i = 1.0 - Fr_n

        return Fr_n, Fr_i



    def species_fraction_solid(self, pH):
        '''
        :param pH: the pH in the compartment
        :return:

        Source: Franco and Trapp (2008)
        applicapable in: sediment solid and soil solid

        '''

        Fr_n = 0.0 # netural fraction (Qn)
        Fr_i = 0.0 # anionic fraction (Q-, ionic acid)

        # Henderson-Hasselbalch equation (Henderson, 1908)
        if self.chem_type == 'organic acid': # cation
            Fr_n = 1.0/(1.0 + (10.0 ** (pH - 0.6 - self.pKa)))
            Fr_i = 1.0 - Fr_n
        elif self.chem_type == 'organic base': # anion
            Fr_n = 1.0/(1.0 + (10.0 ** (self.pKa - 4.5)))
            Fr_i = 1.0 - Fr_n

        return Fr_n, Fr_i