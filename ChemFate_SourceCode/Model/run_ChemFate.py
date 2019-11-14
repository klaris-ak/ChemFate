import os

from Code.model_setup import Model_SetUp

chem_type = 'Nanomaterial' # can only be 'NonionizableOrganic', 'IonizableOrganic', 'Metal', or 'Nanomaterial'
CUR_PATH = os.path.dirname(os.path.abspath(__file__))

# chem_file = CUR_PATH + '/../Input/ChemParam_ethylbenzene.xlsx'

if chem_type == 'NonionizableOrganic':
    chem_file = CUR_PATH + '/../Input/ChemParam_nonionizableOrganic.xlsx'
elif chem_type == 'IonizableOrganic':
    chem_file = CUR_PATH + '/../Input/ChemParam_ionizableOrganic.xlsx'
elif chem_type == 'Metal':
    chem_file = CUR_PATH + '/../Input/ChemParam_metal.xlsx'
else:
    chem_file = CUR_PATH + '/../Input/ChemParam_nanomaterial.xlsx'


region_file = CUR_PATH + '/../Input/Region.xlsx'
release_file = CUR_PATH + '/../Input/ChemRelease.xlsx'
run_option = 1 # can be 1 or 2:
bgPercOption2 = 10 # can be anywhere between 0-100
output_file_path = CUR_PATH + '/../Output/'
file_name = 'nanoCopper' # pyrimethanil, cyprodinil, copper, nanoCopper

# year month day
start_date = '2005 1 1'
end_date = '2005 1 31'


model = Model_SetUp(start_date, end_date, run_option, bgPercOption2,
                    chem_type, chem_file, region_file, release_file, output_file_path, file_name)
model.run_model()

