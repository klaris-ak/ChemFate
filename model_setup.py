from __future__ import division
from datetime import datetime
import pandas as pd
from collections import OrderedDict
from load_data import LoadData
from load_data_nano import load_data
from model_solver import org_solver, ion_solver, nano_solver
from generate_result import GenerateResult


class Model_SetUp:

    def __init__(self, start_date, end_date, run_option, bgPercOption2,
                 chem_type, chem_file, region_file, release_file, output_file_path, file_name):
        # start date and end date need to be in the format of "%Y %m %d", eg:'2005 2 3'
        # option contains two options
        # option 1 - set background concentration to 0 or front end replace the concentration sheet data directly
        # option 2 - set background concentration to 0 first, and then run the model;
        # and then calculate the average concentrations and set it to the background concentration
        # and then run the model again

        self.start_date = start_date
        self.end_date = end_date
        self.run_option = run_option
        self.bgPercOption2 = bgPercOption2/100
        self.chem_type = chem_type
        self.chem_file = chem_file
        self.region_file = region_file
        self.release_file = release_file
        self.output_file_path = output_file_path
        self.file_name = file_name

    def simulation_days(self):
        start_day = datetime.strptime(self.start_date, "%Y %m %d")
        end_day = datetime.strptime(self.end_date, "%Y %m %d")
        sim_days = (end_day - start_day).days + 1
        return sim_days


    def bgConc_new_cal(self, input_file_list, HL_air, HL_fWater, HL_fSedS):
        # TODO: do we also need to consider the ion's bgconc?
        sim_days = self.simulation_days()
        rows_to_skip = int(sim_days * (1 - self.bgPercOption2))
        bgConc_new = {}
        conc_data_n = pd.read_csv(input_file_list[0], skiprows = range(1, rows_to_skip+1), skipinitialspace=True)

        bgConc_new['air'] = conc_data_n.air_conc.mean()
        bgConc_new['aer'] = conc_data_n.aerosol_conc.mean()
        bgConc_new['fw'] = conc_data_n.fw_conc.mean()
        bgConc_new['fSS'] = conc_data_n.fw_sus_sed_conc.mean()
        bgConc_new['fSedS'] = conc_data_n.fw_sed_solid.mean()
        bgConc_new['fSedW'] = conc_data_n.fw_sed_water.mean()
        bgConc_new['sw'] = conc_data_n.sw_conc.mean()
        bgConc_new['sSS'] = conc_data_n.sw_sus_sed_conc.mean()
        bgConc_new['sSedS'] = conc_data_n.sw_sed_solid.mean()
        bgConc_new['sSedW'] = conc_data_n.sw_sed_water.mean()
        bgConc_new['soilA1'] = conc_data_n.undeveloped_soil_air.mean()
        bgConc_new['soilW1'] = conc_data_n.undeveloped_soil_water.mean()
        bgConc_new['soilS1'] = conc_data_n.undeveloped_soil_solid.mean()
        bgConc_new['dsoil1'] = conc_data_n.deep_undeveloped_soil.mean()
        bgConc_new['soilA2'] = conc_data_n.urban_soil_air.mean()
        bgConc_new['soilW2'] = conc_data_n.urban_soil_water.mean()
        bgConc_new['soilS2'] = conc_data_n.urban_soil_solid.mean()
        bgConc_new['dsoil2'] = conc_data_n.deep_urban_soil.mean()
        bgConc_new['soilA3'] = conc_data_n.agricultural_soil_air.mean()
        bgConc_new['soilW3'] = conc_data_n.agricultural_soil_water.mean()
        bgConc_new['soilS3'] = conc_data_n.agricultural_soil_solid.mean()
        bgConc_new['dsoil3'] = conc_data_n.deep_agricultural_soil.mean()
        bgConc_new['soilA4'] = conc_data_n.biosolids_soil_air.mean()
        bgConc_new['soilW4'] = conc_data_n.biosolids_soil_water.mean()
        bgConc_new['soilS4'] = conc_data_n.biosolids_soil_solid.mean()
        bgConc_new['dsoil4'] = conc_data_n.deep_biosolids_soil.mean()

        if HL_air <= 24:
            # if the halflife of air is less than 1 day (24hr)
            # set the external air/aerosol to 20% of the compartment air/aerosol concentration
            bgConc_new['gairc_n'] = bgConc_new['air'] * 0.2
        elif HL_air > 24 and HL_air <= 672:
            # if the halflife of air is more than 1 day and less than 4 weeks
            # set the external air/aerosol to 50% of the compartment air/aerosol concentration
            bgConc_new['gairc_n'] = bgConc_new['air'] * 0.5
        else:
            # if the halflife of air is more than 4 weeks (1 month)
            # set the external air/aerosol to 80% of the compartment air/aerosol concentration
            bgConc_new['gairc_n'] = bgConc_new['air'] * 0.8

        # freshwater external concentration
        if HL_fWater <= 24:
            bgConc_new['gfreshwc_n'] = bgConc_new['fWater'] * 0.2
        elif HL_fWater > 24 and HL_fWater <= 672:
            bgConc_new['gfreshwc_n'] = bgConc_new['fWater'] * 0.5
        else:
            bgConc_new['gfreshwc_n'] = bgConc_new['fWater'] * 0.8

        # freshwater sediment external concentraton
        if HL_fSedS <= 24:
            bgConc_new['gfSedc_n'] = bgConc_new['fSedS'] * 0.2
        elif HL_fWater > 24 and HL_fWater <= 672:
            bgConc_new['gfSedc_n'] = bgConc_new['fSedS'] * 0.5
        else:
            bgConc_new['gfSedc_n'] = bgConc_new['fSedS'] * 0.8

        if self.chem_type != 'NonionizableOrganic':
            bgConc_new_i = {}
            conc_data_i = pd.read_csv(input_file_list[1], skiprows=range(1, rows_to_skip + 1), skipinitialspace=True)
            bgConc_new_i['air'] = conc_data_i.air_conc.mean()
            bgConc_new_i['aer'] = conc_data_i.aerosol_conc.mean()
            bgConc_new_i['fWater'] = conc_data_i.fw_conc.mean()
            bgConc_new_i['fSS'] = conc_data_i.fw_sus_sed_conc.mean()
            bgConc_new_i['fSedS'] = conc_data_i.fw_sed_solid.mean()
            bgConc_new_i['fSedW'] = conc_data_i.fw_sed_water.mean()
            bgConc_new_i['sWater'] = conc_data_i.sw_conc.mean()
            bgConc_new_i['sSS'] = conc_data_i.sw_sus_sed_conc.mean()
            bgConc_new_i['sSedS'] = conc_data_i.sw_sed_solid.mean()
            bgConc_new_i['sSedW'] = conc_data_i.sw_sed_water.mean()
            bgConc_new_i['soilA1'] = conc_data_i.undeveloped_soil_air.mean()
            bgConc_new_i['soilW1'] = conc_data_i.undeveloped_soil_water.mean()
            bgConc_new_i['soilS1'] = conc_data_i.undeveloped_soil_solid.mean()
            bgConc_new_i['soilDeep1'] = conc_data_i.deep_undeveloped_soil.mean()
            bgConc_new_i['soilA2'] = conc_data_i.urban_soil_air.mean()
            bgConc_new_i['soilW2'] = conc_data_i.urban_soil_water.mean()
            bgConc_new_i['soilS2'] = conc_data_i.urban_soil_solid.mean()
            bgConc_new_i['soilDeep2'] = conc_data_i.deep_urban_soil.mean()
            bgConc_new_i['soilA3'] = conc_data_i.agricultural_soil_air.mean()
            bgConc_new_i['soilW3'] = conc_data_i.agricultural_soil_water.mean()
            bgConc_new_i['soilS3'] = conc_data_i.agricultural_soil_solid.mean()
            bgConc_new_i['soilDeep3'] = conc_data_i.deep_agricultural_soil.mean()
            bgConc_new_i['soilA4'] = conc_data_i.biosolids_soil_air.mean()
            bgConc_new_i['soilW4'] = conc_data_i.biosolids_soil_water.mean()
            bgConc_new_i['soilS4'] = conc_data_i.biosolids_soil_solid.mean()
            bgConc_new_i['soilDeep4'] = conc_data_i.deep_biosolids_soil.mean()

            compart_list = ['air', 'aer', 'fWater', 'fSS', 'fSedS', 'fSedW', 'sWater', 'sSS', 'sSedS', 'sSedW',
                            'soilA1', 'soilW1', 'soilS1', 'soilDeep1', 'soilA2', 'soilW2', 'soilS2', 'soilDeep2',
                            'soilA3', 'soilW3', 'soilS3', 'soilDeep3', 'soilA4', 'soilW4', 'soilS4', 'soilDeep4']
            for compart in compart_list:
                bgConc_new[compart] = bgConc_new[compart] + bgConc_new_i[compart]


            # freshwater external concentration
            if HL_fWater <= 24:
                bgConc_new['gfreshwc_i'] = bgConc_new_i['fWater'] * 0.2
            elif HL_fWater > 24 and HL_fWater <= 672:
                bgConc_new['gfreshwc_i'] = bgConc_new_i['fWater'] * 0.5
            else:
                bgConc_new['gfreshwc_i'] = bgConc_new_i['fWater'] * 0.8

            # freshwater sediment external concentraton
            if HL_fSedS <= 24:
                bgConc_new['gfSedc_i'] = bgConc_new_i['fSedS'] * 0.2
            elif HL_fWater > 24 and HL_fWater <= 672:
                bgConc_new['gfSedc_i'] = bgConc_new_i['fSedS'] * 0.5
            else:
                bgConc_new['gfSedc_i'] = bgConc_new_i['fSedS'] * 0.8

        return bgConc_new

    def bgConc_new_cal_nano(self, input_file, time):
        rows_to_skip = int(time * (1 - self.bgPercOption2))
        bgConc_new = OrderedDict()
        conc_data = pd.read_csv(input_file, skiprows = range(1, rows_to_skip+1), skipinitialspace=True)
        # does not require unit conversions
        bgConc_new['A'] = conc_data.Air.mean()
        bgConc_new['Aer'] = conc_data.Aerosols.mean()
        bgConc_new['fW'] = conc_data.Freshwater.mean()
        bgConc_new['fSS'] = conc_data.Freshwater_SS.mean()
        bgConc_new['fwSed'] = conc_data.Freshwater_Sed.mean()
        bgConc_new['sW'] = conc_data.Marine.mean()
        bgConc_new['sSS'] = conc_data.Marine_SS.mean()
        bgConc_new['swSed'] = conc_data.Marine_Sed.mean()
        bgConc_new['S1'] = conc_data.UndevSoil_Solid.mean()
        bgConc_new['soilW1'] = conc_data.UndevSoil_W.mean()
        bgConc_new['dsoil1'] = conc_data.UndevDeepSoil.mean()
        bgConc_new['S2'] = conc_data.UrbanSoil_Solid.mean()
        bgConc_new['soilW2'] = conc_data.UrbanSoil_W.mean()
        bgConc_new['dsoil2'] = conc_data.UrbanDeepSoil.mean()
        bgConc_new['S3'] = conc_data.AgSoil_Solid.mean()
        bgConc_new['soilW3'] = conc_data.AgSoil_W.mean()
        bgConc_new['dsoil3'] = conc_data.AgDeepSoil.mean()
        bgConc_new['S4'] = conc_data.BioSoil_Solid.mean()
        bgConc_new['soilW4'] = conc_data.BioSoil_W.mean()
        bgConc_new['dsoil4'] = conc_data.BioDeepSoil.mean()
        bgConc_new['fWdis'] = conc_data.Freshwater_dis.mean()
        bgConc_new['fWSeddis'] = conc_data.FreshwaterSed_Dis.mean()
        bgConc_new['Swdis'] = conc_data.Marine_dis.mean()
        bgConc_new['swSeddis'] = conc_data.MarineSed_dis.mean()
        bgConc_new['soilW1dis'] = conc_data.UndevSoilW_dis.mean()
        bgConc_new['soilW2dis'] = conc_data.UrbanSoilW_dis.mean()
        bgConc_new['soilW3dis'] = conc_data.AgSoilW_dis.mean()
        bgConc_new['soilW4dis'] = conc_data.BioSolidW_dis.mean()
        bgConc_new['gairc'] = conc_data.Air.mean()*0.05
        bgConc_new['gaerc'] =  conc_data.Aerosols.mean()*0.05

        return bgConc_new


    def run_model(self):
        if self.chem_type != 'Nanomaterial':
            # load data
            data = LoadData(self.chem_type, self.chem_file, self.region_file, self.release_file, self.start_date,
                            self.end_date)
            chemParams, presence, env, climate, bgConc, release, release_scenario = data.run_loadData()
            V_bulk_list = [env['airV'], env['fwV'], env['sedFWV'], env['swV'], env['sedSWV'], env['soilV1'],
                           env['soilV2'], env['soilV3'], env['soilV4']]
            sim_days = self.simulation_days()
            funC_df_list = []
            funM_df_list = []
        else:
            # load data
            time, presence, env, climate, bgConc, chemParams, release, release_scenario = load_data(self.region_file, self.release_file,
                                                                           self.chem_file, self.start_date,
                                                                           self.end_date)

            V_bulk_list = [env['airV'], env['freshwV'], env['sedFWV'], env['seawV'], env['sedSWV'], env['soilV1'],
                           env['soilV2'], env['soilV3'], env['soilV4']]

        if self.run_option == 1:
            if self.chem_type == 'NonionizableOrganic':
                date_array, process_array, funC_kg_1, funC_kg_1_sub, funM_kg_1, funM_kg_1_sub = \
                    org_solver(self.start_date, sim_days, presence, env, climate, chemParams, bgConc, release)
                funC_df_list = [funC_kg_1, funC_kg_1_sub]
                funM_df_list = [funM_kg_1, funM_kg_1_sub]

            elif self.chem_type == 'IonizableOrganic' or self.chem_type=='Metal':
                date_array, process_array, funC_kg_1, funC_kg_2, funC_kg_3, funC_kg_1_sub, funC_kg_2_sub, funC_kg_3_sub, \
                funM_kg_1, funM_kg_2, funM_kg_3, funM_kg_1_sub, funM_kg_2_sub, funM_kg_3_sub = \
                    ion_solver(self.chem_type, self.start_date, sim_days, presence, env, climate, chemParams, bgConc, release)
                funC_df_list = [funC_kg_1, funC_kg_1_sub, funC_kg_2, funC_kg_2_sub, funC_kg_3, funC_kg_3_sub]
                funM_df_list = [funM_kg_1, funM_kg_1_sub, funM_kg_2, funM_kg_2_sub, funM_kg_3, funM_kg_3_sub]

            elif self.chem_type == 'Nanomaterial':
                # run option 1 is for a single run
                date_array, process_array, funC_kg, funC_kg_sub, funM_kg, funM_kg_sub, \
                funC_kg_1, funC_kg_2, funC_kg_3, funM_kg_1, funM_kg_2, funM_kg_3 = \
                    nano_solver(self.start_date, time, presence, env, climate, chemParams, bgConc, release)
                funC_df_list = [funC_kg_1, funC_kg_2, funC_kg_3]
                funM_df_list = [funM_kg_1, funM_kg_2, funM_kg_3]

        # generate results and plots
        result = GenerateResult()
        result.store_output(self.chem_type, chemParams['name'], env['name'], release_scenario, release,
                            date_array, process_array, V_bulk_list, funC_df_list, funM_df_list,
                            self.output_file_path, self.file_name)


