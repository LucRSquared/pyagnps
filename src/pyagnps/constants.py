# Defining default variables and dictionaries particularly for AnnAGNPS output files
DEFAULT_EI_NUMBER = 100 # Value obtained doing a weighted average of values for the east side of the country


#--- AnnAGNPS Control Files ---

DEFAULT_OUTPUT_OPTIONS_AA = {
    'AA_Gullies_Erosion': '',
    'AA_N_Ld_Mass': '',
    'AA_N_Ld_Ratio': '',
    'AA_N_Ld_UA': '',
    'AA_N_Yld_Mass': '',
    'AA_N_Yld_Ratio': '',
    'AA_N_Yld_UA': '',
    'AA_OC_Ld_Mass': '',
    'AA_OC_Ld_Ratio': '',
    'AA_OC_Ld_UA': '',
    'AA_OC_Yld_Mass': '',
    'AA_OC_Yld_Ratio': '',
    'AA_OC_Yld_UA': '',
    'AA_P_Ld_Mass': '',
    'AA_P_Ld_Ratio': '',
    'AA_P_Ld_UA': '',
    'AA_P_Yld_Mass': '',
    'AA_P_Yld_Ratio': '',
    'AA_P_Yld_UA': '',
    'AA_Sed_Eros_Mass': '',
    'AA_Sed_Eros_Ratio': '',
    'AA_Sed_Eros_UA': 'T',
    'AA_Sed_Ld_Mass': '',
    'AA_Sed_Ld_Ratio': '',
    'AA_Sed_Ld_UA': '',
    'AA_Sed_Yld_Mass': '',
    'AA_Sed_Yld_Ratio': '',
    'AA_Sed_Yld_UA': 'T',
    'AA_Wtr_Ld_Mass': '',
    'AA_Wtr_Ld_Ratio': '',
    'AA_Wtr_Ld_UA': '',
    'AA_Wtr_Yld_Mass': '',
    'AA_Wtr_Yld_Ratio': '',
    'AA_Wtr_Yld_UA': 'T'
}

DEFAULT_OUTPUT_OPTIONS_GLOBAL = {
    'Glbl_All_V3_csv': '',
    'Glbl_All_V3_dpp': '',
    'Glbl_All_V3_npt': '',
    'Glbl_All_V3_sim': '',
    'Glbl_All_V3_txt': '',
    'Log_to_File': '',
    'Log_to_Screen': '',
    'Warning_File': '',
    'V1/2_Output_Files': '',
    'Glbl_All_Cells': '',
    'Glbl_All_Feedlots': '',
    'Glbl_All_Fld_Ponds': '',
    'Glbl_All_Gullies': '',
    'Glbl_All_Pt_Srcs': '',
    'Glbl_All_Reaches': '',
    'Glbl_All_Impound': '',
    'Glbl_All_Wetlands': '',
    'Glbl_All_AA_Nutr': '',
    'Glbl_All_AA_Pest': '',
    'Glbl_All_AA_Sed': '',
    'Glbl_All_AA_Wtr': '',
    'Glbl_All_EV_Nutr': 'F',
    'Glbl_All_EV_Pest': 'F',
    'Glbl_All_EV_Sed': 'T',
    'Glbl_All_EV_Wtr': 'T',
    'Glbl_All_V2/3_Mass': '',
    'Glbl_All_V2/3_Ratio': '',
    'Glbl_All_V2/3_UA': '',
    'V2_Concepts': 'T',
    'V2_AA': '',
    'V2_EV': '',
    'V1_AA': '',
    'V1_EV': 'T'
 }

# As given by Ron
DEFAULT_OUTPUT_OPTIONS_TBL = {
    'CCHE1D': 'T',
    'CONCEPTS_XML': '',
    'Gaging_Station_Hyd': 'F',
    'REMM': '',
    'Gaging_Station_Evt': ''
}

# For AIMS this is what I think we should do with just the outlet reach
# DEFAULT_OUTPUT_OPTIONS_TBL = {
#     'CCHE1D': 'F',
#     'CONCEPTS_XML': '',
#     'Gaging_Station_Hyd': 'T',
#     'REMM': '',
#     'Gaging_Station_Evt': ''
# }

DEFAULT_ANNAGNPS_ID = {
    'Version': 6.0,
    'Input_Units': 1,
    'Output_Units': 1,
    'CCHE1D_Output_Units': float('nan'),
    'Screen_Output_Units': float('nan')
}

DEFAULT_GLOBAL_FACTORS_FLAGS = {
    'Hdct_Detachment_Coef_a': '',
    'Hdct_Detachment_Exp_Coef_b': '',
    'Urban_Repair_Month': '',
    'Urban_Repair_Day': '',
    'Urban_Repair_Year': '',
    'Cropland_Repair_Month': '',
    'Cropland_Repair_Day': '',
    'Cropland_Repair_Year': '',
    'Forest_Repair_Month': '',
    'Forest_Repair_Day': '',
    'Forest_Repair_Year': '',
    'Pasture_Repair_Month': '',
    'Pasture_Repair_Day': '',
    'Pasture_Repair_Year': '',
    'Rangeland_Repair_Month': '',
    'Rangeland_Repair_Day': '',
    'Rangeland_Repair_Year': '',
    'Hdct_Erodibility_Coef_a': '',
    'Hdct_Erodibility_Exp_Coef_b': '',
    'Width_Nachtergaele': '',
    'Width_Hydraulic_Geometry': '',
    'Width_Non-submerging_Tailwater': '',
    'Width_Woodwards_Equilibrium': '',
    'Width_Woodwards_Ultimate': '',
    'Width_Wells_Eq.9': '',
    'Erosion_Vrfy': '',
    'Hydrograph_Vrfy': '',
    'Nickpoint_Vrfy': '',
    'Repair_Dates_Vrfy': '',
    'Sed_Yield_to_Gully_Mouth_Vrfy': '',
    'Sed_Yield_to_Rcvg_Reach_Vrfy': '',
    'Min_Interception_Evaporation': '',
    'Max_Interception_Evaporation': '',
    'Detention_Coef_a': '',
    'Detention_Coef_b': '',
    'RCN_Convergence_Tolerance': '',
    'RCN_Max_Iterations': '',
    'Avbl_Soil_Moist_Ratio_AMC_II': '',
    'Max_Avbl_Sed_Conc_for_Sht_Flw': '',
    'Max_Avbl_Sed_Conc_for_Conc_Flw': '',
    'AA_Unit_Area_Baseflow': '',
    'RCN_Calib_Only': '',
    'Calculate_Baseflow': '',
    'FAO_ET_Enhancement': '',
    'Basal_Crop_Coef_Climate_Adjust': '',
    'Wshd_Storm_Type_ID': 'Std. SCS Type II',
    'Dflt_Geology_ID': '',
    'Dflt_Hydraulic_Geom_ID': '',
    'Dflt_Init_Soil_Conditions_ID': '',
    'Dflt_Crop_RCN_ID': '',
    'Dflt_Non-Crop_RCN_ID': '',
    'Width_Wells_Eq.8': '',
    'Width_Reserved_i': '',
    'Width_Reserved_j': '',
    'Width_Reserved_k': '',
    'Critical_Shear_Stress': '',
    'RUSLE2_Flag': '',
    'Dflt_RUSLE2_ID': '',
    'Reach_Routing_Flag':'',
    'Input_Units_Code': 0
}

DEFAULT_SIM_PERIOD_DATA = {
    'Simulation_Begin_Month': '',
    'Simulation_Begin_Day': '',
    'Simulation_Begin_Year': '',
    'Simulation_End_Month': '',
    'Simulation_End_Day': '',
    'Simulation_End_Year': '',
    'Rainfall_Fctr': '',
    '10-Year_EI': '',
    'EI_Number': DEFAULT_EI_NUMBER,
    'Irrigation_Climate_Code': '',
    'Soil_Moisture_Steps': float('nan'),
    'Annual_K_Fctr_Code': '',
    'Variable_K_Fctr_Code': '',
    'Number_Init_Years': 2,
    'Init_Method_Code': float('nan'),
    'Winter_Bouts': float('nan'),
    'Input_Units_Code': 0
}


# NLDAS-2

# GRIB format
_NLDAS_PRODUCTS_002 = { 'NLDAS_FORA0125_H.002':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_FORA0125_H.A',
                             'description': 'NLDAS-2 hourly primary forcing dataset, File A'},
                        'NLDAS_FORA0125_M.002':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_FORA0125_M.A',
                             'description': 'NLDAS-2 monthly primary forcing dataset, File A'},
                        'NLDAS_FORA0125_MC.002':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_FORA0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology primary forcing dataset, File A (30y averaged, 1980-2009)'},
                        'NLDAS_FORB0125_H.002':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_FORB0125_H.A',
                             'description': 'NLDAS-2 hourly secondary forcing dataset, File B'},
                        'NLDAS_FORB0125_M.002':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_FORB0125_M.A',
                             'description': 'NLDAS-2 monthly secondary forcing dataset, File B'},
                        'NLDAS_FORB0125_MC.002':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_FORB0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology secondary forcing dataset, File B (30y averaged, 1980-2009)'},
                        'NLDAS_MOS0125_H.002':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_MOS0125_H.A',
                             'description': 'NLDAS-2 hourly MOSAIC model data'},
                        'NLDAS_MOS0125_M.002':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_MOS0125_M.A',
                             'description': 'NLDAS-2 monthly MOSAIC model data'},
                        'NLDAS_MOS0125_MC.002':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_MOS0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology MOSAIC model dataset (30y averaged, 1980-2009)'},
                        'NLDAS_NOAH0125_H.002':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_NOAH0125_H.A',
                             'description': 'NLDAS-2 hourly NOAH model data'},
                        'NLDAS_NOAH0125_M.002':
                            {'type': 'monthly',
                              'fileroot': 'NLDAS_NOAH0125_M.A',
                                'description': 'NLDAS-2 monthly NOAH model data'},
                        'NLDAS_NOAH0125_MC.002':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_NOAH0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology NOAH model dataset (30y averaged, 1980-2009)'},
                        'NLDAS_VIC0125_H.002':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_VIC0125_H.A',
                             'description': 'NLDAS-2 hourly VIC model data'},
                        'NLDAS_VIC0125_M.002':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_VIC0125_M.A',
                             'description': 'NLDAS-2 monthly VIC model data'},
                        'NLDAS_VIC0125_MC.002':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_VIC0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology VIC model dataset (30y averaged, 1980-2009)'},
}

# NetCDF format
_NLDAS_PRODUCTS_20 = { 'NLDAS_FORA0125_H.2.0':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_FORA0125_H.A',
                             'description': 'NLDAS-2 hourly primary forcing dataset, File A'},
                        'NLDAS_FORA0125_D.2.0': # Non-NLDAS product, this one was aggregated in hourse
                            {'type': 'daily',
                             'fileroot': 'NLDAS_FORA0125_D.A',
                             'description': 'NLDAS-2 daily primary forcing dataset, File A'},
                        'NLDAS_FORA0125_M.2.0':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_FORA0125_M.A',
                             'description': 'NLDAS-2 monthly primary forcing dataset, File A'},
                        'NLDAS_FORA0125_MC.2.0':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_FORA0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology primary forcing dataset, File A (30y averaged, 1980-2009)'},
                        'NLDAS_FORB0125_H.2.0':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_FORB0125_H.A',
                             'description': 'NLDAS-2 hourly secondary forcing dataset, File B'},
                        'NLDAS_FORB0125_M.2.0':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_FORB0125_M.A',
                             'description': 'NLDAS-2 monthly secondary forcing dataset, File B'},
                        'NLDAS_FORB0125_MC.2.0':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_FORB0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology secondary forcing dataset, File B (30y averaged, 1980-2009)'},
                        'NLDAS_MOS0125_H.2.0':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_MOS0125_H.A',
                             'description': 'NLDAS-2 hourly MOSAIC model data'},
                        'NLDAS_MOS0125_M.2.0':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_MOS0125_M.A',
                             'description': 'NLDAS-2 monthly MOSAIC model data'},
                        'NLDAS_MOS0125_MC.2.0':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_MOS0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology MOSAIC model dataset (30y averaged, 1980-2009)'},
                        'NLDAS_NOAH0125_H.2.0':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_NOAH0125_H.A',
                             'description': 'NLDAS-2 hourly NOAH model data'},
                        'NLDAS_NOAH0125_M.2.0':
                            {'type': 'monthly',
                              'fileroot': 'NLDAS_NOAH0125_M.A',
                                'description': 'NLDAS-2 monthly NOAH model data'},
                        'NLDAS_NOAH0125_MC.2.0':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_NOAH0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology NOAH model dataset (30y averaged, 1980-2009)'},
                        'NLDAS_VIC0125_H.2.0':
                            {'type': 'hourly',
                             'fileroot': 'NLDAS_VIC0125_H.A',
                             'description': 'NLDAS-2 hourly VIC model data'},
                        'NLDAS_VIC0125_M.2.0':
                            {'type': 'monthly',
                             'fileroot': 'NLDAS_VIC0125_M.A',
                             'description': 'NLDAS-2 monthly VIC model data'},
                        'NLDAS_VIC0125_MC.2.0':
                            {'type': 'monthly_climatology',
                             'fileroot': 'NLDAS_VIC0125_MC.ACLIM',
                             'description': 'NLDAS-2 monthly climatology VIC model dataset (30y averaged, 1980-2009)'},
}

_BASE_URL_NLDAS = 'https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/'
