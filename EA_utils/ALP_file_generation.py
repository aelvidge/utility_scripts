def generate_ALP_area_file_from_ArcGIS_generated_site_coord_txt_file(site_name, site_name_varname, coordinate_system_name, x_coord_varname, y_coord_varname, input_file_path, output_file_directory_path, buildable_area_or_exclusion_zone, text_to_append_to_filename):
    #usage example:
    # import EA_utils.ALP_file_generation as ALPfg
    #####################################################
    # site_name = "N9p3_v1" #If assigned, then site_name_varname (below) is ignored. If unassigned (left empty), then the site name is extracted from site_name_varname
    # site_name_varname = "" #If assigned, then the site name is extracted from this variable, and site_name (above) is ignored. If unassigned (left empty), then site_name is used
    # coordinate_system_name = "UTM32N"
    # x_coord_varname = f"X_{coordinate_system_name}"
    # y_coord_varname = f"Y_{coordinate_system_name}"
    # """The input file must be tab-delimited!"""
    # input_file_path = fr"C:\Support\ArcGIS\Projects\Germany\outputs\target_site_buildable_areas\V1\{site_name}_buildable_area_vertex_coords_UTM32N.txt"
    # """Best to have all directory paths within C:\Support\ALP_V5_bugfix\Projects\ as this would be a sensible path for others to us on their own computers (whereas a path that has "AE78333" in it would not be)"""
    # output_file_directory_path = fr"C:\Support\ALP_V5_bugfix\Projects\Germany\inputs\site_buildable_area_coordinates\V1"
    # buildable_area_or_exclusion_zone = 'buildable_area'
    # text_to_append_to_filename = '_v1'
    #####################################################
    # ALPfg.ALP_area_file_from_ArcGIS_generated_site_coord_txt_file(site_name, site_name_varname, coordinate_system_name, x_coord_varname, y_coord_varname, input_file_path, output_file_directory_path, buildable_area_or_exclusion_zone, text_to_append_to_filename)
    
    import pandas as pd
    import os
    from unidecode import unidecode

    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    df = pd.read_csv(input_file_path,sep="\t")

    if site_name_varname:
        site_names = df.loc[:,site_name_varname].unique()
        site_names_to_output = [unidecode(site_name) for site_name in site_names]
        site_names_to_output = [site_name.replace("- ", "") for site_name in site_names_to_output]
        site_names_to_output = [site_name.replace(" ", "_") for site_name in site_names_to_output]
        print("\nInput site names: " + str(list(site_names)))
        print("Output site names: " + str(site_names_to_output) +"\n")
    else:
        site_names = [site_name]
        site_names_to_output = [site_name]

    for site_name_i,site_name in enumerate(site_names):
        site_name_to_output = site_names_to_output[site_name_i]
        if buildable_area_or_exclusion_zone == 'buildable_area':
            output_file_path = output_file_directory_path + "\\" + site_name_to_output + "_" + coordinate_system_name + "_ALPboundaryfile_BuildableArea" + text_to_append_to_filename + ".txt"
        elif buildable_area_or_exclusion_zone == 'exclusion_zone':
            output_file_path = output_file_directory_path + "\\" + site_name_to_output + "_" + coordinate_system_name + "_ALPboundaryfile_ExclusionZone" + text_to_append_to_filename + ".txt"
        if site_name_varname:
            x_data = df.loc[df[site_name_varname] == site_name,x_coord_varname]
            y_data = df.loc[df[site_name_varname] == site_name,y_coord_varname]
        else:
            x_data = df.loc[:,x_coord_varname]
            y_data = df.loc[:,y_coord_varname]
        print("Writing to: " + output_file_path)
        with open(output_file_path, "w") as f:
            f.write(str(len(x_data)) + "\t" + str(len(x_data)))
            for i in range(2): #ALP requires coordinates to be repeated
                for point_i,_ in enumerate(x_data):
                    f.write("\n" + str(x_data.iloc[point_i]) + "\t" + str(y_data.iloc[point_i]))

def ALP_neighbour_array_file_from_ArcGIS_generated_turbine_coord_txt_file(region_name, input_file_directory_path, input_filenames, turbine_file_directory_path, turbine_filenames, output_file_directory_path, coordinate_system_name, x_coord_varname, y_coord_varname):
    #usage example:
    #####################################################
    # region_name = "Denmark"
    # """Best to have all directory paths within C:\Support\ALP_V5_bugfix\Projects\ as this would be a sensible path for others to us on their own computers (whereas a path that has "AE78333" in it would not be)"""
    # input_file_directory_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\Denmark\Data\ArcGIS\output"
    # input_filenames = ["WindFarms20231001_KriegersFlak2_withbuffer_turbine_layout.txt",
    #                 "WindFarms20231001_Thor_withbuffer_turbine_layout.txt",
    #                 "WindFarms20231001_VesterhavNordSud_withbuffer_turbine_layout.txt"]
    # turbine_file_directory_path = r"C:\Support\ALP_V5_bugfix\Projects\Denmark\inputs\turbine_files\neighbour_sites"
    # turbine_filenames = ["KriegersFlak2_20210122_DNV_Arch_15MW230.pow","Thor_20210504_SG_14MW-222_HWRT+LEP_AD1.230.pow","VesterhavNordSud_20180403_Siemens_SWT-8.0-167_General_PB+HWRT.pow"]
    # output_file_directory_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\Denmark\Data\ArcGIS\output"
    # coordinate_system_name = "UTM32N"
    # x_coord_varname = "X_UTM32N"
    # y_coord_varname = "Y_UTM32N"
    #####################################################
    import pandas as pd
    import os

    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    input_file_paths = [input_file_directory_path + "\\" + input_filename for input_filename in input_filenames]
    turbine_file_paths = [turbine_file_directory_path + "\\" + turbine_filename for turbine_filename in turbine_filenames]

    output_file_path = output_file_directory_path + "\\" + region_name + "_" + coordinate_system_name + "_ALPneighbourturbines_rough_layouts.csv"
    print("Writing to: " + output_file_path)
    with open(output_file_path, "w") as f:
        for input_file_path_i,input_file_path in enumerate(input_file_paths):
            df = pd.read_csv(input_file_path,sep="\t")
            if input_file_path_i>0:
                f.write("\n")
            f.write("Turbine file," + turbine_file_paths[input_file_path_i])
            f.write("\nTurbine ID,Eastings (m),Northings (m)")
            for turbine_i, turbine_data in df.iterrows():
                f.write("\n" + str(turbine_i+1) + "," + str(round(turbine_data[x_coord_varname])) + "," + str(round(turbine_data[y_coord_varname])))
                
def ALP_neighbour_array_file_from_ArcGIS_generated_turbine_coord_xlsx_file(region_name, input_file_path, turbine_file_directory_path, output_file_directory_path, coordinate_system_name, site_name_varname, x_coord_varname, y_coord_varname):
    #usage example:
    #####################################################
    # region_name = "Denmark"
    # input_file_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\Denmark\Data\ArcGIS\output\20231009 Denmark Offshore WTG Coords Digitised 4c.xlsx"
    # """Best to have all directory paths within C:\Support\ALP_V5_bugfix\Projects\ as this would be a sensible path for others to us on their own computers (whereas a path that has "AE78333" in it would not be)"""
    # turbine_file_directory_path = r"C:\Support\ALP_V5_bugfix\Projects\Denmark\inputs\turbine_files\neighbour_sites"
    # output_file_directory_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\Denmark\Data\ArcGIS\output"
    # coordinate_system_name = "UTM32N"
    # site_name_varname = "Project"
    # x_coord_varname = "Easting (m) UTM32N"
    # y_coord_varname = "Northing (m) UTM32N"
    #####################################################
    import pandas as pd
    import os
    from unidecode import unidecode
    
    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    df = pd.read_excel(input_file_path)
    print(df)
    site_names = df.loc[:,site_name_varname].dropna().unique()
    site_names_to_output = [unidecode(site_name) for site_name in site_names]
    site_names_to_output = [site_name.replace("- ", "") for site_name in site_names_to_output]
    site_names_to_output = [site_name.replace(" ", "_") for site_name in site_names_to_output]
    site_names_for_turbine_files = [site_name.replace("_", "") for site_name in site_names_to_output]

    print("\nInput site names: " + str(list(site_names)))

    output_file_path = output_file_directory_path + "\\" + region_name + "_" + coordinate_system_name + "_ALPneighbourturbines_specified_layouts.csv"
    print("Writing to: " + output_file_path)
    with open(output_file_path, "w") as f:
        for site_name_i,site_name in enumerate(site_names):
            site_name_to_output = site_names_to_output[site_name_i]
            site_name_for_turbine_file = site_names_for_turbine_files[site_name_i]
            turbine_file_name_list = [fname for fname in os.listdir(turbine_file_directory_path) if site_name_for_turbine_file + "_" in fname]
            if len(turbine_file_name_list)>1:
                raise ValueError('More than one turbine file matches the site name ' + site_name_for_turbine_file)
            elif not turbine_file_name_list[0]:
                raise ValueError('A turbine file for this site does not exist in the specified directory')
            else:
                turbine_file_name = turbine_file_name_list[0]
                print("Turbine file used: " + turbine_file_name)
            turbine_file_path = turbine_file_directory_path + "\\" + turbine_file_name
            if site_name_i>0:
                f.write("\n")
            f.write("Turbine file," + turbine_file_path)
            f.write("\nTurbine ID,Eastings (m),Northings (m)")
            x_data = round(df.loc[df[site_name_varname] == site_name,x_coord_varname]).astype(int)
            y_data = round(df.loc[df[site_name_varname] == site_name,y_coord_varname]).astype(int)
            for point_i,_ in enumerate(x_data):
                f.write("\n" + str(point_i+1) + "," + str(x_data.iloc[point_i]) + "," + str(y_data.iloc[point_i]))