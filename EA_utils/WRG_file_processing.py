def change_height_values_in_WRG_file(input_file_paths, hub_height):
    #usage example:
    # import EA_utils.WRG_file_processing as WRGfp
    #####################################################
    # input_file_paths = [
    #     r"C:\Support\WindFarmerAnalyst\Projects\Celtic_Sea\20241111_Tonn_Nua_rev4_layouts_and_EYAs\inputs\wind_data\maps\Vortex\VortexMap697473_160m_12WDSectors.wrg",
    #     r"C:\Support\WindFarmerAnalyst\Projects\Celtic_Sea\20241111_Tonn_Nua_rev4_layouts_and_EYAs\inputs\wind_data\maps\Vortex\at_timeseries_locations\VortexMap697473_160m_12WDSectors_at_mast.wrg"]
    # hub_height = 155
    #####################################################
    # WRGfp.change_height_values_in_WRG_file(input_file_paths, hub_height)

    import pandas as pd
    from pathlib import Path

    for idx, input_file_path in enumerate(input_file_paths):
        output_file_path = f'{input_file_path[:-4] }_at_hub_height.wrg'

        wrg_first_line = pd.read_csv(input_file_path, sep=" ", skipinitialspace=True, header=None, nrows=1).to_string(header=0, index=0, float_format='%.0f')
        map_df = pd.read_csv(input_file_path, sep=" ", skiprows=[0], skipinitialspace=True, header=None)
        map_df[4] = hub_height
        
        print("Writing to: " + output_file_path)
        with open(output_file_path, "w") as f:
            f.write(wrg_first_line + "\n")
            f.write(map_df.to_string(header=0, index=0))
    
    
def generate_WRG_at_mast_location_from_Map_WRG(input_file_directory_path, input_filenames, output_file_directory_path, mast_easting, mast_northing):
    #usage example:
    # import EA_utils.WRG_file_processing as WRGfp
    #####################################################
    # """The input file must be tab-delimited!"""
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Belgium\inputs\wind_data\maps\Vortex"
    # input_filenames = ["vortex.*.12_*.wrg"] #Takes wildcards for processing multiple files in one go
    # output_file_directory_path = input_file_directory_path + r"\at_timeseries_locations"
    # mast_easting = 469311
    # mast_northing = 5718990
    #####################################################
    # WRGfp.generate_WRG_at_mast_location_from_Map_WRG(input_file_directory_path, input_filenames, output_file_directory_path, mast_easting, mast_northing)
    
    import pandas as pd
    import numpy as np
    import os
    import glob
    from pathlib import Path
    
    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    for input_filename in input_filenames:
        input_file_paths = glob.glob(fr'{input_file_directory_path}\{input_filename}')

        for _, input_file_path in enumerate(input_file_paths):
            output_file_path = fr'{output_file_directory_path}/{Path(input_file_path).stem}_at_mast.wrg'

            map_df = pd.read_csv(input_file_path, sep=" ", skiprows=[0], skipinitialspace=True, header=None)
            map_eastings = map_df.iloc[:,1].values
            map_northings = map_df.iloc[:,2].values

            distances = ((map_eastings-mast_easting)**2+(map_northings-mast_northing)**2)**(1/2)
            index_of_first_occurrence_of_minimum_distance = np.argmin(distances)
            mast_df = map_df.iloc[index_of_first_occurrence_of_minimum_distance,:]
            mast_df_string = ' '.join(mast_df.to_string(index=False).split())
            closest_map_easting_to_mast = map_eastings[index_of_first_occurrence_of_minimum_distance]
            closest_map_northing_to_mast = map_northings[index_of_first_occurrence_of_minimum_distance]

            wrg_gridspacing = pd.read_csv(input_file_path, sep=" ", skipinitialspace=True, header=None, nrows=1, usecols=[4]).values[0][0]
            wrg_first_line = f"  1 1 {mast_easting} {mast_northing} {wrg_gridspacing}"

            print("Writing to: " + output_file_path)
            with open(output_file_path, "w") as f:
                f.write(wrg_first_line + "\n")
                f.write(mast_df_string)
    
    return output_file_path
    
def generate_WRG_from_ERA5_timeseries(site_name, data_height, input_TXT_file_path, output_WRG_file_path, minmax_easting, minmax_northing, grid_spacing, number_of_WD_sectors, trim_Windographer_TS_file):
    #usage example:
    #####################################################
    # """The input file must be tab-delimited!"""
    # site_name = "Nordsoen_I_A1"
    # data_height = 162 #m
    # input_TXT_file_path = rf"C:\Support\WindFarmerAnalyst\Projects\Denmark\inputs\wind_data\timeseries\Era5\10-year_TS\Era5_{site_name}_TS_20130101_to_20231011_{data_height}m.txt"
    # output_WRG_file_path = rf"C:\Support\WindFarmerAnalyst\Projects\Denmark\inputs\wind_data\maps\Era5\Era5_{site_name}_map_20130101_to_20231011_{data_height}m.wrg"
    # minmax_easting = [329604, 725604] #m
    # minmax_northing = [6109355, 6370355] #m
    # grid_spacing = 3000 #m
    # number_of_WD_sectors = 12
    # trim_Windographer_TS_file = True
    #####################################################
    import pandas as pd
    import numpy as np
    from scipy.stats import weibull_min
    import os
    from EA_utils.timeseries_data_extraction_generation_processing import trim_Windographer_TS_file as trim_TS_file
    
    if trim_Windographer_TS_file:
        input_TXT_file_path = trim_TS_file(site_name,data_height,input_TXT_file_path,12)

    eastings = np.arange(minmax_easting[0],minmax_easting[1]+grid_spacing,grid_spacing)
    northings = np.arange(minmax_northing[0],minmax_northing[1]+grid_spacing,grid_spacing)

    output_top_line = f"  {len(eastings)} {len(northings)} {minmax_easting[0]} {minmax_northing[0]} {grid_spacing}"

    WD_sector_midpoints = np.arange(0,360,360/number_of_WD_sectors)
    print(f"\nWD sector midpoints: {WD_sector_midpoints}")

    ts_df = pd.read_csv(input_TXT_file_path, sep="\t")
    M = f"Speed_{data_height}m"
    D = f"Direction_{data_height}m"

    output_column_widths=[10,10,10,8,5,5,6,15,3]
    for i in range(number_of_WD_sectors): output_column_widths.extend([4,4,5])

    [Weibull_k_all,_,Weibull_A_all] = weibull_min.fit(ts_df[M].values, floc = 0)
    rho = 1.225 #density of air in standard atmosphere
    avgP = np.mean(0.5*rho*ts_df[M]**3) #average power density in Wm-2

    ts_df_WD_sector_idxs = [[]]*len(ts_df)
    for df_i, WD in enumerate(ts_df[D].values):
        if WD<=WD_sector_midpoints[-1]+360/number_of_WD_sectors/2:
            ts_df_WD_sector_idxs[df_i] = (np.abs(WD_sector_midpoints - WD)).argmin()
        else:
            ts_df_WD_sector_idxs[df_i] = 0

    freq_of_occurrence = [[]]*number_of_WD_sectors
    Weibull_A_inWDsector = [[]]*number_of_WD_sectors
    Weibull_k_inWDsector = [[]]*number_of_WD_sectors
    for WD_sector_idx in range(number_of_WD_sectors):
        freq_of_occurrence[WD_sector_idx] = ts_df_WD_sector_idxs.count(WD_sector_idx)/len(ts_df_WD_sector_idxs)
        idx = [i for i,ts_df_WD_sector_idx in enumerate(ts_df_WD_sector_idxs) if ts_df_WD_sector_idx==WD_sector_idx]
        [Weibull_k_inWDsector[WD_sector_idx],_,Weibull_A_inWDsector[WD_sector_idx]] = weibull_min.fit(ts_df[M].values[idx], floc = 0)

    with open(output_WRG_file_path, 'w') as file:
        file.write(f'{output_top_line}\n')
        for northing in northings:
            items = [[]]*len(output_column_widths)
            for easting in eastings:
                items[0] = "GridPoint "
                items[1] = round(easting)
                items[2] = round(northing)
                items[3] = 0 #elevation of site - safe to call this 0 for offshore projects
                items[4] = round(data_height)
                items[5] = round(Weibull_A_all,1)
                items[6] = round(Weibull_k_all,2)
                items[7] = int(avgP)
                items[8] = round(number_of_WD_sectors)
                for WD_sector_idx in range(number_of_WD_sectors):
                    items[9+3*WD_sector_idx] = round(freq_of_occurrence[WD_sector_idx]*1000)
                    items[10+3*WD_sector_idx] = round(Weibull_A_inWDsector[WD_sector_idx]*10)
                    items[11+3*WD_sector_idx] = round(Weibull_k_inWDsector[WD_sector_idx]*100)
                output_line = "".join("%*s" % item for item in zip(tuple(output_column_widths), tuple(items)))
                file.write(f'{output_line}\n')


def convert_wrg_file_from_one_UTM_zone_to_another(file_directory_path, input_wrg_filename, output_wrg_filename, input_UTM_zone_number, output_UTM_zone_number):
    
    #!!!!!!!!!!!!DO NOT USE - WINDFARMER DOES NOT ACCEPT THE WRG FILES GENERATED BY THIS SCRIPT AS THEY ARE NOT ON A STANDARD GRID!!!!!!!!!!!!!
    
    #usage example:
    ###########################################################################
    # file_directory_path = r"C:\Support\ALP_V5_bugfix\Projects\Denmark\inputs\wind_data"
    # input_wrg_filename = r"DenmarkEastGrid_vortex_670615_152m_UTM33N.wrg"
    # output_wrg_filename = r"DenmarkEastGrid_vortex_670615_152m_UTM32N.wrg"
    # input_UTM_zone_number = 33
    # output_UTM_zone_number = 32
    ###########################################################################
    import pandas as pd
    from EA_utils.coordinate_translation import utm_to_latlon
    from EA_utils.coordinate_translation import latlon_to_utm
    
    input_wrg_file_path = fr"{file_directory_path}\{input_wrg_filename}"
    output_wrg_file_path = fr"{file_directory_path}\{output_wrg_filename}"

    df_line1 = pd.read_csv(input_wrg_file_path, delimiter=r"\s+", header=None, nrows=1)

    df = pd.read_csv(input_wrg_file_path, delimiter=r"\s+", header=None, skiprows=1)
    X_in = df[1]
    Y_in = df[2]

    lat, lon = utm_to_latlon(X_in, Y_in, input_UTM_zone_number, "N")
    X_out, Y_out, _, _ = latlon_to_utm(lat, lon, output_UTM_zone_number)
    
    df_line1[2] = int(min(X_out))
    df_line1[3] = int(min(Y_out))

    df[1] = X_out.astype(int)
    df[2] = Y_out.astype(int)

    print(rf"Saving to {output_wrg_file_path}")
    df.to_csv(output_wrg_file_path, sep=' ', index=False, header=False)

    with open(output_wrg_file_path, 'r') as original: data = original.read()
    with open(output_wrg_file_path, 'w') as modified: modified.write(f"{df_line1.to_string(header=False,index=False,index_names=False)}\n" + data)