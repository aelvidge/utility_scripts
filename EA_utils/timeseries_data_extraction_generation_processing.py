def prepare_VortexTS_files_for_import_into_Windographer(input_file_directory_path, output_file_directory_path, process_files_that_include_this_text_in_filename, append_to_filenames, number_of_lines_at_top_of_each_file_to_delete):
    #usage example:
    # import EA_utils.timeseries_data_extraction_generation_processing as TS_processing
    #####################################################
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Germany\inputs\wind_data\timeseries\Vortex"
    # output_file_directory_path = input_file_directory_path
    # process_files_that_include_this_text_in_filename = 'vortex.serie' #file extension must be omitted
    # append_to_filenames = ""
    # number_of_lines_at_top_of_each_file_to_delete = 3 # use 3 as default
    #####################################################
    # TS_processing.prepare_VortexTS_files_for_import_into_Windographer(input_file_directory_path, output_file_directory_path, process_files_that_include_this_text_in_filename, append_to_filenames, number_of_lines_at_top_of_each_file_to_delete)
    import os
    import glob
    
    if process_files_that_include_this_text_in_filename[-4:]=='.txt':
        raise ValueError('File extionsion must be omitted from process_files_that_include_this_text_in_filename')
    
    input_filenames = [os.path.basename(x) for x in glob.glob(input_file_directory_path + '\\*' + process_files_that_include_this_text_in_filename + '*.txt')]
    if not input_filenames:
        print(input_file_directory_path + '\\*' + process_files_that_include_this_text_in_filename + '*.txt')
        raise ValueError('No matching files')
    
    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    output_file_path = [[]]*len(input_filenames)
    for input_filename_i,input_filename in enumerate(input_filenames):
        input_file_path = input_file_directory_path + "\\" + input_filename
        output_file_name = input_filename[:-4] + append_to_filenames
        output_file_name = output_file_name.replace(".","_")
        output_file_name = output_file_name.replace(" ","_")
        output_file_name = output_file_name.replace("+","plus")
        output_file_name = output_file_name + ".txt"
        output_file_path[input_filename_i] = output_file_directory_path + "\\" + output_file_name
        
        i_end = input_filename.find("m UTC")+1
        i_start = input_filename[0:i_end].rfind(" ")+1
        height_str = input_filename[i_start:i_end]
        
        with open(input_file_path, 'r') as fin:
            data = fin.read().splitlines(True)
            data = data[number_of_lines_at_top_of_each_file_to_delete:]
        if data[0][:-1] == "YYYYMMDD HHMM  M(m/s) D(deg)  T(C)  De(k/m3) PRE(hPa)      RiNumber  RH(%)   RMOL(1/m)":
            data[0] = "YYYYMMDD HHMM  M_" + height_str + "(m/s) D_" + height_str + "(deg)  T_" + height_str + "(C)  De_" + height_str + "(k/m3) PRE_" + height_str + "(hPa)      RiNumber_" + height_str + "  RH_" + height_str + "(%)   RMOL_" + height_str + "(1/m)\n"
        elif data[0][:-1] == "YYYYMMDD HHMM  M(m/s) D(deg)  T(C)  De(k/m3) PRE(hPa)      RiNumber  RH(%)":
            data[0] = "YYYYMMDD HHMM  M_" + height_str + "(m/s) D_" + height_str + "(deg)  T_" + height_str + "(C)  De_" + height_str + "(k/m3) PRE_" + height_str + "(hPa)      RiNumber_" + height_str + "  RH_" + height_str + "(%)\n"
        else:
            raise ValueError('The header line (line ' + str(number_of_lines_at_top_of_each_file_to_delete+1) + ' in the original file) is not as expected - check the input file is correct')
        
        with open(output_file_path[input_filename_i], 'w') as fout:
            print("Writing to: " + output_file_path[input_filename_i])
            fout.writelines(data)
    
    return output_file_path

def trim_Windographer_TS_file(site_name,data_height,input_file_path,number_of_lines_at_top_of_each_file_to_delete):
    #usage example:
    #####################################################
    # site_name = "DanishKriegersFlakII"
    # data_height = 78 #in m
    # input_file_path = rf"C:\Support\WindFarmerAnalyst\Projects\Denmark\inputs\wind_data\timeseries\Era5\10-year_TS\Era5_{site_name}_TS_20130101_to_20231011_{data_height}m.txt"
    # number_of_lines_at_top_of_each_file_to_delete = 12
    #####################################################
    from pathlib import Path
    
    header_line = f"DateTime	Speed_{data_height}m	Direction_{data_height}m\n"
    output_file_path = rf"{Path(input_file_path).parent}\{Path(input_file_path).stem}_trimmed.txt"

    with open(input_file_path, 'r') as fin:
        data = fin.read().splitlines(True)
        data = data[number_of_lines_at_top_of_each_file_to_delete:]
        data[0] = header_line
    print("Writing to: " + output_file_path)
    with open(output_file_path, 'w') as fout:
        fout.writelines(data)

    return output_file_path


def ts_dataset_from_era5_data_from_CDS_API(input_filepath, output_directory, U10m_varname, V10m_varname, output_filename_10m, U100m_varname, V100m_varname, output_filename_100m):
    #usage example:
    # import EA_utils.timeseries_data_extraction_generation_processing as ts_processing
    # #############################################################################################
    # input_filepath = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\R&D\tool_development\Climate Intelligence Tool\trial data - Creyab Windfarm\ERA-5 data\Era5_Creyab_Windfarm_20000101_to_20231231.nc"
    # output_directory = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\R&D\tool_development\Climate Intelligence Tool\trial data - Creyab Windfarm\ERA-5 data"

    # U10m_varname = "u10" #set to empty string if no 10-m wind variable
    # V10m_varname = "v10" #set to empty string if no 10-m wind variable
    # output_filename_10m = "Era5_Creyab_Windfarm_20000101_to_20231231.txt" #set to empty string if no 10-m wind variable

    # U100m_varname = "u100" #set to empty string if no 100-m wind variable
    # V100m_varname = "v100" #set to empty string if no 100-m wind variable
    # output_filename_100m = "Era5_Creyab_Windfarm_20000101_to_20231231.txt" #set to empty string if no 100-m wind variable
    # #############################################################################################
    # ts_processing.ts_dataset_from_era5_data_from_CDS_API(input_filepath, output_directory, U10m_varname, V10m_varname, output_filename_10m, U100m_varname, V100m_varname, output_filename_100m)
    import xarray as xr
    import metpy.calc as mpcalc
    from metpy.units import units
    import pandas as pd
    import numpy as np
    
    ds = xr.open_dataset(input_filepath)

    print("\nLatitude = " + str(ds["latitude"].data[0]))
    print("Longitude = " + str(ds["longitude"].data[0]))
    times = ds["time"].data
    Y, M, D, h, m, s = [times.astype('datetime64[%s]' % kind) for kind in 'YMDhms']
    years = Y.astype(int) + 1970
    months = M.astype(int) % 12 + 1
    days = (D - M).astype(int) + 1
    hours = (h - D).astype(int)
    minutes = (m - h).astype(int)
    seconds = (s - m).astype(int)
    years_str = [(str(year)) for year in years]
    months_str = [("0" + str(month))[-2:] for month in months]
    days_str = [("0" + str(day))[-2:] for day in days]
    hours_str = [("0" + str(hour))[-2:] for hour in hours]
    minutes_str = [("0" + str(minute))[-2:] for minute in minutes]
    datesYYYYMMDD_str = [years_str[i] + months_str[i] + days_str[i] for i,_ in enumerate(times)]
    timesHHMM_str = [hours_str[i] + minutes_str[i] for i,_ in enumerate(times)]

    if U10m_varname!='':
        ws10m = np.squeeze(mpcalc.wind_speed(ds[U10m_varname].data * units('m/s'),ds[V10m_varname].data * units('m/s')).magnitude)
        wd10m = np.squeeze(mpcalc.wind_direction(ds[U10m_varname].data * units('m/s'),ds[V10m_varname].data * units('m/s')).magnitude)
    if U100m_varname!='':
        ws100m = np.squeeze(mpcalc.wind_speed(ds[U100m_varname].data * units('m/s'),ds[V100m_varname].data * units('m/s')).magnitude)
        wd100m = np.squeeze(mpcalc.wind_direction(ds[U100m_varname].data * units('m/s'),ds[V100m_varname].data * units('m/s')).magnitude)

    if U10m_varname!='':
        df10m_out = pd.DataFrame({"YYYYMMDD" : datesYYYYMMDD_str, 'HHMM' : timesHHMM_str, 'M_10m(m/s)' : ws10m, "D_10m(deg)" : wd10m})
    if U100m_varname!='':
        df100m_out = pd.DataFrame({"YYYYMMDD" : datesYYYYMMDD_str, 'HHMM' : timesHHMM_str, 'M_100m(m/s)' : ws100m, "D_100m(deg)" : wd100m})

    if U10m_varname!='':
        df10m_out.to_csv(output_directory + "\\" + output_filename_10m, sep='\t', index=False, float_format='%.3f')
    if U100m_varname!='':
        df100m_out.to_csv(output_directory + "\\" + output_filename_100m, sep='\t', index=False, float_format='%.3f')

    print("Saving to " + output_directory + "\\")
    if U10m_varname!='':
        with open(output_directory + "\\" + output_filename_10m, 'r') as original: data = original.read()
        with open(output_directory + "\\" + output_filename_10m, 'w') as modified: modified.write("Lat=51.857293  Lon=-7.052537  Height=10  Timezone=00.0   ASL-Height=0\n" + data)
    if U100m_varname!='':
        with open(output_directory + "\\" + output_filename_100m, 'r') as original: data = original.read()
        with open(output_directory + "\\" + output_filename_100m, 'w') as modified: modified.write("Lat=51.857293  Lon=-7.052537  Height=10  Timezone=00.0   ASL-Height=0\n" + data)
    
def ts_dataset_from_era5_data_from_CDS_toolbox(input_directory, output_directory, U10m_filename, V10m_filename, U10m_varname, V10m_varname, output_filename_10m, U100m_filename, V100m_filename, U100m_varname, V100m_varname, output_filename_100m):
    #usage example:
    ##############################################################################################
    # input_directory = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\CelticSeas\ConstrainingWindUncertainty\Data\Era5_data"
    # output_directory = input_directory

    # U10m_filename = "2004-23_singlelevs_1hrly_U10m_combined.nc"
    # V10m_filename = "2004-23_singlelevs_1hrly_V10m_combined.nc"
    # U10m_varname = "u10"
    # V10m_varname = "v10"
    # output_filename_10m = "Era5TS_1hrly_at_M5buoy_10m.txt"

    # U100m_filename = "2004-23_singlelevs_1hrly_U100m_combined.nc"
    # V100m_filename = "2004-23_singlelevs_1hrly_V100m_combined.nc"
    # U100m_varname = "u100"
    # V100m_varname = "v100"
    # output_filename_100m = "Era5TS_1hrly_at_M5buoy_100m.txt"
    ##############################################################################################
    import xarray as xr
    import metpy.calc as mpcalc
    from metpy.units import units
    import pandas as pd

    ds_U10m = xr.open_dataset(input_directory + "\\" + U10m_filename)
    ds_V10m = xr.open_dataset(input_directory + "\\" + V10m_filename)
    ds_U100m = xr.open_dataset(input_directory + "\\" + U100m_filename)
    ds_V100m = xr.open_dataset(input_directory + "\\" + V100m_filename)

    print("\nLatitude = " + str(ds_U10m["lat"].data))
    print("Longitude = " + str(ds_U10m["lon"].data))
    times = ds_U10m["time"].data
    Y, M, D, h, m, s = [times.astype('datetime64[%s]' % kind) for kind in 'YMDhms']
    years = Y.astype(int) + 1970
    months = M.astype(int) % 12 + 1
    days = (D - M).astype(int) + 1
    hours = (h - D).astype(int)
    minutes = (m - h).astype(int)
    seconds = (s - m).astype(int)
    years_str = [(str(year)) for year in years]
    months_str = [("0" + str(month))[-2:] for month in months]
    days_str = [("0" + str(day))[-2:] for day in days]
    hours_str = [("0" + str(hour))[-2:] for hour in hours]
    minutes_str = [("0" + str(minute))[-2:] for minute in minutes]
    datesYYYYMMDD_str = [years_str[i] + months_str[i] + days_str[i] for i,_ in enumerate(times)]
    timesHHMM_str = [hours_str[i] + minutes_str[i] for i,_ in enumerate(times)]

    ws10m = mpcalc.wind_speed(ds_U10m[U10m_varname].data * units('m/s'),ds_V10m[V10m_varname].data * units('m/s')).magnitude
    wd10m = mpcalc.wind_direction(ds_U10m[U10m_varname].data * units('m/s'),ds_V10m[V10m_varname].data * units('m/s')).magnitude
    ws100m = mpcalc.wind_speed(ds_U100m[U100m_varname].data * units('m/s'),ds_V100m[V100m_varname].data * units('m/s')).magnitude
    wd100m = mpcalc.wind_direction(ds_U100m[U100m_varname].data * units('m/s'),ds_V100m[V100m_varname].data * units('m/s')).magnitude

    df10m_out = pd.DataFrame({"YYYYMMDD" : datesYYYYMMDD_str, 'HHMM' : timesHHMM_str, 'M_10m(m/s)' : ws10m, "D_10m(deg)" : wd10m})
    df100m_out = pd.DataFrame({"YYYYMMDD" : datesYYYYMMDD_str, 'HHMM' : timesHHMM_str, 'M_100m(m/s)' : ws100m, "D_100m(deg)" : wd100m})

    df10m_out.to_csv(output_directory + "\\" + output_filename_10m, sep='\t', index=False, float_format='%.3f')
    df100m_out.to_csv(output_directory + "\\" + output_filename_100m, sep='\t', index=False, float_format='%.3f')

    print("Saving to " + output_directory + "\\")
    with open(output_directory + "\\" + output_filename_10m, 'r') as original: data = original.read()
    with open(output_directory + "\\" + output_filename_10m, 'w') as modified: modified.write("Lat=51.41  Lon=-6.42  Height=10  Timezone=00.0   ASL-Height=0\n" + data)
    with open(output_directory + "\\" + output_filename_100m, 'r') as original: data = original.read()
    with open(output_directory + "\\" + output_filename_100m, 'w') as modified: modified.write("Lat=51.41  Lon=-6.42  Height=10  Timezone=00.0   ASL-Height=0\n" + data)
    
def get_mean_density_from_VortexTS(input_file_path):
    #Run on file generated from prepare_VortexTS_files_for_import_into_Windographer.py
    #usage example:
    # import EA_utils.timeseries_data_extraction_generation_processing as TS_processing
    ##############################################################################################
    # input_file_path = r"C:\Support\WindFarmerAnalyst\Projects\AR9\inputs\wind_data\timeseries\Vortex\Bowdun_vortex_serie_344255_20y_175m_UTCplus00_0_ERA5.txt"
    ##############################################################################################
    # _ = TS_processing.get_mean_density_from_VortexTS(input_file_path)
    
    import numpy as np

    with open(input_file_path, 'r') as fin:
        lines=fin.readlines()[1:]
        densities = []
        for line in lines:
            split_line = line.split()
            densities.append(split_line[5])
        fin.close()
    mean_density = np.mean(np.asarray(densities).astype(float))

    print(f'Mean air density = {mean_density:.5f} kg m-3')
    return mean_density
    
def get_proportion_of_VortexTS_meeting_a_criteria(input_file_path, variables, criteria):
    #Run on file generated from prepare_VortexTS_files_for_import_into_Windographer.py
    #usage example:
    # import EA_utils.timeseries_data_extraction_generation_processing as TS_processing
    #############################################################################################
    # input_file_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\Celtic_Sea\02 Wind Data\timeseries\Models\Vortex\near_site\ERA5\vortex_serie_397209_20y_155m_UTCplus01_0_ERA5.txt"
    # variables = ["temperature","wind_speed"] # list of variables: "wind_speed", "wind_direction", "temperature", "density", "pressure", "Richardson_number", or "relative_humidity"
    # criteria = [">15",">14"] # list of criteria, each of which must be either '>X', '<X', '>=X', '<=X', or '==X', where X is a threshold value (float) of the chosen variable
    #############################################################################################
    # TS_processing.get_proportion_of_VortexTS_meeting_a_criteria(input_file_path, variable, criteria)
    import numpy as np

    var_list_in_column_order = ["date", "time", "wind_speed", "wind_direction", "temperature", "density", "pressure", "Richardson_number", "relative_humidity"]

    for variable_idx, variable in enumerate(variables):
        variable = variables[variable_idx]
        var_idx = var_list_in_column_order.index(variable)
        with open(input_file_path, 'r') as fin:
            lines=fin.readlines()[1:]
            if variable_idx == 0:
                var_values = np.empty(shape=(len(variables),len(lines)))
            for line_idx, line in enumerate(lines):
                split_line = line.split()
                var_values[variable_idx,line_idx] = split_line[var_idx]
            fin.close()

    var_values_processing = var_values
    for criterion_idx, criterion in enumerate(criteria):
        
        if criterion[1] == '=':
            var_threshold_value = float(criterion[2:])
            if criterion[0] == '>':
                var_values_processed = var_values[:,var_values[criterion_idx,:]>=var_threshold_value]
                var_values_processing = var_values_processing[:,var_values_processing[criterion_idx,:]>=var_threshold_value]
            elif criterion[0] == '<':
                var_values_processed = var_values[:,var_values[criterion_idx,:]<=var_threshold_value]
                var_values_processing = var_values_processing[:,var_values_processing[criterion_idx,:]<=var_threshold_value]
            elif criterion[0] == '=':
                var_values_processed = var_values[:,var_values[criterion_idx,:]==var_threshold_value]
                var_values_processing = var_values_processing[:,var_values_processing[criterion_idx,:]==var_threshold_value]
        else:
            var_threshold_value = float(criterion[1:])
            if criterion[0] == '>':
                var_values_processed = var_values[:,var_values[criterion_idx,:]>var_threshold_value]
                var_values_processing = var_values_processing[:,var_values_processing[criterion_idx,:]>var_threshold_value]
            elif criterion[0] == '<':
                var_values_processed = var_values[:,var_values[criterion_idx,:]<var_threshold_value]
                var_values_processing = var_values_processing[:,var_values_processing[criterion_idx,:]<var_threshold_value]
        
        proportion_of_var_values_that_meet_criterion = len(var_values_processed[0,:]) / len(lines)
        print(f'Percentage of var values that meet criteria #{criterion_idx+1} = {100*proportion_of_var_values_that_meet_criterion:.3f}%')

    proportion_of_var_values_that_meet_criteria = len(var_values_processing[0,:]) / len(lines)
    print(f'Percentage of var values that meet all criteria = {100*proportion_of_var_values_that_meet_criteria:.3f}%')