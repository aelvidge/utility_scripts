def concatenate_NetCDF_files(input_directory, string_globs_in_filenames_of_files_to_be_combined, output_file_name):
    #usage example:
    #####################################################
    # input_directory = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\CelticSeas\ConstrainingWindUncertainty\Data\Era5_data\1hrly"
    # string_globs_in_filenames_of_files_to_be_combined = "*V100m.nc"
    # output_file_name = "2004-23_singlelevs_1hrly_V100m_combined.nc"
    #####################################################
    import xarray as xr

    ds = xr.open_mfdataset(input_directory + "\\" + string_globs_in_filenames_of_files_to_be_combined, combine = "by_coords")
    ds.to_netcdf(input_directory + "\\" + output_file_name) # Export netcdf file
    print("New NetCDF file saved here: " + input_directory + "\\" + output_file_name)