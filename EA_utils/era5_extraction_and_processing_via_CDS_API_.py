def era5_extraction_API_singlelevels_yearatatime(site_name, site_latitude, site_longitude, snap_extraction_location_to_nearest_Era5_native_gridpoint, start_date, end_date, variables, hour_interval, allow_inclusion_of_unverified_initial_release_data, output_file_directory_path, Era5_native_grid_latlon_txtfile_path, do_map_plot, do_data_extraction):
    #usage example:
    # import era5_extraction_and_processing_via_CDS_API.era5_extraction_API_singlelevels_yearatatime as era5_extract
    #####################################################################################################
    # site_name = 'Creyab_Windfarm'
    # site_latitude = 54.114
    # site_longitude = -3.917
    # snap_extraction_location_to_nearest_Era5_native_gridpoint = false
    # start_date = "2000-01-01"
    # end_date = "2023-12-31"
    # variables = ['100m_u_component_of_wind', '100m_v_component_of_wind']
    # hour_interval = 1
    # allow_inclusion_of_unverified_initial_release_data = False
    # output_file_directory_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\R&D\tool_development\Climate Intelligence Tool\trial data - Creyab Windfarm\ERA-5 data"
    # Era5_native_grid_latlon_txtfile_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Resources\Reanalyses\Era5\era5_winds10m_global_native_grid_latlon.txt"
    # do_map_plot = True
    # do_data_extraction = True
    #####################################################################################################
    # era5_extract(site_name, site_latitude, site_longitude, snap_extraction_location_to_nearest_Era5_native_gridpoint, start_date, end_date, variables, hour_interval, allow_inclusion_of_unverified_initial_release_data, output_file_directory_path, Era5_native_grid_latlon_txtfile_path, do_map_plot, do_data_extraction)
    
    # Extracts and concatenates ERA5 single-level data using the CDS API and saves it as a netcdf file
    # For information on CDS API, see https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
    # For information on the ERA5 grid, see https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    
    # Parameters:
    
    # site_name : string
    #     e.g. "M5 Buoy"
    
    # site_latitude : float
    #     e.g. 51.69
    #     Extracted data will be interpolated to this latitude from the native ERA5 grid
    
    # site_longitude : float
    #     e.g. -6.70
    #     Extracted data will be interpolated to this longitude from the native ERA5 grid
    
    # snap_extraction_location_to_nearest_Era5_native_gridpoint : boolean
    #     True or False
    
    # start_date : string
    #     As "YYYY-MM-DD"
    #     e.g. "2003-02-01"
    
    # end_date : string
    #     As "YYYY-MM-DD"
    #     e.g. "2024-03-16"
    
    # variables : list of strings
    #     e.g. ['10m_u_component_of_wind', '10m_v_component_of_wind','100m_u_component_of_wind', '100m_v_component_of_wind']
    #     For the full list of available variables, see https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
    
    # hour_interval : integer
    #     e.g. 1, which would give hourly output
    
    # allow_inclusion_of_unverified_initial_release_data : boolean
    #     True or False
    
    # output_file_directory_path : string
    #     e.g. r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects\CelticSea\01 Analysis\20231214_ConstrainingWindUncertainty\data\Era5_data\at_buoy_location"
    
    # Era5_native_grid_latlon_txtfile_path : string
    #     e.g. r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Resources\Reanalyses\Era5\era5_winds10m_global_native_grid_latlon.txt"
    #     Gives location of era5_winds10m_global_native_grid_latlon.txt
    
    # do_map_plot : boolean
    #     True or False
    
    # do_data_extraction : boolean
    #     True or False
    
    import cdsapi # The CDS API client is a python based library, used to request data from the datasets listed in the CDS catalogue
    import numpy as np
    import xarray as xr
    import pandas as pd
    from geopy import distance
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import datetime
    
    """Load Era5 native grid coordinates"""
    df = pd.read_csv(Era5_native_grid_latlon_txtfile_path, delimiter='\t')
    native_grid_lats = df["grid_latitude"].values
    native_grid_lons = df["grid_longitude"].values

    """If site longitude is closer to 180 degrees than 0 degrees, add 360 to negative longitudes"""
    if site_longitude <= -90 or site_longitude > 90:
        if site_longitude <= -90: site_longitude = site_longitude + 360
        native_grid_lons[native_grid_lons<0] = native_grid_lons[native_grid_lons<0] + 360

    """Find, and derive distance from, nearest native Era5 grid point to site location"""
    dlat = np.abs(native_grid_lats-site_latitude)
    dlon = np.abs(native_grid_lons-site_longitude)
    nearest_native_gridpoints_ids = np.argsort(np.maximum(dlon,dlat))[0:4] #finds the four nearest grid points
    dist_to_nearest_native_gridpoints = [None]*4
    for idx, gridpoints_idx in enumerate(nearest_native_gridpoints_ids):
        #calculate distance to each of the four nearest grid points
        dist_to_nearest_native_gridpoints[idx] = distance.distance((native_grid_lats[gridpoints_idx],native_grid_lons[gridpoints_idx]),(site_latitude,site_longitude)).km
    nearest_native_gridpoint_idx = nearest_native_gridpoints_ids[np.argmin(dist_to_nearest_native_gridpoints)]
    if snap_extraction_location_to_nearest_Era5_native_gridpoint:
        site_latitude, site_longitude = native_grid_lats[nearest_native_gridpoint_idx], native_grid_lons[nearest_native_gridpoint_idx]
        dist_to_nearest_native_gridpoint = 0
    else:
        dist_to_nearest_native_gridpoint = distance.distance((native_grid_lats[nearest_native_gridpoint_idx],native_grid_lons[nearest_native_gridpoint_idx]),(site_latitude,site_longitude)).km
    print(f"Nearest native Era5 gridpoint to site is {dist_to_nearest_native_gridpoint:.2f} km away")

    if do_map_plot:
        """Generate and save figure of site location over Era5 native grid and coastline"""
        site_colour = [0,0,1]
        site_symbol = 'o'
        native_grid_colour = [0.8,0.4,0.4]
        native_grid_symbol = 'x'

        plt.figure(figsize=(7,7))
        ax = plt.axes(projection=ccrs.Stereographic(central_latitude=site_latitude, central_longitude=site_longitude))
        #ax.add_feature(cfeature.OCEAN) #commented out as computationally expensive
        ax.add_feature(cfeature.LAND,facecolor=[0.85,0.85,0.85])
        ax.coastlines(linewidth=0.75, color='k')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle=":", color=[0.7,0.7,0.7], draw_labels=True, x_inline=False, y_inline=False)
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180,1))
        gl.ylocator = mticker.FixedLocator(np.arange(-90,90,1))
        gl.xlabel_style = {'fontsize': 8, 'color': [0.5,0.5,0.5]}
        gl.ylabel_style = {'fontsize': 8, 'color': [0.5,0.5,0.5]}
        ax.set_extent([site_longitude-1.9, site_longitude+1.9, site_latitude-1.1, site_latitude+1.1], crs=ccrs.PlateCarree())
        plt.scatter(site_longitude, site_latitude, s=80, color=site_colour, marker=site_symbol, transform=ccrs.PlateCarree())
        plt.scatter(native_grid_lons, native_grid_lats, s=20, color=native_grid_colour, marker=native_grid_symbol, transform=ccrs.PlateCarree())
        plt.scatter(native_grid_lons[nearest_native_gridpoint_idx], native_grid_lats[nearest_native_gridpoint_idx], s=60, edgecolor=native_grid_colour, facecolor='none', marker='o', transform=ccrs.PlateCarree())
        plt.text(site_longitude, site_latitude+0.07, site_name, color='b', horizontalalignment='center', transform=ccrs.PlateCarree(), fontsize=10)
        plt.title(f"Nearest native Era5 gridpoint is {dist_to_nearest_native_gridpoint:.2f} km away", fontsize=8, color=native_grid_colour)
        plt.savefig(f"{output_file_directory_path}\\map_{site_name}_and_nearby_Era5_gridpoints.png", dpi=150)
        print(f"Map of site and Era5 gridpoints saved here: {output_file_directory_path}\\map_{site_name}_and_nearby_Era5_gridpoints.png")

    if do_data_extraction:
        """Ensure requested period falls within period of Era5 data availability"""
        Era5_start_date = "1940-01-01"
        Era5_end_date = (datetime.date.today() - datetime.timedelta(days = 6)).strftime('%Y-%m-%d')
        if datetime.datetime.strptime(start_date, '%Y-%m-%d').date() < datetime.datetime.strptime(Era5_start_date, '%Y-%m-%d').date():
            start_date = Era5_start_date
            print(f"Start date requested is before the start of the Era-5 record, so has been updated accordingly, to {start_date}")
        if datetime.datetime.strptime(end_date, '%Y-%m-%d').date() > datetime.datetime.strptime(Era5_end_date, '%Y-%m-%d').date():
            end_date = Era5_end_date
            print(f"End date requested is after the date of the most recently available Era-5 data, so has been updated accordingly, to {end_date}")
        
        """Define CDS API data-retieve request"""
        def CDS_retrieve_era5_data(client, variables, extraction_year, extraction_months, extraction_days, hour_interval, site_latitude, site_longitude, output_file_path):
            client.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'variable': variables,
                    'year': str(extraction_year),
                    'month': [str(extraction_month) for extraction_month in extraction_months],
                    'day': [str(extraction_day) for extraction_day in extraction_days],
                    'time': [f'{hr:02}:00' for hr in np.arange(0,24,hour_interval)],
                    'area': [site_latitude, site_longitude, site_latitude, site_longitude],
                    'format': 'netcdf',
                    },
                output_file_path)
            print(f"Saved Era5 data to {output_file_path}\n")

        """Characterise extraction start and end dates"""
        start_year = int(start_date[0:4])
        end_year = int(end_date[0:4])
        start_month = int(start_date[5:7])
        end_month = int(end_date[5:7])
        start_day = int(start_date[8:10])
        end_day = int(end_date[8:10])
        extraction_years = np.arange(start_year, end_year+1, 1)

        """Place Era-5 data extraction requests via CDS API. This is done year-at-a-time (longer periods can cause timeout failure), with the exception of the final year of extractions that include the most recent month for which Era5 data is available, for which month-at-a-time extraction is necessary."""
        CDS_API_client = cdsapi.Client()
        output_file_paths = []
        for extraction_year in extraction_years:
            # define months of year to request data for
            if extraction_year == start_year: extraction_month_start = start_month
            else: extraction_month_start = 1
            if extraction_year == end_year: extraction_month_end = end_month
            else: extraction_month_end = 12
            extraction_months = np.arange(extraction_month_start, extraction_month_end+1, 1)
            
            # if extraction includes data from the most recent month for which Era5 data is available, need to extract month-at-a-time with days defined in order to avoid requesting data that is not yet available, which will cause extraction failure
            # otherwise, extract data for the full year
            if extraction_year == Era5_end_date[0:4] and end_month == Era5_end_date[0:4]:
                for extraction_month in extraction_months:
                    # define days of month to request data for
                    if extraction_year == start_year and extraction_month == start_month: extraction_day_start = start_day
                    else: extraction_day_start = 1
                    if extraction_year == end_year and extraction_month == end_month: extraction_day_end = end_day
                    else: extraction_day_end = 31
                    extraction_days = np.arange(extraction_day_start, extraction_day_end+1, 1)
                    
                    # define export location for NetCDF files
                    output_file_path = fr"{output_file_directory_path}\Era5_{site_name}_{extraction_year}_{extraction_month}.nc"
                    output_file_paths.append(output_file_path)
                    
                    # place CDS API request
                    print(f"Extracting data for location [{site_latitude}, {site_longitude}] and period [{extraction_year}-{extraction_month:02d}-{extraction_day_start:02d} to {extraction_year}-{extraction_month:02d}-{extraction_day_end:02d} inclusive]")
                    CDS_retrieve_era5_data(CDS_API_client,variables,extraction_year,[extraction_month],extraction_days,hour_interval,site_latitude,site_longitude,output_file_path)
                    
            else:
                extraction_days = np.arange(1, 32, 1)
                
                # define export location for NetCDF files
                output_file_path = fr"{output_file_directory_path}\Era5_{site_name}_{extraction_year}.nc"
                output_file_paths.append(output_file_path)
                
                # place CDS API request
                print(f"Extracting data for the period {extraction_year}-{extraction_month_start:02d} to {extraction_year}-{extraction_month_end:02d} inclusive")
                CDS_retrieve_era5_data(CDS_API_client,variables,extraction_year,extraction_months,extraction_days,hour_interval,site_latitude,site_longitude,output_file_path)

        """Load all extracted data from NetCDF files"""
        datasets = [[]]*len(output_file_paths)
        for output_file_path_idx, output_file_path in enumerate(output_file_paths):
            ds = xr.open_dataset(output_file_path)
            if "expver" in ds.coords: # The possible occurrence of an additional dimension, "expver", needs accounting for. This parameter specifies whether the data is verified or initial-release Era5 data. Initial release Era5 is only available up to 5 days prior to the most recently available data.
                print(output_file_path)
                var_names = list(ds.keys())
                for time_idx in range(len(ds.coords['time'])):
                    if np.isnan(ds[var_names[0]].values[time_idx,0,0,0]):
                        if allow_inclusion_of_unverified_initial_release_data:
                            for var_name in var_names:
                                ds[var_name][time_idx,0,0,0] = ds[var_name][time_idx,1,0,0]
                            print(f"Warning: data for {ds.coords['time'].values[time_idx]} is unverified 'initial release' data")
                        else:
                            print(f"Warning: removing unverified 'initial release' data for {ds.coords['time'].values[time_idx]}")
                datasets[output_file_path_idx] = ds.sel(expver=1, drop=True)
            else:
                datasets[output_file_path_idx] = ds

        """Concatenate all data along the "time" dimension and export as one NetCDF file for the full period"""
        ds = xr.concat(datasets, dim="time")
        ds = ds.sel(time=slice(start_date, end_date))
        combined_output_file_path = f"{output_file_directory_path}" + fr"\Era5_{site_name}_{start_date}_to_{end_date}.nc".replace("-","")
        ds.to_netcdf(combined_output_file_path) # Export netcdf file
        print(f"New combined NetCDF file saved here: {combined_output_file_path}")
        
        
def era5_extraction_API_modellevels_yearatatime(site_name, site_latitude, site_longitude, site_height, start_date, end_date, variables, hour_interval, output_file_directory_path, Era5_native_grid_latlon_txtfile_path):
    """
    Extracts ERA5 model-level data using the CDS API
    For information on CDS API, see https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
    For information on the ERA5 grid, see https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    """
    #usage example:
    ################################################################################
    # site_name = "HesseloKattegat_TS"
    # site_latitude = 56.381479
    # site_longitude = 11.475470
    # site_height = 180
    # start_date = "2013-01-01" #must have format YYYY-MM-DD
    # end_date = "2023-10-11" #must have format YYYY-MM-DD
    # variables = ['10m_u_component_of_wind', '10m_v_component_of_wind','100m_u_component_of_wind', '100m_v_component_of_wind']
    # hour_interval = 1
    # output_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Denmark\inputs\wind_data\timeseries\Era5\10-year_TS\on_model_levels"
    # Era5_native_grid_latlon_txtfile_path = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Resources\Reanalyses\Era5\era5_winds10m_global_native_grid_latlon.txt"
    ################################################################################
    import cdsapi
    import numpy as np
    import xarray as xr
    import pandas as pd
    from geopy import distance
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import datetime
    import eccodes
    import cfgrib
    
    df = pd.read_csv(Era5_native_grid_latlon_txtfile_path, delimiter='\t')
    grid_lat = df["grid_latitude"].values
    grid_lon = df["grid_longitude"].values
    abslat = np.abs(grid_lat-site_latitude)
    abslon = np.abs(grid_lon-site_longitude)
    nearest_gridpoint_i = np.argmin(np.maximum(abslon,abslat))
    dist = distance.distance((grid_lat[nearest_gridpoint_i],grid_lon[nearest_gridpoint_i]),(site_latitude,site_longitude)).km
    print(f"Nearest Era5 gridpoint to site is {dist:.2f} km away")
    #dist = [[]]*len(grid_lat)
    #for i in range(len(grid_lat)):
    #    dist[i] = distance.distance((grid_lat[nearest_gridpoint_i],grid_lon[nearest_gridpoint_i]),(site_latitude,site_longitude))
    #    print(i/len(grid_lat))

    ax = plt.axes(projection=ccrs.Stereographic(central_latitude=site_latitude, central_longitude=site_longitude))
    #ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND,facecolor=[0.8,0.8,0.8])
    ax.coastlines(linewidth=0.75, color='k')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle="--", color=[0.5,0.5,0.5], draw_labels=True, x_inline=False, y_inline=False)
    gl.xlocator = mticker.FixedLocator(np.arange(-180,180,1))
    gl.ylocator = mticker.FixedLocator(np.arange(-90,90,1))
    gl.xlabel_style = {'fontsize': 8, 'color': [0.5,0.5,0.5]}
    gl.ylabel_style = {'fontsize': 8, 'color': [0.5,0.5,0.5]}
    ax.set_extent([site_longitude-1.9, site_longitude+1.9, site_latitude-1.1, site_latitude+1.1], crs=ccrs.PlateCarree())
    plt.scatter(site_longitude, site_latitude, s=80, color='b', marker='o', transform=ccrs.PlateCarree())
    plt.scatter(grid_lon, grid_lat, s=20, color=[1,0.5,0.5], marker='x', transform=ccrs.PlateCarree())
    plt.scatter(grid_lon[nearest_gridpoint_i], grid_lat[nearest_gridpoint_i], s=60, color=[1,0.5,0.5], marker='x', transform=ccrs.PlateCarree())
    plt.text(site_longitude, site_latitude+0.07, site_name, color='b', horizontalalignment='center', transform=ccrs.PlateCarree())
    plt.title(f"Nearest Era5 gridpoint is {dist:.2f} km away", fontsize=8, color=[1,0.5,0.5])
    plt.savefig(f"{output_file_directory_path}\\map_site_and_nearby_Era5_gridpoints.png", dpi=150)
    print(f"Map of site and Era5 gridpoints saved here: {output_file_directory_path}\\map_site_and_nearby_Era5_gridpoints.png")

    model_levels = np.arange(127,137,1)
    heights_of_model_levels = np.array([334.24, 287.52, 244.69, 205.44, 169.51, 136.62, 106.54, 79.04, 53.92, 30.96, 10])
    i = np.absolute(np.asarray(heights_of_model_levels)-site_height).argmin()
    c = int(site_height<heights_of_model_levels[i])
    r = np.arange(i-2+c,i+2+c,1)
    model_levels_required = "/".join(str(ele) for ele in model_levels[r])
    print(f"Heights of model levels to be requested: {heights_of_model_levels[r]}")

    date_of_earliest_Era5_data = datetime.datetime.strptime("1940-01-01", '%Y-%m-%d').date()
    if datetime.datetime.strptime(end_date, '%Y-%m-%d').date() < date_of_earliest_Era5_data:
        start_date = date_of_earliest_Era5_data

    date_of_latest_Era5_data = datetime.date.today() - datetime.timedelta(days = 6)
    if datetime.datetime.strptime(end_date, '%Y-%m-%d').date() > date_of_latest_Era5_data:
        end_date = date_of_latest_Era5_data

    def CDS_retrieve_era5_data(variables,extraction_year,extraction_months,extraction_days,hour_interval,site_latitude,site_longitude,model_levels_required,output_file_path):
        print(f'{extraction_year}-{extraction_months[0]:0>2}-{extraction_days[0]:0>2}/to/{extraction_year}-{extraction_months[-1]:0>2}-{extraction_days[-1]:0>2}')
        print(f'{extraction_year}-{extraction_months[0]:0>2}-{extraction_days[0]:0>2}/{extraction_year}-{extraction_months[0]:0>2}-{extraction_days[0]+1:0>2}')
        c.retrieve('reanalysis-era5-complete', {# Requests follow MARS syntax
                                                # Keywords 'expver' and 'class' can be dropped. They are obsolete
                                                # since their values are imposed by 'reanalysis-era5-complete'
            'levelist': model_levels_required,  # 1 is top level, 137 the lowest model level in ERA5. Use '/' to separate values. See https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
            'levtype' : 'ml',
            'param'   : '131/132/156',          # Full information at https://apps.ecmwf.int/codes/grib/param-db/
                                                # The native representation for temperature is spherical harmonics
            'stream'  : 'oper',                 # Denotes ERA5. Ensemble members are selected by 'enda'
            'date'    : f'{extraction_year}-{extraction_months[0]:0>2}-{extraction_days[0]:0>2}/to/{extraction_year}-{extraction_months[-1]:0>2}-{extraction_days[-1]:0>2}',
            'time'    : f'00/to/23/by/{hour_interval}',
            'area'    : [site_latitude, site_longitude, site_latitude, site_longitude],
            'type'    : 'an',
            'format'  : 'nc',
            }, output_file_path)           # Output file.
        print(f"Saved Era5 data to {output_file_path}\n")
        
    year_start = int(start_date[0:4])
    year_end = int(end_date[0:4])
    month_start = int(start_date[5:7])
    month_end = int(end_date[5:7])
    day_start = int(start_date[8:10])
    day_end = int(end_date[8:10])

    extraction_years = np.arange(year_start, year_end+1, 1)

    output_file_paths = []
    c = cdsapi.Client()
    for extraction_year in extraction_years:
        if extraction_year == year_start: extraction_month_start = month_start
        else: extraction_month_start = 1
        if extraction_year == year_end: extraction_month_end = month_end
        else: extraction_month_end = 12
        extraction_months = np.arange(extraction_month_start, extraction_month_end+1, 1)
        
        if extraction_year == date_of_latest_Era5_data.year and month_end == date_of_latest_Era5_data.month:
            for extraction_month in extraction_months:
                if extraction_year == year_start and extraction_month == month_start: extraction_day_start = day_start
                else: extraction_day_start = 1
                if extraction_year == year_end and extraction_month == month_end: extraction_day_end = day_end
                else: extraction_day_end = 31
                extraction_days = np.arange(extraction_day_start, extraction_day_end+1, 1)
                output_file_path = fr"{output_file_directory_path}\Era5_{site_name}_{extraction_year}_{extraction_month}.grib"
                output_file_paths.append(output_file_path)
                print(f"Extracting data for the period {extraction_year}-{extraction_month:02d}-{extraction_day_start:02d} to {extraction_year}-{extraction_month:02d}-{extraction_day_end:02d} inclusive")
                CDS_retrieve_era5_data(variables,extraction_year,[extraction_month],extraction_days,hour_interval,site_latitude,site_longitude,model_levels_required,output_file_path)
                
        else:
            extraction_days = np.arange(1, 31+1, 1)
            output_file_path = fr"{output_file_directory_path}\Era5_{site_name}_{extraction_year}.grib"
            output_file_paths.append(output_file_path)
            print(f"Extracting data for the period {extraction_year}-{extraction_month_start:02d} to {extraction_year}-{extraction_month_end:02d} inclusive")
            CDS_retrieve_era5_data(variables,extraction_year,extraction_months,extraction_days,hour_interval,site_latitude,site_longitude,model_levels_required,output_file_path)

    datasets = [[]]*len(output_file_paths)
    for output_file_path_i,output_file_path in enumerate(output_file_paths):
        ds = xr.open_dataset(output_file_path, engine="cfgrib")
        if "expver" in ds.coords:
            for i in range(len(ds.coords["time"])):
                if np.isnan(ds.u10.values[i,0,0,0]):
                    ds.u10[i,0,0,0] = ds.u10[i,1,0,0]
                    ds.v10[i,0,0,0] = ds.v10[i,1,0,0]
                    ds.u100[i,0,0,0] = ds.u100[i,1,0,0]
                    ds.v100[i,0,0,0] = ds.v100[i,1,0,0]
            datasets[output_file_path_i] = ds.sel(expver=1, drop=True)
        else:
            datasets[output_file_path_i] = ds

    ds = xr.concat(datasets, dim="time")
    ds = ds.sel(time=slice(start_date, end_date))
    combined_output_file_path = f"{output_file_directory_path}" + fr"\Era5_{site_name}_{start_date}_to_{end_date}_.nc".replace("-","")
    ds.to_netcdf(combined_output_file_path) # Export netcdf file
    print(f"New combined NetCDF file saved here: {combined_output_file_path}")
    
    
def era5_get_native_latlon(output_file_path_grib, output_file_path_txt_latlon):
    #usage example:
    ##############################################################################################
    # output_file_path_grib = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Resources\Reanalyses\Era5\era5_winds10m_global_native_grid.grib"
    # output_file_path_txt_latlon = r"C:\Users\AE78333\OneDrive - SSE PLC\Documents\Resources\Reanalyses\Era5\era5_winds10m_global_native_grid_latlon.txt"
    ##############################################################################################
    import cdsapi
    
    c = cdsapi.Client()
    c.retrieve('reanalysis-era5-complete', {# Requests follow MARS syntax
                                            # Keywords 'expver' and 'class' can be dropped. They are obsolete
                                            # since their values are imposed by 'reanalysis-era5-complete'
    'date'    : '2023-01-01',            # The hyphens can be omitted
    'levtype' : 'sfc',
    'param'   : '165/166',               # Full information at https://apps.ecmwf.int/codes/grib/param-db/
                                            # The native representation for temperature is spherical harmonics
    'stream'  : 'oper',                  # Denotes ERA5. Ensemble members are selected by 'enda'
    'time'    : '00',                    # You can drop :00:00 and use MARS short-hand notation, instead of '00/06/12/18'
    'type'    : 'an',
    'format'  : 'grib',
    }, output_file_path_grib)               # Output file.

    import xarray as xr
    import cfgrib
    ds = xr.open_dataset(output_file_path_grib, engine="cfgrib")
    print(ds)
    grid_latitude = ds.coords['latitude'].values
    grid_longitude = ds.coords['longitude'].values
    grid_longitude[grid_longitude>180] = grid_longitude[grid_longitude>180] - 360

    import pandas as pd
    df = pd.DataFrame({'grid_latitude':grid_latitude, 'grid_longitude':grid_longitude})
    df.to_csv(output_file_path_txt_latlon, index=False, sep="\t", float_format='%1.14f')