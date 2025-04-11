def wind_rose_from_Vortex_timeseries(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, number_of_wind_direction_sectors, wind_direction_frequency_bins_for_rose, wind_speed_bins_for_rose, max_number_of_subplot_columns):
    #For plotting windroses on subplot axes as done here, you need to edit windrose.py: see  https://stackoverflow.com/questions/42733194/subplot-of-windrose-in-matplotlib
    #usage example:
    # import EA_utils.wind_data_visualisation as wind_data_vis
    # import numpy as np
    #####################################################
    # timeseries_names = ["Bowdun-175m", "MachairWind-175m", "SpioradNaMara-175m"]
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\AR9\inputs\wind_data\timeseries\Vortex"
    # input_filenames = ["Bowdun_vortex_serie_344255_20y_175m_UTCplus00_0_ERA5.txt", "MachairWind_vortex_serie_757643_20y_175m_UTCplus01_0_ERA5.txt","SpioradNaMara_vortex_serie_757641_20y_175m_UTCplus01_0_ERA5.txt"]
    # output_file_directory_path = fr"{input_file_directory_path}\wind_roses"
    # output_filename = "AR9_TS_windroses_175m.png"
    # number_of_wind_direction_sectors = 12
    # wind_direction_frequency_bins_for_rose = np.arange(5, 25, step=5) #in % 
    # wind_speed_bins_for_rose = np.arange(0, 25, step=5) #in m/s
    # max_number_of_subplot_columns = 1
    #####################################################
    # wind_data_vis.wind_rose_from_Vortex_timeseries(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, number_of_wind_direction_sectors, wind_direction_frequency_bins_for_rose, wind_speed_bins_for_rose, max_number_of_subplot_columns)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import windrose
    
    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    n_subplot_cols = int(np.min([max_number_of_subplot_columns, len(timeseries_names)]))
    n_subplot_rows = int(np.ceil(len(timeseries_names)/n_subplot_cols))
    fig = plt.figure(figsize=(n_subplot_cols*10, n_subplot_rows*10))

    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    LARGE_SIZE = 20
    plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE) # fontsize of the legend
    plt.rc('figure', titlesize=LARGE_SIZE) # fontsize of the figure title

    ws_data_column_number = 2
    wd_data_column_number = 3
    for plot_idx, timeseries_name in enumerate(timeseries_names):
        input_filename = input_filenames[plot_idx]
        input_file_path = fr"{input_file_directory_path}\{input_filename}"
        df = pd.read_csv(input_file_path, delimiter=r"\s+", header=None, skiprows=1)
        ws = df[ws_data_column_number]
        wd = df[wd_data_column_number]

        ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, plot_idx+1, projection="windrose")
        ax.bar(wd, ws, nsector=number_of_wind_direction_sectors, normed=True, bins=wind_speed_bins_for_rose, opening=0.8, edgecolor="white")
        ax.set_yticks(wind_direction_frequency_bins_for_rose)
        ax.set_yticklabels([f"{wd}%" for wd in wind_direction_frequency_bins_for_rose])
        ax.legend(loc="lower right", title="Wind speed [$m \cdot s^{-1}$]")
        plt.title(timeseries_name,fontsize=LARGE_SIZE)

    output_file_path = fr"{output_file_directory_path}\{output_filename}"
    print("Writing to: " + output_file_path)
    plt.savefig(output_file_path)
    
def wind_rose_from_TAB_file(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, number_of_wind_direction_sectors, wind_direction_frequency_bins_for_rose, wind_speed_bins_for_rose, max_number_of_subplot_columns):
    #For plotting windroses on subplot axes as done here, you need to edit windrose.py: see  https://stackoverflow.com/questions/42733194/subplot-of-windrose-in-matplotlib
    #usage example:
    # import EA_utils.wind_data_visualisation as wind_data_vis
    # import numpy as np
    #####################################################
    # timeseries_names = ["TonnNua-155m", "TonnNua-175m"]
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Celtic_Sea\20250120_Tonn_Nua_rev5_layouts_and_EYAs\inputs\wind_data\timeseries"
    # input_filenames = ["155m_12sector_K2.tab", "175m_12sector_K2.tab"]
    # output_file_directory_path = fr"{input_file_directory_path}\wind_roses"
    # output_filename = "Tonn_Nua_Vortex_vs_K2_TS_windroses.png"
    # wind_direction_frequency_bins_for_rose = np.arange(5, 25, step=5) #in % 
    # wind_speed_bins_for_rose = np.arange(0, 25, step=5) #in m/s
    # max_number_of_subplot_columns = 1
    #####################################################
    # wind_data_vis.wind_rose_from_TAB_file(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, number_of_wind_direction_sectors, wind_direction_frequency_bins_for_rose, wind_speed_bins_for_rose, max_number_of_subplot_columns)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import windrose

    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    n_subplot_cols = int(np.min([max_number_of_subplot_columns, len(timeseries_names)]))
    n_subplot_rows = int(np.ceil(len(timeseries_names)/n_subplot_cols))
    fig = plt.figure(figsize=(n_subplot_cols*10, n_subplot_rows*10))

    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    LARGE_SIZE = 20
    plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE) # fontsize of the legend
    plt.rc('figure', titlesize=LARGE_SIZE) # fontsize of the figure title

    for plot_idx, timeseries_name in enumerate(timeseries_names):
        input_filename = input_filenames[plot_idx]
        input_file_path = fr"{input_file_directory_path}\{input_filename}"
        df = pd.read_csv(input_file_path, delimiter=r"\s+", header=0, index_col=0, skiprows=3)
        ws_bin_edges = np.append(0, df.index)
        ws_bin_centres = (ws_bin_edges[1:]+ws_bin_edges[0:-1])/2
        wd_sector_frequencies = df.columns
        number_of_wind_direction_sectors = len(wd_sector_frequencies)
        wd_sector_centres = pd.array(np.arange(0,360,360/number_of_wind_direction_sectors))
        df_wswdfreq = (df * wd_sector_frequencies.astype(float).to_list()).astype(int)
        wd = []
        ws = []
        for wd_sector_idx, wd_sector in enumerate(wd_sector_centres):
            for ws_bin_idx, ws_bin_centre in enumerate(ws_bin_centres):
                freq = df_wswdfreq.iat[ws_bin_idx, wd_sector_idx]
                for _ in range(freq):
                    wd = np.append(wd, wd_sector)
                    ws = np.append(ws, ws_bin_centre)

        ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, plot_idx+1, projection="windrose")
        ax.bar(wd, ws, nsector=number_of_wind_direction_sectors, normed=True, bins=wind_speed_bins_for_rose, opening=0.8, edgecolor="white")
        ax.set_yticks(wind_direction_frequency_bins_for_rose)
        ax.set_yticklabels([f"{wd}%" for wd in wind_direction_frequency_bins_for_rose])
        ax.legend(loc="lower right", title="Wind speed [$m \cdot s^{-1}$]")
        plt.title(timeseries_name,fontsize=LARGE_SIZE)

    output_file_path = fr"{output_file_directory_path}\{output_filename}"
    print("Writing to: " + output_file_path)
    plt.savefig(output_file_path)
    
def wind_speed_frequency_distribution_from_Vortex_timeseries(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, max_number_of_subplot_columns, plot_ws_lims, plot_freq_lims):
    #usage example:
    # import EA_utils.wind_data_visualisation as wind_data_vis
    # import numpy as np
    # ####################################################
    # timeseries_names = ["TonnNua-155m_Vortex", "TonnNua-175m_Vortex"]
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Celtic_Sea\20250120_Tonn_Nua_rev5_layouts_and_EYAs\inputs\wind_data\timeseries"
    # input_filenames = ["vortex_serie_800841_20y_155m_UTCplus00_0_ERA5.txt", "vortex_serie_800841_20y_175m_UTCplus00_0_ERA5.txt"]
    # output_file_directory_path = fr"{input_file_directory_path}\wind_freq_dist"
    # output_filename = "Tonn_Nua_Vortex_TS_windfreqdist.png"
    # max_number_of_subplot_columns = 2
    # plot_ws_lims = [-2,38] #leave empty for auto
    # plot_freq_lims = [0,9] #leave empty for auto
    # ####################################################
    # wind_data_vis.wind_speed_frequency_distribution_from_Vortex_timeseries(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, max_number_of_subplot_columns, plot_ws_lims, plot_freq_lims)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import EA_utils.statistics as stats

    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    n_subplot_cols = int(np.min([max_number_of_subplot_columns, len(timeseries_names)]))
    n_subplot_rows = int(np.ceil(len(timeseries_names)/n_subplot_cols))
    fig = plt.figure(figsize=(n_subplot_cols*10, n_subplot_rows*10))

    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    LARGE_SIZE = 20
    plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE) # fontsize of the legend
    plt.rc('figure', titlesize=LARGE_SIZE) # fontsize of the figure title

    ws_data_column_number = 2
    wd_data_column_number = 3
    for plot_idx, timeseries_name in enumerate(timeseries_names):
        input_filename = input_filenames[plot_idx]
        input_file_path = fr"{input_file_directory_path}\{input_filename}"
        df = pd.read_csv(input_file_path, delimiter=r"\s+", header=None, skiprows=1)
        ws = df[ws_data_column_number]
        [weib_ws_bins, weib_dist, weib_A, weib_k] = stats.weibull(ws)
        wind_speed_bins_for_plot = np.arange(np.min(ws), np.max(ws), step=1) #in m/s #wind speed bin interval must equal 1, as this is assumed in derivation of Weibull distribution
        
        bin_match_ids=np.digitize(ws, wind_speed_bins_for_plot)-1
        ws_freqs = np.bincount(bin_match_ids)/len(ws)
            
        ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, plot_idx+1)
        ax.bar(wind_speed_bins_for_plot, 100*ws_freqs)
        plt.plot(weib_ws_bins, 100*weib_dist, color='r')
        if plot_ws_lims:
            plt.xlim(plot_ws_lims)
        if plot_freq_lims:
            plt.ylim(plot_freq_lims)
        plt.title(timeseries_name,fontsize=LARGE_SIZE)
        plt.xlabel('Wind Speed (m/s)')
        plt.ylabel('Frequency (%)')
        plt.text(0.99, 0.99, f'Mean: {np.nanmean(ws):.3f}, Weibull A: {weib_A:.3f}, Weibull k: {weib_k:.3f}', horizontalalignment='right', verticalalignment='top', transform = ax.transAxes)

    output_file_path = fr"{output_file_directory_path}\{output_filename}"
    print("Writing to: " + output_file_path)
    plt.savefig(output_file_path)
    
def wind_speed_frequency_distribution_from_TAB_file(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, wind_speed_bins_for_plots, max_number_of_subplot_columns, plot_ws_lims, plot_freq_lims):
    #usage example:
    # import EA_utils.wind_data_visualisation as wind_data_vis
    # import numpy as np
    # ####################################################
    # timeseries_names = ["TonnNua-155m_Vortex", "TonnNua-175m_Vortex", "TonnNua-155m_K2", "TonnNua-175m_K2"]
    # input_file_directory_path = r"C:\Support\WindFarmerAnalyst\Projects\Celtic_Sea\20250120_Tonn_Nua_rev5_layouts_and_EYAs\inputs\wind_data\timeseries"
    # input_filenames = ["vortex_serie_800841_20y_155m_UTCplus00_0_ERA5-Exported.tab", "vortex_serie_800841_20y_175m_UTCplus00_0_ERA5-Exported.tab", "155m_12sector_K2.tab", "175m_12sector_K2.tab"]
    # output_file_directory_path = fr"{input_file_directory_path}\wind_freq_dist"
    # output_filename = "Tonn_Nua_Vortex_vs_K2_TS_windfreqdist.png"
    # wind_speed_bins_for_plots = [] #List of values [e.g. np.arange(0, 35, step=1).tolist()], or empty list to use bins from TAB file. For a direct comparison between data and Weibull distribution, each wind speed bin should be 1 m/s apart.
    # max_number_of_subplot_columns = 2
    # plot_ws_lims = [-2,38] #leave empty for auto
    # plot_freq_lims = [0,9] #leave empty for auto
    # ####################################################
    # wind_data_vis.wind_speed_frequency_distribution_from_TAB_file(timeseries_names, input_file_directory_path, input_filenames, output_file_directory_path, output_filename, wind_speed_bins_for_plots, max_number_of_subplot_columns, plot_ws_lims, plot_freq_lims)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import EA_utils.statistics as stats

    if not os.path.exists(output_file_directory_path):
        os.makedirs(output_file_directory_path)

    n_subplot_cols = int(np.min([max_number_of_subplot_columns, len(timeseries_names)]))
    n_subplot_rows = int(np.ceil(len(timeseries_names)/n_subplot_cols))
    fig = plt.figure(figsize=(n_subplot_cols*10, n_subplot_rows*10))

    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    LARGE_SIZE = 20
    plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE) # fontsize of the legend
    plt.rc('figure', titlesize=LARGE_SIZE) # fontsize of the figure title

    for plot_idx, timeseries_name in enumerate(timeseries_names):
        input_filename = input_filenames[plot_idx]
        input_file_path = fr"{input_file_directory_path}\{input_filename}"
        df = pd.read_csv(input_file_path, delimiter=r"\s+", header=0, index_col=0, skiprows=3)
        ws_bin_edges = np.append(0, df.index)
        ws_bin_centres = (ws_bin_edges[1:]+ws_bin_edges[0:-1])/2
        wd_sector_frequencies = df.columns
        number_of_wind_direction_sectors = len(wd_sector_frequencies)
        wd_sector_centres = pd.array(np.arange(0,360,360/number_of_wind_direction_sectors))
        df_wswdfreq = (df * wd_sector_frequencies.astype(float).to_list()).astype(int)
        ws = []
        for wd_sector_idx, wd_sector in enumerate(wd_sector_centres):
            for ws_bin_idx, ws_bin_centre in enumerate(ws_bin_centres):
                freq = df_wswdfreq.iat[ws_bin_idx, wd_sector_idx]
                for _ in range(freq):
                    ws = np.append(ws, ws_bin_centre)
        [weib_ws_bins, weib_dist, weib_A, weib_k] = stats.weibull(ws)
        if not wind_speed_bins_for_plots:
            wind_speed_bins_for_plot = ws_bin_centres
        else:
            wind_speed_bins_for_plot = wind_speed_bins_for_plots
        
        df_wsfreq = df_wswdfreq.sum(axis=1)
        bin_match_ids=np.digitize(ws_bin_centres, wind_speed_bins_for_plot)-1
        ws_freqs = np.empty(len(wind_speed_bins_for_plot))
        for ws_bin_for_plot_idx, ws_bin_for_plot in enumerate(wind_speed_bins_for_plot):
            ws_bins_to_include = [bin_ids for bin_ids, bin_match_idx in enumerate(bin_match_ids) if bin_match_idx==ws_bin_for_plot_idx]
            ws_freqs[ws_bin_for_plot_idx] = df_wsfreq.iloc[ws_bins_to_include].sum()/df_wsfreq.sum()
            
        ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, plot_idx+1)
        ax.bar(wind_speed_bins_for_plot, 100*ws_freqs)
        plt.plot(weib_ws_bins, 100*weib_dist, color='r')
        if plot_ws_lims:
            plt.xlim(plot_ws_lims)
        if plot_freq_lims:
            plt.ylim(plot_freq_lims)
        plt.title(timeseries_name,fontsize=LARGE_SIZE)
        plt.xlabel('Wind Speed (m/s)')
        plt.ylabel('Frequency (%)')
        plt.text(0.99, 0.99, f'Mean: {np.nanmean(ws):.3f}, Weibull A: {weib_A:.3f}, Weibull k: {weib_k:.3f}', horizontalalignment='right', verticalalignment='top', transform = ax.transAxes)

    output_file_path = fr"{output_file_directory_path}\{output_filename}"
    print("Writing to: " + output_file_path)
    plt.savefig(output_file_path)