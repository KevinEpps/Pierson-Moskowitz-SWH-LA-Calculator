# CDIP Buoy Processor

'''
This program was developed by Kevin Epps from the Engineering section 
at the Amphibious Vehicle Test Branch (AVTB), based on a MATLAB 
program developed by Mike Slivka.  Its purpose is to download 
wave buoy data from the Coastal Data Information Program (CDIP) and 
display it in a way that is useful for test planning and analysis.
Specifically, it normlaizes the significant wave height of the 
localized sea state based to a 3' Pierson Moskowitz equivalent.
It also provides a forecast of buoy height, as well as a prediction
of the normalized significant wave height using an XGBRegressor trained
model.

CDIP is run by the Scripps Institution of Oceanography
(SIO) at the University of California San Diego (UCSD).

https://cdip.ucsd.edu/
'''

import netCDF4
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz
import cftime
import math
from xgboost import XGBRegressor
import joblib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# List of buoys
buoys = ['045', '264', '271']

# Create a PdfPages object
pdf_pages = PdfPages('wave_heights.pdf')

for stn in buoys:

    pi_value = math.pi
    gravity = 32.2

    # Load the prediction model from disk
    filename = 'finalized_xgb_model.sav'
    loaded_model = joblib.load(filename)

    station_names = {
        '045': 'Oceanside Harbor',
        '179': 'Astoria Canyon',
        '162': 'Clatsop Spit',
        '264': 'Red Beach',
        '271': 'Green Beach'
        # Add more stations as needed
    }

    #Download and extract necessary variables from the netCDF file using the specificed timeframe
    def download_wave_data(stn, start_time, end_time):
        
        data_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/{stn}p1_rt.nc'
        nc = netCDF4.Dataset(data_url)
        time_var = nc.variables['waveTime']
        start_index = netCDF4.date2index(start_time, time_var, select='nearest')
        end_index = netCDF4.date2index(end_time, time_var, select='nearest')
        waveHs = nc.variables['waveHs'][start_index:end_index+1]  * 3.281  # Convert to feet
        waveFrequency = nc.variables['waveFrequency'][:]
        waveEnergyDensity = nc.variables['waveEnergyDensity'][start_index:end_index+1,:].flatten() # Flattens the array, converting it from a multi-dimensional array (m^2/Hz) to a one-dimensional array (m^2)
        waveBandwidth = nc.variables['waveBandwidth'][:]
        waveDirection = nc.variables['waveMeanDirection'][start_index:end_index+1]
        time_array = [cftime.num2pydate(t, time_var.units).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('US/Pacific')) for t in time_var[start_index:end_index+1]]
        nc.close()
        waveEnergyDensity_ft = waveEnergyDensity * 10.764 #converts from m^2 to ft^2
        return waveHs, waveFrequency, waveEnergyDensity_ft, waveBandwidth, waveDirection, time_array, waveEnergyDensity

    #Download and extract necessary variables for forecast
    def download_wave_data_forecast(stn, start_time_f, end_time_f):
        
        #forecast_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_validation/BP{stn}_forecast.nc' #Old way of forecasting on CDIP
        forecast_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_validation/BP{stn}_ecmwf_fc.nc' #Using ECMWF (European Centre for Medium-Range Weather Forecasts)
        nc = netCDF4.Dataset(forecast_url)
        time_var_f = nc.variables['waveTime']
        start_index_f = netCDF4.date2index(start_time_f, time_var_f, select='nearest')
        end_index_f = netCDF4.date2index(end_time_f, time_var_f, select='nearest')
        waveHs_f = nc.variables['waveHs'][start_index_f:end_index_f+1]  * 3.281  # Convert to feet
        waveTa_f = nc.variables['waveTa'][start_index_f:end_index_f+1]
        time_array_f = time_array = [cftime.num2pydate(t, time_var_f.units).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('US/Pacific')) for t in time_var_f[start_index_f:end_index_f+1]]
        nc.close()
        return waveHs_f, waveTa_f, time_array_f

    #Set the buoy station and timeframe
    stn = stn

    #Set start and end times for realtime data
    end_time_utc = datetime.utcnow()
    start_time_utc = end_time_utc - timedelta(hours=48)

    #Set start and end times for forecast data
    #start_time_f_utc = datetime.utcnow() - timedelta(hours=1)
    start_time_f_utc = start_time_utc
    #end_time_f_utc = start_time_f_utc + timedelta(hours=48)
    end_time_f_utc = end_time_utc + timedelta(hours=48)

    utc_tz = pytz.timezone('UTC')
    pst_tz = pytz.timezone('US/Pacific')

    #Convert realtime to PST
    end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
    start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

    #Convert forecast to PST
    end_time_f = utc_tz.localize(end_time_f_utc).astimezone(pst_tz)
    start_time_f = utc_tz.localize(start_time_f_utc).astimezone(pst_tz)

    #Execute the download data function for realtime and forecast timeframes
    waveHs, waveFrequency, waveEnergyDensity_ft, waveBandwidth, waveDirection, time_array, waveEnergyDensity = download_wave_data(stn, start_time, end_time)
    waveHs_f, waveTa_f, time_array_f = download_wave_data_forecast(stn, start_time_f, end_time_f)

    # Implement forecasting using CDIP provided values for Hs and Ta
    '''
    Uses a prediction model trained on 20+ years of historic data from Oceanside buoy to give the length adjusted value based off of Hs and Ta

    Model: XGBoost
    Mean absolute error: 0.042923021030845346
    Mean squared error: 0.00467019630602591
    R^2 score: 0.9973058696336262

    '''
    # Create forecast table list
    forecast = []

    for i in range(len(time_array_f)):
        # Take the values for the current time step
        forecast_waveHs = waveHs_f[i]
        forecast_waveTa = waveTa_f[i]
        
        # Create a DataFrame for the current time step
        df_forecast = pd.DataFrame({'Buoy Wave Height (ft)': [forecast_waveHs], 'Wave Period': [forecast_waveTa]})
        
        # Make length adjusted wave height predictions on the forecast data
        predictions = loaded_model.predict(df_forecast)
        
        # Append the current time, prediction and forecast_waveHs to the forecast list
        forecast.append((time_array_f[i], forecast_waveHs, predictions[0]))

    # Convert the forecast list to a DataFrame for easy handling
    df_forecast = pd.DataFrame(forecast, columns=['Time', 'Forecast Buoy Wave Height (ft)', 'Forecast LA SWH (ft)'])

    #Finds the swell separation frequency based on the maximum energy density with the frequency range of 0.075 to 0.10125 Hz which represent swell energy
    #Returns a dataframe for both swell and sea energy densities
    def apply_correction_factor(df, freq_range_min=0.075, freq_range_max=0.10125):
        freq_range = df[(df['Frequency'] >= freq_range_min) & (df['Frequency'] <= freq_range_max)]
        min_density_idx = freq_range['Energy Density (ft^2/Hz)'].idxmin()
        swell_separation_frequency = df.loc[min_density_idx, 'Frequency']
        swell_separation_density = df.loc[min_density_idx, 'Energy Density (ft^2/Hz)']

        #Uncomment the lines below for troubleshooting of this function
        #print(f'Minimum Swell sep density within the range of {freq_range_min} to {freq_range_max} Hz is ', swell_separation_density)
        #print('The corresponding swell separation frequency is', swell_separation_frequency)

        sea_energy_density = df['Energy Density (ft^2/Hz)'].copy()
        min_frequency = df['Frequency'].min()
        for i, freq in enumerate(df['Frequency']):
            if freq < swell_separation_frequency:
                factor = ((freq - min_frequency) / (swell_separation_frequency - min_frequency)) ** 8
                sea_energy_density[i] *= factor
        return sea_energy_density, swell_separation_density, swell_separation_frequency

    #Creates an empty list to be filled with results from the calculations below
    results = []

    #Steps through the time array, performs calculations to determine the normalization factor that transforms the localized sea state
    #to one that is normalized to a 3' Pierson Moskowitz, and adds the results to a list
    for i in range(len(time_array)):
        # Take the values for the current time step
        current_waveHs = waveHs[i]
        current_waveEnergyDensity_ft = waveEnergyDensity_ft[i * len(waveFrequency):(i + 1) * len(waveFrequency)]
        #current_waveBandwidth = waveBandwidth[i]
        current_waveDirection = waveDirection[i]
        current_waveEnergyDensity = waveEnergyDensity[i]

        df = pd.DataFrame({'Frequency': waveFrequency, 'Bandwidth': waveBandwidth, 'Energy Density (ft^2/Hz)': current_waveEnergyDensity_ft})

        #print('Raw Wave Energy Density', waveEnergyDensity) #Use as a check with the Datawell Spreadsheet

        # Apply the correction factor to get the 'Sea Energy Density (ft^2/Hz)'
        df['Sea Energy Density (ft^2/Hz)'], swell_separation_density, swell_separation_frequency = apply_correction_factor(df)

        # Calculate swell energy density
        df['Swell Energy Density (ft^2/Hz)'] = df['Energy Density (ft^2/Hz)'] - df['Sea Energy Density (ft^2/Hz)']

        # Calculate the 'Sea M0 Ordinate (ft^2)' column
        df['Sea M0 Ordinate (ft^2)'] = df['Sea Energy Density (ft^2/Hz)'] * df['Bandwidth']

        # Calculate the 'Sea M1 Ordinate (ft^2)' column
        df['Sea M1 Ordinate (ft^2-Hz)'] = df['Sea M0 Ordinate (ft^2)'] * df['Frequency']

        # Return the most recent direction for the 'Mean Direction (deg)' column
        df['Mean Direction (deg)'] = waveDirection[0]

        # Finds the total area of the sea energy density and converts it to significant height
        sea_area = df['Sea M0 Ordinate (ft^2)'].sum()
        sea_std_dev = sea_area ** 0.5
        sea_sig_height = 4 * sea_std_dev

        # Find the frequency corresponding to the maximum sea energy density
        max_sea_energy_density_freq = df.loc[df['Sea Energy Density (ft^2/Hz)'].idxmax(), 'Frequency']

        # Compute the sea ordinate and corresponding period and length
        sea_M1 = df['Sea M1 Ordinate (ft^2-Hz)'].sum()

        sea_mean_period = 1/sea_M1

        sea_avg_length = 5.12 * sea_mean_period**2 #5.12 is gravity/2PI

        sea_direction = df.loc[df['Sea Energy Density (ft^2/Hz)'].idxmax(), 'Mean Direction (deg)'] + 11.53199 #From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14
        swell_direction = df.loc[df['Swell Energy Density (ft^2/Hz)'].idxmax(), 'Mean Direction (deg)'] + 11.53199 #From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14

        # PM calculations and equivalent energy density
        # Uses the significant sea height to calculate a Pierson-Moskowitz equivalent height
        w_modal = 0.4 * (gravity/sea_sig_height)**0.5
        f_modal = w_modal/(2*pi_value)
        t_modal = 1/f_modal
        A = 0.0081*gravity**2
        B = -0.032*(gravity/sea_sig_height)**2
        rps = 0.545*((-B)**0.5)**0.5
        pm_area = -A/(4*B)
        pm_std_dev = pm_area**0.5
        pm_sig_height = 4 * pm_std_dev

        # Create a new DataFrame for PM values
        pm_df = pd.DataFrame()
        pm_df['Frequency'] = df['Frequency']
        pm_df['Bandwidth'] = df['Bandwidth']
        pm_df['Frequency (rad/sec)'] = df['Frequency'] * 2 * np.pi

        # Add the 'ft^2-sec' column to the DataFrame
        pm_df['ft^2-sec'] = A / (pm_df['Frequency (rad/sec)'] + 1e-4)**5 * np.exp(B / (pm_df['Frequency (rad/sec)'] + 1e-4)**4)

        # Replace any infinite or NaN values with 0
        pm_df['ft^2-sec'] = pm_df['ft^2-sec'].replace([np.inf, -np.inf, np.nan], 0)

        pm_df['PM Energy Density (ft^2/Hz)'] = pm_df['ft^2-sec'] * 2*pi_value

        # Calculate the 'PM M0 Ordinate (ft^2)' column
        pm_df['PM M0 Ordinate (ft^2)'] = pm_df['PM Energy Density (ft^2/Hz)'] * pm_df['Bandwidth']

        # Calculate the 'Sea M1 Ordinate (ft^2)' column
        pm_df['PM M1 Ordinate (ft^2-Hz)'] = pm_df['PM M0 Ordinate (ft^2)'] * pm_df['Frequency']

        pm_area = pm_df['PM M0 Ordinate (ft^2)'].sum()
        pm_std_dev = pm_area ** 0.5
        pm_sig_height = 4 * pm_std_dev

        # Find the frequency corresponding to the maximum sea energy density
        max_pm_energy_density_freq = pm_df.loc[pm_df['PM Energy Density (ft^2/Hz)'].idxmax(), 'Frequency']

        # Compute the modal period
        pm_modal_period = 1 / max_pm_energy_density_freq

        # Calculates the sum of the M1 ordinates
        pm_M1 = pm_df['PM M1 Ordinate (ft^2-Hz)'].sum()

        # Finds the mean period of the M1 ordinate
        pm_mean_period = 1/pm_M1

        # Calculates the average length based on the period
        pm_avg_length = 5.12 * pm_mean_period**2 #5.12 is gravity/2PI

        # Determines the normalization factor based on the ratio of the average lengths
        cm = (pm_avg_length/sea_avg_length)**0.5

        # Since the worst case scenario sea state is a 3' Pierson-Moskowitz, the normalization factor will never be greater than 1 (which represents a fully developed 3' PM)
        if cm > 1:
            cm = 1
        else:
            cm *= 1

        # Calculates additional correction factors based on sea modal period and wave direction
        # Additional height is added when the period is shorter than 7 seconds and when the sea and swell directions are within 40 degrees of each other
        # This is based off of subject matter expertise
        df['Measured M0 (ft^2)'] = df['Energy Density (ft^2/Hz)'] * df['Bandwidth']
        measured_sig_height = 4 * (df['Measured M0 (ft^2)'].sum()**0.5)
        length_adjusted_height = sea_sig_height * cm

        results.append([time_array[i], current_waveHs, cm, measured_sig_height, length_adjusted_height, sea_direction])

    # Determine the sea state based on the length adjusted PM SWH
    # Sea state SWH limits are based off of the MCTP SUROB table for ACV
    def sea_state_from_swh(length_adjusted_height):
        sea_state_limits_swh_ft = [0, 0.299, 0.99, 3.99, 7.99, 13.00]

        sea_state_swh = 0
        for i, limit in enumerate(sea_state_limits_swh_ft):
            if length_adjusted_height < limit:
                sea_state_swh = i
                break

        return sea_state_swh

    sea_state = sea_state_from_swh(length_adjusted_height)
    #print(f"Length Adjusted PM Sea state: {sea_state}") #Uncomment if output_df method is not used

    # Create a dataframe using the calculated values for each time step
    results_df = pd.DataFrame(results, columns=['Time', 'Buoy Wave Height (ft)', 'CM Normalization Factor', 'Calculated Total Height', 'Length Adjusted Height', 'Sea Direction'])

    # Create a dictionary to hold the output values
    output_dict = {
        "Length Adjusted PM Sea state": [sea_state],
        "Calculated PM Total Height": [round(measured_sig_height, 2)],
        "PM 3ft Length Adjusted Height": [round(length_adjusted_height, 2)]
    }

    # Create a dataframe to hold the output values
    #output_df = pd.DataFrame(output_dict) #Uncomment to generate table of values rather than have them displayed on plots

    # Display the output values dataframe
    #print(output_df)

    # Use for troubleshooting
    #print('CM Normalization Factor', cm)

    # Prints the 'Calculated Total Height', which is also the reported significant wave height form the buoy
    #print('Calculated PM Total Height', round(measured_sig_height, 2), 'ft') #Uncomment if output_df method is not used

    # Prints the 'Length Adjusted Height', which is the 3' normalized PM SWH
    #print('PM 3ft Length Adjusted Height', round(length_adjusted_height, 2), 'ft') #Uncomment if output_df method is not used

    # Plots the buoy reported SWH vs the calculated SWH vs the normalized SWH
    plt.figure(figsize=(12, 6))
    plt.plot(results_df['Time'], results_df['Buoy Wave Height (ft)'], label='Buoy Wave Height (ft)')
    #plt.plot(results_df['Time'], results_df['Calculated Total Height'], label='Calculated Buoy Height (ft)') #Use for troubleshooting
    plt.plot(results_df['Time'], results_df['Length Adjusted Height'], label='Length Adjusted Height (ft)')

    # Add forecasts for buoy and length adjusted
    plt.plot(df_forecast['Time'], df_forecast['Forecast Buoy Wave Height (ft)'], label='Forecast Buoy Wave Height (ft)', linestyle='--')
    plt.plot(df_forecast['Time'], df_forecast['Forecast LA SWH (ft)'], label='Forecast Length Adjusted Height (ft)', linestyle='--')

    plt.legend(loc='upper left')
    plt.title(f'Wave Heights over Time for {station_names.get(stn, "Unknown Station")}')
    plt.xlabel('Time')
    plt.ylabel('Wave Height (ft)')
    plt.xticks(rotation=45)

    # Add text on the plot
    plt.text(0.95, 0.09, f'Buoy Wave Height: {round(measured_sig_height, 2)} ft',
            verticalalignment='bottom', horizontalalignment='right',
            transform=plt.gca().transAxes,
            color='blue', fontsize=10)

    plt.text(0.95, 0.05, f'PM 3ft Length Adjusted Height: {round(length_adjusted_height, 2)} ft',
            verticalalignment='bottom', horizontalalignment='right',
            transform=plt.gca().transAxes,
            color='orange', fontsize=10)
    
    plt.text(0.95, 0.01, f'Length Adjusted PM Sea state: {sea_state}',
            verticalalignment='bottom', horizontalalignment='right',
            transform=plt.gca().transAxes,
            color='black', fontsize=10)

    plt.tight_layout()
    # Instead of plt.show(), use:
    pdf_pages.savefig()

pdf_pages.close()  # Close the PdfPages object
