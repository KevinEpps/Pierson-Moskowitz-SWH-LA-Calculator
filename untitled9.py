# -*- coding: utf-8 -*-
"""Untitled9.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1TWv62cSdki9K8yrIam6tiG0GqJcExFfo
"""

import numpy as np
import xarray as xr
import pandas as pd
from time import time
from tqdm import tqdm
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz

stn = input("Enter the station ID: Oceanside = 045, Red Beach = 264 ")

# set the date and time range for the latest 12 hour period
end_time_utc = datetime.utcnow()
start_time_utc = end_time_utc - timedelta(hours=12)

# convert UTC time to PST time
utc_tz = pytz.timezone('UTC')
pst_tz = pytz.timezone('US/Pacific')
end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

# set UTC to PST offset
utc_offset = (pst_tz.utcoffset(datetime.now()) - utc_tz.utcoffset(datetime.now())).seconds / 3600

# CDIP Realtime Dataset URL
address1 = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/' + stn + 'p1_rt.nc'

# open the NetCDF dataset
nc = netCDF4.Dataset(address1)

# get the time variable
time_var = nc.variables['waveTime']

# find the indices for the time range
start_index = netCDF4.date2index(start_time, time_var, select='nearest')
end_index = netCDF4.date2index(end_time, time_var, select='nearest')

# get the necessary variables
time_var = nc.variables['waveTime'][:]
hs_var = nc.variables['waveHs'][:] #Significant Wave Height in Meters
hsf_var = hs_var * 3.28084
t_var = nc.variables['waveTp'][:] #Peak wave period
fq_var = nc.variables['waveFrequency'][:] #Wave frequency
wave_energy = nc.variables['waveEnergyDensity'][:] #Wave energy density

# computes the mean of the wave_energy array across the second axis (axis=1), which is equivalent to taking the mean across all the columns (i.e., for each row). This results in a 1-dimensional array wave_energy_mean with length equal to the number of rows in the original wave_energy array. Each element in this new array represents the mean wave energy density for a single time step.
wave_energy_mean = np.mean(wave_energy, axis=1)

# selects the last column of the wave_energy array using indexing. This results in a 1-dimensional array wave_energy_last with length equal to the number of rows in the original wave_energy array. Each element in this new array represents the wave energy density for a single time step, but only at the last frequency band (i.e., the last column in the original wave_energy array).
wave_energy_last = wave_energy[:, -1]

#This will compute the mean wave energy across all 19531 rows for each of the 64 columns, resulting in an array of shape (64,)
mean_wave_energy = np.mean(wave_energy, axis=0)

waveDp_var = nc.variables['waveDp'][:] #Peak wave direction
waveDm_var = nc.variables['waveMeanDirection'][:] #Mean wave direction #not used - unable to figure out smoothed direction
waveDm_last = waveDm_var[:, -1]
waveDm_mean = np.mean(waveDm_var, axis=1)
#waveDm_var = waveDm_var.T #Mean wave direction rotated 90 deg #not used - unable to figure out smoothed direction
waveTa_var = nc.variables['waveTa'][:] #Average wave period
sst_var = nc.variables['sstSeaSurfaceTemperature'][:] #Sea surface temperature
sstf_var = sst_var * 9/5 + 32 #Convert sea surface temperature to farenheit
band_var = nc.variables['waveBandwidth'][:] # Wave bandwidth
len_WaveFreq = len(nc.variables['waveFrequency'])

data = {
    'Time': time_var,
    'Average Wave Period': waveTa_var,
    'Mean Wave Direction Last': waveDm_last,
    'Peak Wave Period': t_var,
    'Significant Wave Height (ft)': hsf_var,
    'Wave Energy Density Last': wave_energy_last,
    'Peak Wave Direction': waveDp_var
}

buoy_data = pd.DataFrame(data)

buoy_data_subset = buoy_data[start_index:end_index+1]

for index, row in buoy_data_subset.iterrows():
    significant_wave_height_feet_now = row['Significant Wave Height (ft)']
    significant_wave_height_feet_now_str = f"{round(significant_wave_height_feet_now, 2)}"

    peak_wave_period_seconds_now = row['Peak Wave Period']
    peak_wave_period_seconds_now_str = f"{round(peak_wave_period_seconds_now, 1)}"

    peak_wave_direction_deg_true_now = row['Peak Wave Direction']
    peak_wave_direction_deg_true_now_str = f"{round(peak_wave_direction_deg_true_now, 1)}"

    average_wave_period_seconds_now = row['Average Wave Period']

    wave_time_utc_now = row['Time']
    #wave_time_local_now = wave_time_utc_now - timedelta(hours=utc_offset)
    wave_time_local_now_str = end_time.strftime("%Y-%m-%d %H:%M:%S")

print(fq_var)
print(fq_var.shape)
print(mean_wave_energy)
print(mean_wave_energy.shape)