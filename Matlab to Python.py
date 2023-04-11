 #Current WIP#

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

# set the buoy station number
stn = '045'
#stn = '264'

# set the date range
end_time_utc = datetime.utcnow()
start_time_utc = end_time_utc - timedelta(hours=12)
#start_time_utc = datetime(2023, 2, 28)
#end_time_utc = datetime(2023, 3, 6, 9)

# convert UTC time to PST time
utc_tz = pytz.timezone('UTC')
pst_tz = pytz.timezone('US/Pacific')
end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

# CDIP Realtime Dataset URL
address1 = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/' + stn + 'p1_rt.nc'
address2 = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_validation/BP' + stn + '_forecast.nc'
address3 = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/',stn,'p1/',stn,'p1_historic.nc'

# open the NetCDF dataset
nc = netCDF4.Dataset(address1)

# get the time variable
time_var = nc.variables['waveTime']

# find the indices for the time range
start_index = netCDF4.date2index(start_time, time_var, select='nearest')
end_index = netCDF4.date2index(end_time, time_var, select='nearest')

# get the necessary variables
hs_var = nc.variables['waveHs'][:] #Significant Wave Height in Meters
hsf_var = hs_var * 3.28084
t_var = nc.variables['waveTp'][:] #Peak wave period
fq_var = nc.variables['waveFrequency'][:] #Wave frequency
wave_energy = nc.variables['waveEnergyDensity'][:] #Wave energy density
waveDp_var = nc.variables['waveDp'][:] #Peak wave direction
waveDm_var = nc.variables['waveMeanDirection'][:] #Mean wave direction #not used - unable to figure out smoothed direction
waveDm_var = waveDm_var.T #Mean wave direction rotated 90 deg #not used - unable to figure out smoothed direction
waveTa_var = nc.variables['waveTa'][:] #Average wave period
sst_var = nc.variables['sstSeaSurfaceTemperature'][:] #Sea surface temperature
sstf_var = sst_var * 9/5 + 32 #Convert sea surface temperature to farenheit
band_var = nc.variables['waveBandwidth'][:] # Wave bandwidth
len_WaveFreq = len(nc.variables['waveFrequency'])

print("Sig Wave Height shape:", hsf_var.shape)
print("Peak Wave Period shape:", t_var.shape)
print("Wave frequency shape:", fq_var.shape)
print("Wave energy shape:", wave_energy.shape)
print("Wave bandwidth shape:", band_var.shape)
print("Wave mean direction shape:", waveDm_var.shape)

# extract the data for the time range
hs_data = hsf_var[start_index:end_index+1]
t_data = t_var[start_index:end_index+1]
waveDp_data = waveDp_var[start_index:end_index+1]
waveDm_data = waveDm_var[start_index:end_index+1] #not used - unable to figure out smoothed direction
waveDm_data = waveDm_var.T #rotate 90 deg #not used - unable to figure out smoothed direction
waveTa_data = waveTa_var[start_index:end_index+1]
ed_data = wave_energy[start_index:end_index+1][-1, :]
ed_data = ed_data.T #rotate 90 deg

# create a time array for plotting
time_array = netCDF4.num2date(time_var[start_index:end_index+1], time_var.units).tolist()

# convert the cftime objects to datetime objects
time_array_dt = [datetime.fromisoformat(str(t)) for t in time_array]

# convert datetime objects to strings
time_array_str = [t.strftime("%Y-%m-%d %H:%M:%S") for t in time_array]

# create a dictionary of the data columns
data = {
    "Wave Height": hs_data,
    "Wave Period": t_data,
    "Wave Direction": waveDp_data,
    "Avg Wave Period": waveTa_data,
    "Time": time_array_str
}

# create the dataframe
df = pd.DataFrame(data)

# set the "Time" column as the index
df.set_index("Time", inplace=True)

# Create an array for the Bandwidth
bandwidth = np.zeros(len_WaveFreq)

# Calculate a bandwidth using the formula below.
# The formula uses the midpoint formula to find the bandwidth
for x in range(1, len_WaveFreq - 1):
    bandwidth[x] = ((fq_var[x] - fq_var[x - 1]) / 2
                    + (fq_var[x + 1] - fq_var[x]) / 2)

# The first bandwidth point is equal to the second bandwidth point
bandwidth[0] = bandwidth[1]

# The last bandwidth point is equal the second to last bandwidth point
bandwidth[-1] = bandwidth[-2]

#Wave Period
#Calculate the Wave Period. The wave period is equal to 1/f

Wave_Period = fq_var ** (-1)

#Wave Length    

Wave_Length = 5.12 * Wave_Period ** 2

#Total Energy Density    
#This is done by multiplying the raw energy data by a factor of 10.76

Total_Energy_Density = wave_energy * 10.76391042

#Calculate the Smoothed Total Energy Density.  Use the formula below to
#calculate the Smoothed Total Energy Density from the Total Energy Density

Smoothed_Total_Energy_Density = np.zeros((len_WaveFreq, 1))
for x in range(1, len_WaveFreq-1):
    Smoothed_Total_Energy_Density[x, 0] = Total_Energy_Density[x-1, 0]/4 + Total_Energy_Density[x, 0]/2 + Total_Energy_Density[x+1, 0]/4
Smoothed_Total_Energy_Density[len_WaveFreq-1, 0] = Smoothed_Total_Energy_Density[len_WaveFreq-2, 0]

#Calculate the M0 Ordinate using the smooth total energy density. This is
#done by multiplying the bandwidth and the smoothed total energy density
#together

Smoothed_Total_Energy_Density_M0 = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    Smoothed_Total_Energy_Density_M0[x] = bandwidth[x]*Smoothed_Total_Energy_Density[x]

#Calculate the Smoothed Total Energy Density M1 Ordinate by multiplying the Smoothed 
#M0 Ordinate by the frequency
Smoothed_Total_Energy_M1 = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    Smoothed_Total_Energy_M1[x] = fq_var[x] * Smoothed_Total_Energy_Density_M0[x]


#STOPPED AT LINE 489 IN CDIP_BUOY_PROCESSOR#