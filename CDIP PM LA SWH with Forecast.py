import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz

# set the buoy station number
stn = '045'
#stn = '264'

# set the length of the wave period for length adjustment
#T = 12

# set the date range
#end_time_utc = datetime.utcnow()
#start_time_utc = end_time_utc - timedelta(hours=12)
start_time_utc = datetime(2023, 2, 28)
end_time_utc = datetime(2023, 3, 6, 9)

# convert UTC time to PST time
utc_tz = pytz.timezone('UTC')
pst_tz = pytz.timezone('US/Pacific')
end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

# CDIP Realtime Dataset URL
data_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/' + stn + 'p1_rt.nc'
data_url_forecast = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_validation/BP' + stn + '_forecast.nc'

# open the NetCDF dataset
nc = netCDF4.Dataset(data_url)

# open the forecast NetCDF dataset
nc_f = netCDF4.Dataset(data_url_forecast)

# get the time variable
time_var = nc.variables['waveTime']

# get the forecast time variable
time_var_f = nc_f.variables['waveTime']

# find the indices for the time range
start_index = netCDF4.date2index(start_time, time_var, select='nearest')
end_index = netCDF4.date2index(end_time, time_var, select='nearest')

# find the forecast indices for the time range
start_index_forecast = netCDF4.date2index(start_time, time_var_f, select='nearest')
end_index_forecast = netCDF4.date2index(end_time, time_var_f, select='nearest')

# get the wave height, period, mean direction, frequency, and energy density variables
hs_var = nc.variables['waveHs']
t_var = nc.variables['waveTp']
dmean_var = nc.variables['waveMeanDirection']
fq_var = nc.variables['waveFrequency']
ed_var = nc.variables['waveEnergyDensity']

# get the forecast wave height, period, mean direction, frequency, and energy density variables
hs_var_f = nc_f.variables['waveHs']
t_var_f = nc_f.variables['waveTp']
dmean_var_f = nc_f.variables['waveMeanDirection']
fq_var_f = nc_f.variables['waveFrequency']
ed_var_f = nc_f.variables['waveEnergyDensity']

# extract the wave height and period data for the time range
hs_data = hs_var[start_index:end_index+1]
t_data = t_var[start_index:end_index+1]
dmean_data = dmean_var[start_index:end_index+1]
fq_data = fq_var[:]
ed_data = ed_var[start_index:end_index+1,:]
T = np.mean(t_data)

# extract the forecast wave height and period data for the time range
hs_data_f = hs_var_f[start_index_forecast:end_index_forecast+1]
t_data_f = t_var_f[start_index_forecast:end_index_forecast+1]
dmean_data_f = dmean_var_f[start_index_forecast:end_index_forecast+1]
fq_data_f = fq_var_f[:]
ed_data_f = ed_var_f[start_index_forecast:end_index_forecast+1,:]
T_f = np.mean(t_data_f)

# calculate length adjusted wave height using the Pierson-Moskowitz formula
g = 9.81  # acceleration due to gravity
fp = 1 / T  # peak frequency
spectrum = (5/16) * (hs_data**2) * (fp**4) / (g**2)
lah = np.sqrt(2 * spectrum / fp) * 200  # adjust by a factor of 200

# calculate forecast length adjusted wave height using the Pierson-Moskowitz formula
g = 9.81  # acceleration due to gravity
fp_f = 1 / T_f  # peak frequency
spectrum_f = (5/16) * (hs_data_f**2) * (fp_f**4) / (g**2)
lah_f = np.sqrt(2 * spectrum_f / fp_f) * 200  # adjust by a factor of 200

# convert from meters to feet
lah = lah * 3.281

# convert from meters to feet (forecast)
lah_f = lah_f * 3.281

# create a time array for plotting
time_array = netCDF4.num2date(time_var[start_index:end_index+1], time_var.units).tolist()

# create a forecast time array for plotting
time_array_f = netCDF4.num2date(time_var_f[start_index_forecast:end_index_forecast+1], time_var_f.units).tolist()

# convert the cftime objects to datetime objects and PST timezone
time_array = [datetime.fromisoformat(str(t)).astimezone(pst_tz) for t in time_array]

# convert the forecast cftime objects to datetime objects and PST timezone
time_array_f = [datetime.fromisoformat(str(t)).astimezone(pst_tz) for t in time_array_f]

# adjust the time offset by 6 hours
time_array = [t - timedelta(hours=6) for t in time_array]

# adjust the forecast time offset by 6 hours
time_array_f = [t - timedelta(hours=6) for t in time_array_f]

import matplotlib.pyplot as plt

# create the plot
fig, ax = plt.subplots()

# plot the data and forecast on the same axis object
ax.plot(time_array, lah, label='Length-Adjusted Wave Height')
ax.plot(time_array_f, lah_f, label='Forecasted Length-Adjusted Wave Height')

# format the plot
plt.title('Length-Adjusted Wave Height for Buoy ' + stn)
plt.xlabel('Time (PST)')
plt.ylabel('Wave Height (ft)')
plt.grid()

# set the x-axis format to show date and time in hh:mm format
date_format = '%m/%d/%Y %H:%M'
date_formatter = plt.matplotlib.dates.DateFormatter(date_format)
ax.xaxis.set_major_formatter(date_formatter)
fig.autofmt_xdate()

# add a legend to the plot
plt.legend()

# display the plot
plt.show()

# close the NetCDF dataset
nc.close()
