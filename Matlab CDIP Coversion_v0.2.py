import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pytz
import cftime
import pandas as pd

magnetic_declination = 11.53199; # From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14

# Set the buoy station ID and time range for data download
stn = '045'
end_time_utc = datetime.utcnow()
start_time_utc = end_time_utc - timedelta(hours=24)

# Set the timezones for the start and end times
utc_tz = pytz.timezone('UTC')
pst_tz = pytz.timezone('US/Pacific')

# Convert the start and end times to the Pacific timezone
end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

# Download wave data from CDIP for specified time period and station ID
data_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/{stn}p1_rt.nc'
nc = netCDF4.Dataset(data_url)
time_available_ser = nc.variables['waveTime']
    
# Find the nearest start and end index in the netCDF file based on the start_time and end_time provided
start_index = netCDF4.date2index(start_time, time_available_ser, select='nearest')
end_index = netCDF4.date2index(end_time, time_available_ser, select='nearest')

# Convert time_var to a list of datetime objects in Pacific Standard Time (PST) and adjust for UTC offset
time_available = [cftime.num2pydate(t, time_available_ser.units).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('US/Pacific')) for t in time_available_ser[start_index:end_index+1]]
    
# Extract necessary variables from netCDF file
significant_wave_height_feet = nc.variables['waveHs'][start_index:end_index+1] * 3.28084
peak_wave_period_seconds = nc.variables['waveTp'][start_index:end_index+1]
peak_wave_direction_deg_true = nc.variables['waveDp'][start_index:end_index+1]
average_wave_period_seconds = nc.variables['waveTa'][start_index:end_index+1]
wave_frequency = nc.variables['waveFrequency'][:]
wave_energy = nc.variables['waveEnergyDensity'][start_index:end_index+1,:]
wave_energy_raw = np.rot90(wave_energy)
wave_mean_direction_deg = nc.variables['waveMeanDirection'][start_index:end_index+1]
wave_mean_direction_deg_true = np.rot90(wave_mean_direction_deg)
wave_bandwidth = nc.variables['waveBandwidth'][:]
len_wavefreq = len(wave_frequency)

buoy_data_subset = {
    'time_available': time_available,
    'average_wave_period_seconds': average_wave_period_seconds,
    'peak_wave_direction_deg_true': peak_wave_direction_deg_true,
    'peak_wave_period_seconds': peak_wave_period_seconds,
    'significant_wave_height_feet': significant_wave_height_feet,
    'wave_energy_raw': wave_energy_raw,
    'wave_mean_direction_deg_true': wave_mean_direction_deg_true,
}

length_adjusted_table = pd.DataFrame()
sea_state_table = pd.DataFrame()
estimated_surf_height_table = pd.DataFrame()
estimated_surf_period_table = pd.DataFrame()

for process_number in range(len(buoy_data_subset)):

  #Calculate Wave Energy

  wave_energy_cut = wave_energy_raw
  wave_energy_raw = buoy_data_subset['wave_energy_raw']
  wave_energy = wave_energy_raw[:, process_number].reshape(-1, 1)

  #Create an array for bandwidth
  bandwidth = np.zeros((len_wavefreq, 1))

  for x in range(len_wavefreq-1, 1):
      bandwidth[x, 1] = (wave_frequency[x, 1] - wave_frequency[x-1, 1])/2 + \
                        (wave_frequency[x+1, 1] - wave_frequency[x, 1])/2

  # The first bandwidth point is equal to the second bandwidth point
  bandwidth[0, 0] = bandwidth[1, 0]

  # The last bandwidth point is equal the second to last bandwidth point.
  bandwidth[len_wavefreq-1, 0] = bandwidth[len_wavefreq-2, 0]

  # Calculate the Smoothed direction from the mean direction
  wave_mean_direction_deg_true = buoy_data_subset['wave_mean_direction_deg_true']
  wave_mean_direction_deg_true = wave_mean_direction_deg_true[:, process_number].reshape(-1, 1)

  smoothed_direction = np.zeros((len_wavefreq, 1))
  smoothed_direction[1, 0] = wave_mean_direction_deg_true[0, 0]/4 + \
                            wave_mean_direction_deg_true[1, 0]/2 + \
                            wave_mean_direction_deg_true[2, 0]/4

  for x in range(2, len_wavefreq-2):
      smoothed_direction[x, 0] = (wave_mean_direction_deg_true[x-2, 0] / 9
                                  + 2 * wave_mean_direction_deg_true[x-1, 0] / 9
                                  + 3 * wave_mean_direction_deg_true[x, 0] / 9
                                  + 2 * wave_mean_direction_deg_true[x+1, 0] / 9
                                  + wave_mean_direction_deg_true[x+2, 0] / 9)

  smoothed_direction[len_wavefreq-2, 0] = smoothed_direction[len_wavefreq-3, 0]
  smoothed_direction[len_wavefreq-1, 0] = smoothed_direction[len_wavefreq-3, 0]



  # For the first data point find the smoothed direction point by using interpolation
  smoothed_direction[0] = (smoothed_direction[2] - smoothed_direction[1]) * \
                            (wave_frequency[0] - wave_frequency[1]) / \
                            (wave_frequency[2] - wave_frequency[1]) + \
                            smoothed_direction[1]

  # Wave Period
  # Calculate the Wave Period. The wave period is equal to 1/f
  wave_period = 1/wave_frequency

  #Wave Length (5.12 = 32.2 ft/s^2 / 2PI)
  wave_length = 5.12 * wave_period ** 2
  
  # Total Energy Density    
  # This is done by multiplying the raw energy data by a factor of 10.76 to convert from m^2 to ft^2
  total_energy_density = wave_energy * 10.76391042

  # Calculate the Smoothed Total Energy Density.  Use the formula below to
  # calculate the Smoothed Total Energy Density from the Total Energy Density
  smoothed_total_energy_density = np.zeros((len_wavefreq, 1))

  for x in range(1, len_wavefreq-1):
    smoothed_total_energy_density[x, 0] = (total_energy_density[x-1, 0] / 4
                                            + total_energy_density[x, 0] / 2
                                            + total_energy_density[x+1, 0] / 4)

  smoothed_total_energy_density[len_wavefreq-1, 0] = smoothed_total_energy_density[len_wavefreq-2, 0]

  # Calculate the M0 Ordinate using the smooth total energy density. This is
  # done by multiplying the bandwidth and the smoothed total energy density together
  smoothed_total_energy_density_M0 = np.zeros((len_wavefreq, 1))

  for x in range(0, len_wavefreq):
      smoothed_total_energy_density_M0[x, 0] = bandwidth[x, 0] * smoothed_total_energy_density[x, 0]

  #Calculate the Smoothed Total Energy Density M1 Ordinate by multiplying the Smoothed 
  # M0 Ordinate by the frequency
  smoothed_total_energy_M1 = np.zeros((len_wavefreq, 1))

  for x in range(0, len_wavefreq):
      smoothed_total_energy_M1[x] = wave_frequency[x] * smoothed_total_energy_density_M0[x]
  
  # Calculate the smoothed slope energy density which is the difference
  # between the smoothed total energy density 
  smoothed_slope_energy_density = np.zeros((len_wavefreq, 1))

  for x in range(1, len_wavefreq):
      smoothed_slope_energy_density[x, 0] = smoothed_total_energy_density[x, 0] - smoothed_total_energy_density[x-1, 0]

  # Create the swell seperation Marker
  # Create a zero matrix for the swell seperation marker
  swell_seperation_marker = np.zeros((len_wavefreq, 1))

  # Create another zero matrix to find the minimum point in the Smoothed
  # Total Energy Density
  smoothed_total_energy_density_min_mat = np.zeros((6, 1))
  y = 0 #% y represents the row number in the smoothed total energy density and will be used in the for loop below

  # Find the minimum smoothed total energy density from the 11th - 16th points as the frequency range for swell energy.
  for x in range(10, 16):
    smoothed_total_energy_density_min_mat[y, 0] = smoothed_total_energy_density[x, 0]
    y += 1

  smoothed_total_energy_density_min = np.min(smoothed_total_energy_density_min_mat)

  # The following for loop is to be used to find the swell seperation marker at the minimum value for the smoothed total energy density.
  f_sted = None
  swell_seperation_marker_2 = None
  swell_seperation_marker = np.zeros((len_wavefreq, 1))

  for x in range(0, len_wavefreq):
      if np.any(smoothed_total_energy_density == smoothed_total_energy_density[x, 0]):
          swell_seperation_marker[x, 0] = 2
          f_sted = wave_frequency[x]
          swell_seperation_marker_2 = x

  # Calculate the Swell Seperation Index
  swell_seperation_index = np.zeros((len_wavefreq, 1))

  # The Swell separation index is calculated by adding the swell separation
  # marker to the previous swell separation index value multiplied by a
  # factor of 1.01
  for x in range(1, len_wavefreq):
      swell_seperation_index[x, 0] = swell_seperation_marker[x, 0] + swell_seperation_index[x-1, 0] * 1.01

  # Calculate the Sea Energy Density
  sea_energy_density = np.zeros((len_wavefreq, 1))
  swell_seperation_marker_sum = swell_seperation_marker[x, 0]  # Assumes x is defined and equals len_WaveFreq

  # The Sea Energy density equation is calculated using the frequency values
  # and the smoothed total energy density before the swell separation marker hits 2.0.
  for x in range(0, swell_seperation_marker_2):
      sea_energy_density[x] = smoothed_total_energy_density_min * ((wave_frequency[x] - wave_frequency[0]) / (f_sted - wave_frequency[0])) ** 8

  # After the swell separation marker hits 2.0, the sea energy density is equal to the smoothed total energy density.
  for x in range(swell_seperation_marker_2, len_wavefreq):
      sea_energy_density[x, 0] = smoothed_total_energy_density[x, 0]

  # Calculate the M0 Ordinate using the Sea Energy Density and the bandwidth
  sea_energy_density_M0 = np.zeros((len_wavefreq, 1))
  for x in range(0, len_wavefreq):
      sea_energy_density_M0[x, 0] = bandwidth[x, 0] * sea_energy_density[x, 0]

  # Calculate the M1 Ordinate using the M0 ordinate and the frequency
  sea_energy_density_M1 = np.zeros((len_wavefreq, 1))
  for x in range(0, len_wavefreq):
      sea_energy_density_M1[x] = wave_frequency[x] * sea_energy_density_M0[x]

  # Calculate the Swell Energy Density by taking the difference of the sea
  # energy density and the smoothed total energy density
  swell_energy_density = np.zeros((len_wavefreq, 1))
  for x in range(0, len_wavefreq):
      if smoothed_total_energy_density[x, 0] - sea_energy_density[x, 0] < 0:
          swell_energy_density[x, 0] = 0
      else:
          swell_energy_density[x, 0] = smoothed_total_energy_density[x, 0] - sea_energy_density[x, 0]

  # Calculate the M0 Ordinate using the Swell Energy Density and the bandwidth
  swell_energy_density_M0 = np.zeros((len_wavefreq, 1))
  for x in range(0, len_wavefreq):
      swell_energy_density_M0[x, 0] = bandwidth[x, 0] * swell_energy_density[x, 0]

  # Calculate the M1 Ordinate using the Swell Energy Density M0 ordinate and the frequency
  swell_energy_density_M1 = np.zeros((len_wavefreq, 1))
  for x in range(0, len_wavefreq):
      swell_energy_density_M1[x] = wave_frequency[x] * swell_energy_density_M0[x]

  # Calculate the Bretschneider & Pierson Moskowitz Energy Densities

  # Calculate the frequency in rad/sec
  f_rad = wave_frequency * 2 * np.pi

  # Convert bandwidth to rad/sec
  band_rad = bandwidth * 2 * np.pi

  # Calculate the H1/3 value, this value was not defined in the excel sheet
  H = 4 * np.sum(swell_energy_density_M0) ** (0.5)

  # Calculate the t-modal (wave period @ max swell_ED)
  # Find the max for the Swell Energy Density
  swell_max = np.max(swell_energy_density)
  for x in range(0, len_wavefreq):
      if swell_energy_density[x, 0] == Swell_max:
          t_modal_bs = wave_period[x]

  # Calculate the frequency in rad/sec
  f_rad = wave_frequency*2*np.pi

  # Convert bandwidth to rad/sec
  band_rad = bandwidth*2*np.pi

  # Calculate the H1/3 value
  H = 4*np.sum(swell_energy_density_M0)**0.5

  # Calculate the t-modal (wave period @ max swell_ED)
  # Find the max for the Swell Energy Density
  swell_max = np.max(swell_energy_density)
  t_modal_bs = wave_period[np.argmax(swell_energy_density)]

  # Calculate the W-modal for BretSchneider using the wave period
  w_modal_bs = 2*np.pi/t_modal_bs

  # Calculate the coefficients for the Bretschneider equation
  a_coef_bs = w_modal_bs**4*H**2*1.25/4
  b_coef_bs = -1.25*w_modal_bs**4

  # Calculate the BrettSchneider Energy Density using the BrettSchneider Equation
  bs = a_coef_bs/(f_rad+0.0001)**5*np.exp(b_coef_bs/(f_rad+0.0001)**4)

  # Calculate the M0 Ordinate using the Bretschneider Energy Density
  bs_M0 = bandwidth*bs

  # Calculate the M1 Ordinate using the smooth total energy density
  bs_M1 = wave_frequency*bs_M0

  # Calculate the H1/3 value.  The equation uses the sum of all the Sea
  # Energy Density M0 Ordinate values 

  H_pm = 4*np.sum(np.sqrt(sea_energy_density_M0))


  # Calculate the coefficients for the Pierson_Moskowitz equation

  a_coef_pm = 0.0081*(32.2**2) #32.2 is gravity in ft/s^2
  b_coef_pm = -0.032*(32.2/H_pm)**2 #32.2 is gravity in ft/s^2


  # Calculate the Pierson_Moskowitz Energy Density using the Pierson
  # Moskowitz equation

  pm = np.zeros(len_wavefreq)
  for x in range(len_wavefreq):
      pm[x] = a_coef_pm/(f_rad[x]+0.0001)**5*np.exp(b_coef_pm/(f_rad[x]+0.0001)**4)


  # Calculate the M0 Ordinate using the Pierson-Moskowitz Energy Density

  pm_M0 = np.zeros(len_wavefreq)
  for x in range(len_wavefreq):
      pm_M0[x] = bandwidth[x]*pm[x]


  # Calculate the M1 Ordinate using the Pierson-Moskowitz energy density

  pm_M1 = np.zeros(len_wavefreq)
  for x in range(len_wavefreq):
      pm_M1[x] = wave_frequency[x]*pm_M0[x]

  # Find the Target 3' Sig Pierson-Moskowitz

  a_coef = 0.0081 * 32.2 ** 2  # a_coef equation
  b_coef = -0.032 * (32.2 / 3) ** 2  # b_coef equation
  sigpm = np.zeros((len_wavefreq, 1))


  # Calculate the Target 3' Sig Pierson Moskowitz using the 3' Sig equation

  for x in range(len_wavefreq):
      sigpm[x] = a_coef / (f_rad[x] + 0.0001) ** 5 * np.exp(b_coef / (f_rad[x] + 0.0001) ** 4)


  # Calculate the data for swell height

  swell_height = 2 * (np.sum(swell_energy_density_M0) ** 0.5) / 0.707

  # Calculate the Swell Length

  swell_length = sum(swell_energy_density_M0) / sum(swell_energy_density_M1) # Divide the sums of the swell energy density M0 ordinate by the M1 Ordinate

  swell_mean_period = 1 / swell_length # Get the swell mean period
  swell_mean_period = 1 / swell_mean_period
  swell_length = 5.12 * swell_mean_period ** 2 # Use the equation to get the swell length

  # Calculate the Swell Direction by finding the mean direction at the max
  # swell energy density point

  swell_max = np.max(swell_energy_density)
  for x in range(len_wavefreq):
      if swell_energy_density[x,0] == swell_max:
          swell_direct = wave_mean_direction_deg_true[x,0]

  # Calculate the Swell Period
  swell_max = np.max(swell_energy_density)
  for x in range(len_wavefreq):
      if swell_energy_density[x,0] == swell_max:
          swell_period = wave_period[x]

  # Calculate the Swell 1/10th height
  swell_110_height = 5.1*sum(swell_energy_density_M0)**0.5

  # Calculate the swell 1/100 height
  swell_1100_height = 6.5*sum(swell_energy_density_M0)**0.5

  # Calculate the Estimated Surf Height
  estimated_surf_height = swell_height/0.645
  estimated_surf_height_str = str(round(estimated_surf_height, 1))

  # Calculate the estimated surf period
  estimated_surf_period = swell_period*0.585
  estimated_surf_period_str = str(round(estimated_surf_period, 1))

  # Calculate the Sea Sig Height
  sea_significant_height = ((sum(sea_energy_density_M0))**0.5) * 4

  # Calculate the Average length
  average_length = sum(sea_energy_density_M1)/sum(sea_energy_density_M0)
  average_length = 1/average_length
  average_length = 5.12*average_length**2

  # Calculate the sea modal time
  sea_energy_density_max = max(sea_energy_density)
  for x in range(len_wavefreq):
    if sea_energy_density[x, 0] == sea_energy_density_max:
      sea_modal_period = wave_period[x]

  # Calculate the Sea Direction  
  sea_energy_density_max = max(sea_energy_density)
  for x in range(len_wavefreq):
      if sea_energy_density[x, 0] == sea_energy_density_max:
          sea_direction_true = wave_mean_direction_deg_true[x, 0]

  sea_direction_true_str = str(round(sea_direction_true, 1))

  # Add the magnetic declination
  sea_direction_magnetic = sea_direction_true + magnetic_declination
  sea_direction_magnetic_str = str(round(sea_direction_magnetic, 1))

  #Calculate the Sea Mean Period
  sea_mean_period = sum(sea_energy_density_M0)
  sea_mean_period = sum(sea_energy_density_M1) / sea_mean_period
  sea_mean_period = 1 / sea_mean_period

  #Calculate the Sea 1/10 height.
  sea_110_height = 5.1 * sum(sea_energy_density_M0) ** 0.5

  #Calculate the sea 1/100 height
  sea_1100_height = 6.5 * sum(sea_energy_density_M0) ** 0.5

  #Calculate the total sig height using the Smoothed total energy density M0 ordinate
  calculated_total_significant_height = 4 * sum(smoothed_total_energy_density_M0) ** 0.5

  #Find the Wave Length Warning
  if sea_modal_period < 7:
    wave_length_warning = 'Synchronous Waves'
  else:
    wave_length_warning = 'none'

  # Find the Heading Warning
  if abs(sea_direction_true - swell_direct) < 40:
      heading_warning = 'Swell & Sea Same Dir'
  else:
      heading_warning = 'none'

  # Calculate the length adjusted PM Height
  length_adjusted_pm_height = sum(pm_M1) / sum(pm_M0)
  length_adjusted_pm_height = 1 / length_adjusted_pm_height
  length_adjusted_pm_height = 5.12 * length_adjusted_pm_height ** 2
  length_adjusted_pm_height = (length_adjusted_pm_height / average_length) ** 0.5

  if length_adjusted_pm_height > 1:
      length_adjusted_pm_height = 1

  if wave_length_warning.startswith('none'):
    length_adjusted_pm_height = length_adjusted_pm_height * sea_significant_height
  else:
    length_adjusted_pm_height = length_adjusted_pm_height * sea_significant_height * 1.05

  # If the heading warning is none then the length adjust pm height is
  # multiplied by 1.025

  if 'none' not in heading_warning.lower():
    length_adjusted_pm_height *= 1.025

  # Convert Length Adjusted PM to string for later use
  length_adjusted_pm_height_str = str(np.round(length_adjusted_pm_height, 2)) # Converts to string
  length_adjusted_pm_height_string = length_adjusted_pm_height_str

  # Find the Sea State Level
  # The sea state level is dependent on the length adjusted PM Height

  if 0.1 < length_adjusted_pm_height < 1:
      sea_state_level = 'Low '
  elif 1 < length_adjusted_pm_height < 1.8:
      sea_state_level = 'Mid '
  elif 1.8 < length_adjusted_pm_height < 2.1:
      sea_state_level = 'High '
  elif 2.1 < length_adjusted_pm_height < 3:
      sea_state_level = 'Low '
  elif 3 < length_adjusted_pm_height < 3.7:
      sea_state_level = 'Mid '
  elif 3.7 < length_adjusted_pm_height < 4.5:
      sea_state_level = 'High '
  elif 4.5 < length_adjusted_pm_height < 10:
      sea_state_level = 'Survival Range '
  else:
      sea_state_level = 'Do Not Operate '

  # Find the Sea State number
  # The sea state number is dependent on the significant total height and the
  # the length adjusted PM height
  if calculated_total_significant_height < 0.1:
      sea_state_level_number = '0'
  elif length_adjusted_pm_height < 0.7:
      sea_state_level_number = '1'
  elif 0.7 < length_adjusted_pm_height < 2.1:
      sea_state_level_number = '2'
  elif 2.1 < length_adjusted_pm_height < 4.5:
      sea_state_level_number = '3'
  elif 4.5 < length_adjusted_pm_height < 8.4:
      sea_state_level_number = '4'
  elif 8.4 < length_adjusted_pm_height < 13.8:
      sea_state_level_number = '5'
  elif 13.8 < length_adjusted_pm_height < 21.1:
      sea_state_level_number = '6'
  else:
      sea_state_level_number = 'Greater than 6'

  sea_state = sea_state_level + sea_state_level_number
  sea_state_short = str(sea_state)

  length_adjusted_table['Length_Adjusted_PM_Height'] = np.round(length_adjusted_pm_height, 1)
  sea_state_table['Sea_State'] = sea_state_short
  estimated_surf_height_table['Estimated_Surf_Height'] = np.round(estimated_surf_height, 1)
  estimated_surf_period_table['Estimated_Surf_Period'] = np.round(estimated_surf_period, 1)

  calculations_table = pd.concat([length_adjusted_table, sea_state_table, estimated_surf_height_table, estimated_surf_period_table], axis=1)

  buoy_data_subset.update(calculations_table) # Add to the subset table

  print(calculations_table)