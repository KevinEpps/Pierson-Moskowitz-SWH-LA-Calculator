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
from math import exp

stn = input("Enter the station ID: Oceanside = 045, Red Beach = 264 ")

# set the date and time range for the latest 12 hour period
end_time_utc = datetime.utcnow()
start_time_utc = end_time_utc - timedelta(hours=12)

# From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14
Magnetic_Declination = 11.53199

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

# Create DataFrame if it doesn't exist
import pandas as pd

length_adjusted_table = pd.DataFrame()
Sea_State_table = pd.DataFrame()
Estimated_Surf_Height_table = pd.DataFrame()
Estimated_Surf_Period_table = pd.DataFrame()

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
print(buoy_data_subset.shape)
print(band_var.shape)

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

# Calculate the smoothed slope energy density which is the difference
# between the smoothed total energy density

Smoothed_Slope_Energy_Density = np.zeros(len_WaveFreq)
for x in range(1, len_WaveFreq):
    Smoothed_Slope_Energy_Density[x] = Smoothed_Total_Energy_Density[x] - Smoothed_Total_Energy_Density[x-1]

# Create a zero matrix for the swell separation marker
Swell_Separation_Marker = np.zeros((len_WaveFreq, 1))

# Create another zero matrix to find the minimum point in the Smoothed Total Energy Density
Smoothed_Total_Energy_Density_min_mat = np.zeros((6,1))

y = 0  # y represents the row number in the smoothed total energy density and will be used in the for loop below

# Create another zero matrix to find the minimum point in the Smoothed
# Total Energy Density

Smoothed_Total_Energy_Density_min_mat = np.zeros((6,1))

y = 0 # y represents the row number in the smoothed total energy density and will be used in the for loop below

# Find the minimum smoothed total energy density from the 11th - 16th
# points. This was how this was performed in the datawell buoy document
# and it is still being investigated as to why.

for x in range(10, 16):
    Smoothed_Total_Energy_Density_min_mat[y,0] = Smoothed_Total_Energy_Density[x,0]
    y += 1

Smoothed_Total_Energy_Density_min = np.min(Smoothed_Total_Energy_Density_min_mat)

# Create a zero matrix for the swell separation marker
Swell_Separation_Marker = np.zeros((len_WaveFreq, 1))

# Find the minimum smoothed total energy density from the 11th - 16th points
Smoothed_Total_Energy_Density_min_mat = Smoothed_Total_Energy_Density[10:16, :]
Smoothed_Total_Energy_Density_min = np.min(Smoothed_Total_Energy_Density_min_mat)

# Find the swell separation marker at the minimum value for the smoothed total energy density
for x in range(len_WaveFreq):
    if Smoothed_Total_Energy_Density_min == Smoothed_Total_Energy_Density[x, 0]:
        Swell_Separation_Marker[x, 0] = 2
        f_sted = fq_var[x]
        Swell_Separation_Marker_2 = x

# Calculate the Swell Separation Index
Swell_Separation_Index = np.zeros(len_WaveFreq)

# The Swell separation index is calculated by adding the swell separation
# marker to the previous swell separation index value multiplied by a
# factor of 1.01
for x in range(1, len_WaveFreq):
    Swell_Separation_Index[x] = Swell_Separation_Marker[x] + Swell_Separation_Index[x-1] * 1.01

# Calculate the Sea Energy Density
Sea_Energy_Density = np.zeros(len_WaveFreq)

#Swell_Separation_Marker_sum = Swell_Separation_Marker[x, 1]
Swell_Separation_Marker_sum = Swell_Separation_Marker[x][0]

# The Sea Energy density equation is calculated using the frequency values
# and the smoothed total energy density before the swell separation marker
# hits 2.0.

Sea_Energy_Density = np.zeros(len_WaveFreq)

for x in range(Swell_Separation_Marker_2):
    Sea_Energy_Density[x] = Smoothed_Total_Energy_Density_min * ((fq_var[x] - fq_var[0]) / (f_sted - fq_var[0])) ** 8

# After the swell separation marker hits 2.0, the sea energy density is
# equal to the smoothed total energy density
for x in range(Swell_Separation_Marker_2 + 1, len_WaveFreq):
    Sea_Energy_Density[x] = Smoothed_Total_Energy_Density[x]

# Calculate the M0 Ordinate using the Sea Energy Density and the bandwidth
Sea_Energy_Density_M0 = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    Sea_Energy_Density_M0[x] = bandwidth[x] * Sea_Energy_Density[x]

# Calculate the M1 Ordinate using the M0 ordinate and the frequency
Sea_Energy_Density_M1 = np.zeros((len_WaveFreq, 1))
for x in range(len_WaveFreq):
    Sea_Energy_Density_M1[x] = fq_var[x] * Sea_Energy_Density_M0[x]

# Calculate the Swell Energy Density by taking the difference of the sea
# energy density and the smoothed total energy density

Swell_Energy_Density = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    if Smoothed_Total_Energy_Density[x]-Sea_Energy_Density[x] < 0:
        Swell_Energy_Density[x] = 0
    else:
        Swell_Energy_Density[x] = Smoothed_Total_Energy_Density[x] - Sea_Energy_Density[x]

# Calculate the M0 Ordinate using the Swell Energy Density and the bandwidth
Swell_Energy_Density_M0 = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    Swell_Energy_Density_M0[x] = bandwidth[x] * Swell_Energy_Density[x]

# Calculate the M1 Ordinate using the swell energy density M0 ordinate and the frequency
Swell_Energy_Density_M1 = np.zeros(len_WaveFreq)
for x in range(len_WaveFreq):
    Swell_Energy_Density_M1[x] = fq_var[x] * Swell_Energy_Density_M0[x]

# Calculate the Bretschneider & Pierson Moskowitz Energy Densities

# Calculate the frequency in rad/sec
f_rad = fq_var * 2 * np.pi

# Convert bandwidth to rad/sec
band_rad = bandwidth * 2 * np.pi

# Calculate the H1/3 value, this value was not defined in the excel sheet
H = 4 * np.sum(Swell_Energy_Density_M0)**(.5)

# Calculate the t-modal (wave period @ max swell_ED)
# Find the max for the Swell Energy Density
Swell_max = np.max(Swell_Energy_Density)
for x in range(len_WaveFreq):
    if Swell_Energy_Density[x] == Swell_max:
        T_modal_BS = Wave_Period[x]

# Calculate the W-modal for BretSchneider using the wave period
W_modal_BS = 2 * np.pi / T_modal_BS

# Convert W-modal from rps to Hz
# W_modBS = W_modBS / (2 * np.pi)

# Calculate the coefficients for the Bretschneider equation
a_coef_BS = W_modal_BS ** 4 * H ** 2 * 1.25 / 4
b_coef_BS = -1.25 * W_modal_BS ** 4

# Calculate the BrettSchneider Energy Density using the BretSchneider Equation
BS = np.zeros((len_WaveFreq, 1))
for x in range(len_WaveFreq):
    BS[x] = a_coef_BS / (f_rad[x] + 0.0001) ** 5 * np.exp(b_coef_BS / (f_rad[x] + 0.0001) ** 4)

# Calculate the M0 Ordinate using the Bretschneider Energy Density
BS_M0 = np.zeros((len_WaveFreq, 1))
for x in range(len_WaveFreq):
    BS_M0[x] = bandwidth[x] * BS[x]

# Calculate the M1 Ordinate using the smooth total energy density
BS_M1 = np.zeros((len_WaveFreq,1))
for x in range(len_WaveFreq):
    BS_M1[x] = fq_var[x] * BS_M0[x]

# Calculate the H1/3 value. The equation uses the sum of all the Sea Energy Density M0 Ordinate values
H_PM = 4 * np.sum(Sea_Energy_Density_M0) ** 0.5

# Calculate the coefficients for the Pierson_Moskowitz equation
a_coef_PM = 0.0081 * 32.2 ** 2  # 32.2 is gravity in ft/s^2
b_coef_PM = -0.032 * (32.2 / H_PM) ** 2  # 32.2 is gravity in ft/s^2

# Calculate the Pierson_Moskowitz Energy Density using the Pierson Moskowitz equation
PM = np.zeros((len_WaveFreq, 1))
for x in range(len_WaveFreq):
    PM[x] = a_coef_PM / (f_rad[x] + 0.0001) ** 5 * np.exp(b_coef_PM / (f_rad[x] + 0.0001) ** 4)

# Calculate the M0 Ordinate using the Pierson-Moskowitz Energy Density
PM_M0 = np.zeros((len_WaveFreq,1))
for x in range(len_WaveFreq):
    PM_M0[x] = bandwidth[x] * PM[x]

# Calculate the M1 Ordinate using the Pierson-Moskowitz energy density
PM_M1 = np.zeros((len_WaveFreq,1))
for x in range(len_WaveFreq):
    PM_M1[x] = fq_var[x] * PM_M0[x]

# Find the Target 3' Sig Pierson-Moskowitz
a_coef = 0.0081*32.2**2  # a_coef equation
b_coef = -0.032*(32.2/3)**2  # b_coef equation
sigPM = np.zeros((len_WaveFreq,1))

# Calculate the Target 3' Sig Pierson Moskowitz using the 3' Sig equation
for x in range(len_WaveFreq):
    sigPM[x] = a_coef / (f_rad[x] + 0.0001)**5 * exp(b_coef / (f_rad[x] + 0.0001)**4)

# Calculate the data for swell height
Swell_Height = 2*(sum(Swell_Energy_Density_M0)**0.5) / 0.707

# Calculate the Swell Length
Swell_Length = sum(Swell_Energy_Density_M0) / sum(Swell_Energy_Density_M1) # Divide the sums of the swell energy density M0 ordinate by the M1 Ordinate
Swell_Mean_Period = 1 / Swell_Length # Get the swell mean period
Swell_Mean_Period = 1 / Swell_Mean_Period
Swell_Length = 5.12 * Swell_Mean_Period**2 # Use the equation to get the swell length

# Calculate the Swell Direction by finding the mean direction at the max swell energy density point
Swell_max = max(Swell_Energy_Density)
for x in range(len_WaveFreq):
    if Swell_Energy_Density[x] == Swell_max:
        Swell_Direct = waveDm_var[x]

# Calculate the Swell Period
Swell_max = max(Swell_Energy_Density)
for x in range(len_WaveFreq):
    if Swell_Energy_Density[x] == Swell_max:
        Swell_Period = Wave_Period[x]

# Calculate the Swell 1/10th height
Swell_110_Height = 5.1 * sum(Swell_Energy_Density_M0) ** 0.5

# Calculate the swell 1/100 height 
Swell_1100_Height = 6.5 * sum(Swell_Energy_Density_M0) ** 0.5

# Calculate the Estimated Surf Height
Estimated_Surf_Height = Swell_Height / 0.645
Estimated_Surf_Height_str = str(round(Estimated_Surf_Height, 1)) # Converts to string

# Calculate the estimated surf period
Estimated_Surf_Period = Swell_Period * 0.585
Estimated_Surf_Period_str = str(round(Estimated_Surf_Period, 1)) # Converts to string

# Caluclate the Sea Sig Height
Sea_Significant_Height = (sum(Sea_Energy_Density_M0) ** 0.5) * 4

# Calculate the Average length
Average_Length = sum(Sea_Energy_Density_M1) / sum(Sea_Energy_Density_M0)
Average_Length = 1 / Average_Length
Average_Length = 5.12 * Average_Length ** 2

# Calculate the sea modal time
Sea_Energy_Density_max = max(Sea_Energy_Density)
for x in range(len_WaveFreq):
    if Sea_Energy_Density[x] == Sea_Energy_Density_max:
        Sea_Modal_Period = Wave_Period[x]

# Calculate the Sea Direction
Sea_Energy_Density_max = max(Sea_Energy_Density)
for x in range(len_WaveFreq):
    if Sea_Energy_Density[x] == Sea_Energy_Density_max:
        Sea_Direction_True = waveDm_var[x]
Sea_Direction_True_str = str((Sea_Direction_True, 1))

# Add the magnetic declination

Sea_Direction_Magnetic = Sea_Direction_True + Magnetic_Declination
Sea_Direction_Magnetic_str = str((Sea_Direction_Magnetic, 1))


# Calculate the Sea Mean Period

Sea_Mean_Period = sum(Sea_Energy_Density_M0)
Sea_Mean_Period = sum(Sea_Energy_Density_M1) / Sea_Mean_Period
Sea_Mean_Period = 1 / Sea_Mean_Period


# Calculate the Sea 1/10 height.

Sea_110_Height = 5.1 * sum(Sea_Energy_Density_M0) ** 0.5


# Calculate the sea 1/100 height 

Sea_1100_Height = 6.5 * sum(Sea_Energy_Density_M0) ** 0.5

# Calculate the total sig height using the Smoothed total energy density M0 ordinate
Calculated_Total_Significant_Height = 4 * sum(Smoothed_Total_Energy_Density_M0) ** 0.5

# Find the Wave Length Warning
if Sea_Modal_Period < 7:
    Wave_Length_Warning = 'Synchronous Waves'
else:
    Wave_Length_Warning = 'none'

# Find the Heading Warning
if any(Sea_Direction_True - Swell_Direct) < 40:
    Heading_Warning = 'Swell & Sea Same Dir'
else:
    Heading_Warning = 'none'

# Calculate the length adjusted PM Height

Length_Adjusted_PM_Height = sum(PM_M1) / sum(PM_M0)
Length_Adjusted_PM_Height = 1 / Length_Adjusted_PM_Height
Length_Adjusted_PM_Height = 5.12 * Length_Adjusted_PM_Height**2  # These factors were taken from the excel sheet DataWellV08
Length_Adjusted_PM_Height = (Length_Adjusted_PM_Height / Average_Length)**0.5

if Length_Adjusted_PM_Height > 1:
    Length_Adjusted_PM_Height = 1

if 'none' in Wave_Length_Warning.lower():  # if the wave length warning is none then the following equation is to be used to calculate the length adjusted PM height
    # Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Calculated_Total_Significant_Height
    Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Sea_Significant_Height
else:
    # Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Calculated_Total_Significant_Height * 1.05
    Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Sea_Significant_Height * 1.05

# If the heading warning is none then the length adjust pm height is
# multiplied by 1.025
if Heading_Warning.lower() != 'none':
    Length_Adjusted_PM_Height *= 1.025

# Convert Length Adjusted PM to string for later use
Length_Adjusted_PM_Height_str = str((Length_Adjusted_PM_Height, 2))

# Find the Sea State Level
# The sea state level is dependent on the length adjusted PM Height
if 0.1 < Length_Adjusted_PM_Height < 1:
    Sea_State_Level = 'Low '
elif 1 < Length_Adjusted_PM_Height < 1.8:
    Sea_State_Level = 'Mid '
elif 1.8 < Length_Adjusted_PM_Height < 2.1:
    Sea_State_Level = 'High '
elif 2.1 < Length_Adjusted_PM_Height < 3:
    Sea_State_Level = 'Low '
elif 3 < Length_Adjusted_PM_Height < 3.7:
    Sea_State_Level = 'Mid '
elif 3.7 < Length_Adjusted_PM_Height < 4.5:
    Sea_State_Level = 'High '
elif 4.5 < Length_Adjusted_PM_Height < 10:
    Sea_State_Level = 'Survival Range '
else:
    Sea_State_Level = 'Do Not Operate '

# Find the Sea State number
# The sea state number is dependent on the significant total height and the
# the length adjusted PM height

if Calculated_Total_Significant_Height < 0.1:
    Sea_State_Level_Number = '0'
elif Length_Adjusted_PM_Height < 0.7:
    Sea_State_Level_Number = '1'
elif (0.7 < Length_Adjusted_PM_Height) and (Length_Adjusted_PM_Height < 2.1):
    Sea_State_Level_Number = '2'
elif (2.1 < Length_Adjusted_PM_Height) and (Length_Adjusted_PM_Height < 4.5):
    Sea_State_Level_Number = '3'
elif (4.5 < Length_Adjusted_PM_Height) and (Length_Adjusted_PM_Height < 8.4):
    Sea_State_Level_Number = '4'
elif (8.4 < Length_Adjusted_PM_Height) and (Length_Adjusted_PM_Height < 13.8):
    Sea_State_Level_Number = '5'
elif (13.8 < Length_Adjusted_PM_Height) and (Length_Adjusted_PM_Height < 21.1):
    Sea_State_Level_Number = '6'
else:
    Sea_State_Level_Number = 'Greater than 6'

Sea_State = Sea_State_Level + Sea_State_Level_Number
Sea_State_short = str(Sea_State)

length_adjusted_table = length_adjusted_table.append({'Length_Adjusted_PM_Height': (Length_Adjusted_PM_Height,1)}, ignore_index=True)
Sea_State_table = Sea_State_table.append({'Sea_State': Sea_State_short}, ignore_index=True)
Estimated_Surf_Height_table = Estimated_Surf_Height_table.append({'Estimated_Surf_Height': (Estimated_Surf_Height,1)}, ignore_index=True)
Estimated_Surf_Period_table = Estimated_Surf_Period_table.append({'Estimated_Surf_Period': (Estimated_Surf_Period,1)}, ignore_index=True)


calculations_table = pd.concat([length_adjusted_table, Sea_State_table, Estimated_Surf_Height_table, Estimated_Surf_Period_table], axis=1)

buoy_data_subset = pd.concat([buoy_data_subset, calculations_table], axis=1)

print(calculations_table)
print(buoy_data_subset)
