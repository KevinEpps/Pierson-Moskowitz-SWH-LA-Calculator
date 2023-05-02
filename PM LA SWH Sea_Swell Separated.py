import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pytz
import cftime

# Download wave data from CDIP for specified time period and station ID
def download_wave_data(stn, start_time, end_time):
    data_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/{stn}p1_rt.nc'
    nc = netCDF4.Dataset(data_url)
    time_var = nc.variables['waveTime']
    
    # Find the nearest start and end index in the netCDF file based on the start_time and end_time provided
    start_index = netCDF4.date2index(start_time, time_var, select='nearest')
    end_index = netCDF4.date2index(end_time, time_var, select='nearest')
    
    # Extract necessary variables from netCDF file
    waveHs = nc.variables['waveHs'][start_index:end_index+1]
    waveTp = nc.variables['waveTp'][start_index:end_index+1]
    waveFrequency = nc.variables['waveFrequency'][:]
    waveEnergyDensity = nc.variables['waveEnergyDensity'][start_index:end_index+1,:]
    waveTz = nc.variables['waveTz'][start_index:end_index+1]

    # Convert time_var to a list of datetime objects in Pacific Standard Time (PST) and adjust for UTC offset
    time_array = [cftime.num2pydate(t, time_var.units).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('US/Pacific')) for t in time_var[start_index:end_index+1]]

    nc.close()

    return waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, time_array

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

# Download wave data
waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, time_array = download_wave_data(stn, start_time, end_time)

# Convert waveHs from meters to feet
waveHs_ft = waveHs * 3.28084

# Convert waveEnergyDensity from square meters to square feet
waveEnergyDensity_ft2 = waveEnergyDensity * 10.764

def remove_swell_energy(energy_density, frequency_range):
    modified_energy_density = np.copy(energy_density)
    
    # Determine the frequency range for which to remove swell energy
    freq_range_min = 0.025
    freq_range_max = np.min([0.100, np.max(frequency_range)])
    idx_closest_to_min_freq = np.argmin(np.abs(frequency_range - freq_range_min))
    idx_closest_to_max_freq = np.argmin(np.abs(frequency_range - freq_range_max))
    
    for t in range(energy_density.shape[0]):
        energy_density_0_100_hz = energy_density[t, idx_closest_to_min_freq:idx_closest_to_max_freq+1]
        min_energy_density_index = np.argmin(energy_density_0_100_hz)
        energy_density_min = energy_density_0_100_hz[min_energy_density_index]
        min_freq = frequency_range[idx_closest_to_min_freq+min_energy_density_index]
        for i, freq in enumerate(frequency_range):
            if freq_range_min <= freq <= freq_range_max:
                factor = 1 - energy_density_min * ((freq / freq_range_min) /  (freq_range_max - freq_range_min)) ** 8 / energy_density[t, i]
                factor = max(0, factor)
                modified_energy_density[t, i] = energy_density[t, i] * factor
    
    return modified_energy_density

def calculate_normalized_swh(energy_density, frequency_range, m1_ordinate_measured, m1_ordinate_desired):
    energy_density_no_swell = remove_swell_energy(energy_density, frequency_range)
    m1_ordinate_no_swell = calculate_m1(energy_density_no_swell, frequency_range)
    swh_no_swell = calculate_swh(energy_density_no_swell, frequency_range)
    normalized_swh = normalize_swh(swh_no_swell, m1_ordinate_measured, m1_ordinate_desired)
    return normalized_swh

def calculate_swh(waveEnergyDensity_ft2, waveFrequency):
    # Calculate the frequency intervals (delta_f) for the wave energy density data
    delta_f = np.diff(waveFrequency)
    delta_f = np.append(delta_f, delta_f[-1])  # Repeat the last delta_f for the last frequency

    # Integrate the wave energy density spectrum over the frequency range
    m0 = np.sum(waveEnergyDensity_ft2 * delta_f, axis=1)

    # Calculate the significant wave height (SWH)
    swh = 4 * np.sqrt(m0)

    return swh

swh_calculated_ft = calculate_swh(waveEnergyDensity_ft2, waveFrequency)

def calculate_m1(waveEnergyDensity_ft2, waveFrequency):
    # Calculate the frequency intervals (delta_f) for the wave energy density data
    delta_f = np.diff(waveFrequency)
    delta_f = np.append(delta_f, delta_f[-1])  # Repeat the last delta_f for the last frequency

    # Integrate the wave energy density spectrum over the frequency range
    m0 = np.sum(waveEnergyDensity_ft2 * delta_f, axis=1)

    # Calculate the M1 values by multiplying M0 values by their corresponding frequencies
    m1 = np.sum((waveEnergyDensity_ft2 * waveFrequency) * delta_f, axis=1)

    # Calculate the M1 ordinate by dividing the sum of the M1 values by the sum of the M0 ordinates
    m1_ordinate = np.sum(m1) / np.sum(m0)

    return m1_ordinate

m1_ordinate = calculate_m1(waveEnergyDensity_ft2, waveFrequency)

mean_period_m1 = 1/m1_ordinate
length_m1 = 5.12 * (mean_period_m1 ** 2)

def normalize_swh(swh_measured, m1_ordinate_measured, m1_ordinate_desired):
    # Calculate the average wave lengths for the measured and desired M1 ordinates
    # Wave Length calculation based on 32.2 ft/s^2 / 2PI = 5.12
    length_measured = 5.12 * (1/m1_ordinate_measured)**2
    length_desired = 5.12 * (1/m1_ordinate_desired)**2

    # Calculate the normalization factor
    normalization_factor = np.sqrt(length_desired / length_measured)

    # Normalize the significant wave height
    normalized_swh = swh_measured * normalization_factor

    return normalized_swh

# Calculate the M1 ordinate for the desired 3' PM SWH
# Mean period = 1/M1 ordinate.  # M1 = Sum of ft^2-Hz values / Sum M0 ft^2 values
mean_period_desired = 3.73  # Tp for 3' PM SWH - Based off of 71.18ft average length.  
length_desired = 5.12 * (mean_period_desired ** 2)
m1_ordinate_desired = 1 / mean_period_desired

# Remove the swell energy from the energy density
waveEnergyDensity_ft2_no_swell = remove_swell_energy(waveEnergyDensity_ft2, waveFrequency)

# Subtract modified_energy_density from waveEnergyDensity_ft2 to obtain the swell energy density
swell_energy_density = waveEnergyDensity_ft2 - waveEnergyDensity_ft2_no_swell

# Call remove_swell_energy function to obtain swell energy density
modified_energy_density = remove_swell_energy(waveEnergyDensity_ft2, waveFrequency)
swell_energy_density = waveEnergyDensity_ft2 - modified_energy_density

# Calculate modal period of sea energy density for the most recent row
freq_range_min = 0.025
freq_range_max = np.min([0.580, np.max(waveFrequency)])
idx_closest_to_min_freq = np.argmin(np.abs(waveFrequency - freq_range_min))
idx_closest_to_max_freq = np.argmin(np.abs(waveFrequency - freq_range_max))

sea_energy_density = modified_energy_density[-1, idx_closest_to_min_freq:idx_closest_to_max_freq+1]
max_energy_density_idx = np.argmax(sea_energy_density)
modal_period = 1 / waveFrequency[idx_closest_to_min_freq + max_energy_density_idx]

print(modal_period)

# Calculate the normalized significant wave height using only the sea energy density
normalized_swh_ft = calculate_normalized_swh(waveEnergyDensity_ft2, waveFrequency, m1_ordinate, m1_ordinate_desired)

def plot_swh(time_array, waveHs_ft, swh_calculated, normalized_swh):
    fig, ax = plt.subplots()

    # Plot the measured significant wave height
    ax.plot(time_array, waveHs_ft, label='Measured SWH (ft)', linestyle='-')

    # Plot the calculated significant wave height
    ax.plot(time_array, swh_calculated, label='Calculated SWH (ft)', linestyle='-')

    # Plot the normalized significant wave height
    ax.plot(time_array, normalized_swh, label='Normalized SWH (ft)', linestyle='-')

    # Set the date format for the x-axis
    date_format = mdates.DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(date_format)
    fig.autofmt_xdate()

    # Add labels, title, and legend
    ax.set_xlabel('Time')
    ax.set_ylabel('Significant Wave Height (ft)')
    ax.set_title('Measured vs. Calculated vs. Normalized Significant Wave Height')
    ax.legend()

    # Show the plot
    plt.show()

def plot_energy_density_spectrum(waveEnergyDensity_ft2, waveEnergyDensity_ft2_no_swell, swell_energy_density, waveFrequency, time_index):
    fig, ax = plt.subplots()

    # Plot the original energy density spectrum
    ax.plot(waveFrequency, waveEnergyDensity_ft2[time_index], label='Swell and Sea Energy Density', linestyle='--')

    # Plot the sea energy density spectrum
    ax.plot(waveFrequency, waveEnergyDensity_ft2_no_swell[time_index], label='Sea Energy Density', linestyle='-')

    # Plot the swell energy density spectrum
    ax.plot(waveFrequency, swell_energy_density[time_index], label='Swell Energy Density', linestyle='-')

    # Add labels, title, and legend
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Energy Density (ft^2/Hz)')
    ax.set_title('Energy Density Spectrum Comparison')
    ax.legend()

    # Show the plot
    plt.show()

def calculate_m0(waveEnergyDensity_ft2, waveFrequency):
    # Calculate the frequency intervals (delta_f) for the wave energy density data
    delta_f = np.diff(waveFrequency)
    delta_f = np.append(delta_f, delta_f[-1])  # Repeat the last delta_f for the last frequency

    # Integrate the wave energy density spectrum over the frequency range
    m0 = np.sum(waveEnergyDensity_ft2 * delta_f, axis=1)

    return m0

m0_original = calculate_m0(waveEnergyDensity_ft2, waveFrequency)
m0_modified = calculate_m0(waveEnergyDensity_ft2_no_swell, waveFrequency)

# Call the function to plot the data
plot_swh(time_array, waveHs_ft, swh_calculated_ft, normalized_swh_ft)

print('Significant Wave Height (ft)', np.round(waveHs_ft[-1], 2))
print('Estimated Pierson-Moskowitz Wave Height (ft)', np.round(swh_calculated_ft[-1], 2))
print('Normalized Sea Significant Wave Height (ft)', np.round(normalized_swh_ft[-1], 2))

plot_energy_density_spectrum(waveEnergyDensity_ft2, waveEnergyDensity_ft2_no_swell, swell_energy_density, waveFrequency, -1)

# Call remove_swell_energy function to obtain swell energy density
modified_energy_density = remove_swell_energy(waveEnergyDensity_ft2, waveFrequency)
swell_energy_density = waveEnergyDensity_ft2 - modified_energy_density

# Calculate modal period of sea energy density for the most recent row
freq_range_min = 0.025
freq_range_max = np.min([0.580, np.max(waveFrequency)])
idx_closest_to_min_freq = np.argmin(np.abs(waveFrequency - freq_range_min))
idx_closest_to_max_freq = np.argmin(np.abs(waveFrequency - freq_range_max))

sea_energy_density = modified_energy_density[-1, idx_closest_to_min_freq:idx_closest_to_max_freq+1]
max_energy_density_idx = np.argmax(sea_energy_density)
modal_period = 1 / waveFrequency[idx_closest_to_min_freq + max_energy_density_idx]

print(modal_period)