import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pytz
import cftime
import ipywidgets as widgets


def download_wave_data(stn, start_time, end_time):
    """
    Download wave data from CDIP for the specified time period and station ID.

    :param stn: int, station ID
    :param start_time: datetime, start time for the data extraction
    :param end_time: datetime, end time for the data extraction
    :return: tuple containing wave data variables and time array
    """
    
    # Construct the data URL for the specified station
    data_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/{stn}p1_rt.nc'
    
    # Open the netCDF file at the specified URL
    nc = netCDF4.Dataset(data_url)
    
    # Get the time variable from the netCDF file
    time_var = nc.variables['waveTime']
    
    # Find the nearest start and end index in the netCDF file based on the start_time and end_time provided
    start_index = netCDF4.date2index(start_time, time_var, select='nearest')
    end_index = netCDF4.date2index(end_time, time_var, select='nearest')
    
    # Extract necessary variables from the netCDF file
    waveHs = nc.variables['waveHs'][start_index:end_index+1]
    waveTp = nc.variables['waveTp'][start_index:end_index+1]
    waveFrequency = nc.variables['waveFrequency'][:]
    waveEnergyDensity = nc.variables['waveEnergyDensity'][start_index:end_index+1,:]
    waveTz = nc.variables['waveTz'][start_index:end_index+1]
    waveDirection = nc.variables['waveMeanDirection'][start_index:end_index+1]

    # Convert time_var to a list of datetime objects in Pacific Standard Time (PST) and adjust for UTC offset
    time_array = [cftime.num2pydate(t, time_var.units).replace(tzinfo=pytz.utc).astimezone(pytz.timezone('US/Pacific')) for t in time_var[start_index:end_index+1]]
    
    # Close the netCDF file
    nc.close()

    return waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, waveDirection, time_array

def remove_swell_energy(energy_density, frequency_range):
    """
    Remove the swell energy from the wave energy density spectrum.

    :param energy_density: numpy array, energy density data
    :param frequency_range: numpy array, frequency range data
    :return: numpy array, modified energy density with swell energy removed
    """

    # Create a copy of the input energy_density array to store the modified energy density
    modified_energy_density = np.copy(energy_density)

    # Define the frequency range for which to remove swell energy
    freq_range_min = 0.025
    freq_range_max = np.min([0.100, np.max(frequency_range)])

    # Find the indices of the frequency range in the input data
    idx_closest_to_min_freq = np.argmin(np.abs(frequency_range - freq_range_min))
    idx_closest_to_max_freq = np.argmin(np.abs(frequency_range - freq_range_max))

    # Loop through each time step in the energy_density array
    for t in range(energy_density.shape[0]):
        # Extract the energy_density values within the defined frequency range
        energy_density_0_100_hz = energy_density[t, idx_closest_to_min_freq:idx_closest_to_max_freq+1]

        # Find the index of the minimum energy density within the extracted values
        min_energy_density_index = np.argmin(energy_density_0_100_hz)

        # Get the minimum energy density value and its corresponding frequency
        energy_density_min = energy_density_0_100_hz[min_energy_density_index]
        min_freq = frequency_range[idx_closest_to_min_freq+min_energy_density_index]

        # Loop through the frequency_range and apply the swell energy removal factor
        for i, freq in enumerate(frequency_range):
            if freq_range_min <= freq <= freq_range_max:
                factor = 1 - energy_density_min * ((freq / freq_range_min) /  (freq_range_max - freq_range_min)) ** 8 / energy_density[t, i]
                factor = max(0, factor)
                modified_energy_density[t, i] = energy_density[t, i] * factor

    return modified_energy_density


def apply_modal_period_correction(normalized_swh, modal_period):
    """
    Apply a correction factor to the normalized significant wave height (swh) based on the modal period.

    :param normalized_swh: float, the normalized significant wave height
    :param modal_period: float, the modal period
    :return: float, the corrected significant wave height
    """
    if modal_period < 7:
        return normalized_swh * 1.05
    else:
        return normalized_swh * 1


def calculate_correction_factor(waveDirection, modified_energy_density, swell_energy_density, modal_swh):
    """
    Calculate the correction factor for significant wave height based on wave directions and energy density.

    :param waveDirection: numpy array, wave direction data
    :param modified_energy_density: numpy array, modified energy density data (sea)
    :param swell_energy_density: numpy array, swell energy density data
    :param modal_swh: float, the modal significant wave height
    :return: float, the corrected significant wave height
    """

    # Define magnetic declination constant from NOAA
    magnetic_declination = 11.53199  # From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14

    # Find direction corresponding to the highest swell energy density value
    max_swell_energy_density_index = np.argmax(np.sum(swell_energy_density, axis=1))
    swell_direction = waveDirection[max_swell_energy_density_index] + magnetic_declination

    # Find direction corresponding to the highest sea energy density value
    max_sea_energy_density_index = np.argmax(np.sum(modified_energy_density, axis=1))
    sea_direction = waveDirection[max_sea_energy_density_index] + magnetic_declination

    # Calculate correction factor
    if (swell_direction - sea_direction) < 40:
        return modal_swh * 1.025
    else:
        return modal_swh * 1


def calculate_normalized_swh(energy_density, frequency_range, m1_ordinate_measured, m1_ordinate_desired):
    """
    Calculate the normalized significant wave height (SWH) using energy density data, frequency range,
    and measured and desired M1 ordinates.

    :param energy_density: numpy array, wave energy density data
    :param frequency_range: numpy array, wave frequency data
    :param m1_ordinate_measured: float, measured M1 ordinate value
    :param m1_ordinate_desired: float, desired M1 ordinate value
    :return: float, the normalized significant wave height (SWH)
    """

    # Remove swell energy from energy density data
    energy_density_no_swell = remove_swell_energy(energy_density, frequency_range)

    # Calculate the M1 ordinate with no swell
    m1_ordinate_no_swell = calculate_m1(energy_density_no_swell, frequency_range)

    # Calculate the significant wave height (SWH) with no swell
    swh_no_swell = calculate_swh(energy_density_no_swell, frequency_range)

    # Normalize the SWH based on the measured and desired M1 ordinates
    normalized = normalize_swh(swh_no_swell, m1_ordinate_measured, m1_ordinate_desired)

    # Apply a modal period correction to the normalized SWH
    modal_swh = apply_modal_period_correction(normalized, modal_period)

    # Calculate the normalized SWH using the correction factor
    normalized_swh = calculate_correction_factor(waveDirection_true, modified_energy_density, swell_energy_density, modal_swh)

    return normalized_swh


def calculate_swh(waveEnergyDensity_ft2, waveFrequency):
    """
    Calculate the significant wave height (SWH) from wave energy density and frequency data.

    :param waveEnergyDensity_ft2: numpy array, wave energy density data in ft^2
    :param waveFrequency: numpy array, wave frequency data in Hz
    :return: numpy array, significant wave height (SWH) in ft
    """

    # Calculate the frequency intervals (delta_f) for the wave energy density data
    delta_f = np.diff(waveFrequency)
    delta_f = np.append(delta_f, delta_f[-1])  # Repeat the last delta_f for the last frequency

    # Integrate the wave energy density spectrum over the frequency range
    m0 = np.sum(waveEnergyDensity_ft2 * delta_f, axis=1)

    # Calculate the significant wave height (SWH)
    swh = 4 * np.sqrt(m0)

    return swh


def calculate_m1(waveEnergyDensity_ft2, waveFrequency):
    """
    Calculate the M1 ordinate from wave energy density and frequency data.

    :param waveEnergyDensity_ft2: numpy array, wave energy density data in ft^2
    :param waveFrequency: numpy array, wave frequency data in Hz
    :return: float, M1 ordinate value
    """

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


def normalize_swh(swh_measured, m1_ordinate_measured, m1_ordinate_desired):
    """
    Normalize the significant wave height (SWH) based on measured and desired M1 ordinates.

    :param swh_measured: numpy array, measured significant wave height in ft
    :param m1_ordinate_measured: float, measured M1 ordinate value
    :param m1_ordinate_desired: float, desired M1 ordinate value
    :return: numpy array, normalized significant wave height in ft
    """

    # Calculate the average wave lengths for the measured and desired M1 ordinates
    # Wave Length calculation based on 32.2 ft/s^2 / 2PI = 5.12
    length_measured = 5.12 * (1/m1_ordinate_measured)**2
    length_desired = 5.12 * (1/m1_ordinate_desired)**2

    # Calculate the normalization factor
    normalization_factor = np.sqrt(length_desired / length_measured)

    # Normalize the significant wave height
    normalized_swh = swh_measured * normalization_factor

    return normalized_swh

def plot_swh(time_array, waveHs_ft, swh_calculated, normalized_swh):
    """
    Plot measured, calculated, and normalized significant wave height (SWH) over time.

    :param time_array: list of datetime objects, time values for data points
    :param waveHs_ft: numpy array, measured significant wave height in ft
    :param swh_calculated: numpy array, calculated significant wave height in ft
    :param normalized_swh: numpy array, normalized significant wave height in ft
    """
    fig, ax = plt.subplots()

    # Plot the measured significant wave height
    ax.plot(time_array, waveHs_ft, label='Measured Buoy Significant Wave Height (ft)', linestyle='-')

    # Plot the calculated significant wave height
    ax.plot(time_array, swh_calculated, label='Estimated Pierson-Moskowitz Wave Height (ft)', linestyle='-')

    # Plot the normalized significant wave height
    ax.plot(time_array, normalized_swh, label='3ft Length-Adjusted Significant Wave Height (ft)', linestyle='-')

    # Set the date format for the x-axis
    date_format = mdates.DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(date_format)
    fig.autofmt_xdate()

    # Add labels, title, and legend
    ax.set_xlabel('Time')
    ax.set_ylabel('Significant Wave Height (ft)')
    ax.set_title('Measured vs. Calculated vs. Normalized Significant Wave Height')
    ax.legend()

    # Create the table data
    table_data = [
        ["Measured Buoy Significant Wave Height (ft)", f"{waveHs_ft[-1]:.2f}"],
        ["Estimated Pierson-Moskowitz Wave Height (ft)", f"{swh_calculated_ft[-1]:.2f}"],
        ["3ft Length-Adjusted Significant Wave Height (ft)", f"{normalized_swh_ft[-1]:.2f}"]
    ]

    # Create the table
    table = plt.table(cellText=table_data, colWidths=[0.5, 0.3], cellLoc='left',
                      loc='bottom', bbox=[0.1, -0.65, 0.8, 0.25])

    # Adjust the plot layout to make room for the table
    plt.subplots_adjust(bottom=0.25)

    # Show the plot
    plt.show()

def plot_energy_density_spectrum(waveEnergyDensity_ft2, waveEnergyDensity_ft2_no_swell, swell_energy_density, waveFrequency, time_index):
    """
    Plot the energy density spectrum comparison between swell and sea, sea only, and swell only.

    :param waveEnergyDensity_ft2: numpy array, swell and sea energy density in ft^2/Hz
    :param waveEnergyDensity_ft2_no_swell: numpy array, sea energy density in ft^2/Hz
    :param swell_energy_density: numpy array, swell energy density in ft^2/Hz
    :param waveFrequency: numpy array, frequency values in Hz
    :param time_index: int, time index for the energy density data to be plotted
    """
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
    """
    Calculate the M0 values for the wave energy density data.

    :param waveEnergyDensity_ft2: numpy array, wave energy density in ft^2/Hz
    :param waveFrequency: numpy array, frequency values in Hz
    :return: numpy array, M0 values
    """
    # Calculate the frequency intervals (delta_f) for the wave energy density data
    delta_f = np.diff(waveFrequency)
    delta_f = np.append(delta_f, delta_f[-1])  # Repeat the last delta_f for the last frequency

    # Integrate the wave energy density spectrum over the frequency range
    m0 = np.sum(waveEnergyDensity_ft2 * delta_f, axis=1)

    return m0


def convert_times_to_timezone(start_time_utc, end_time_utc, timezone):
    # Create timezone objects for UTC and the target timezone
    utc_tz = pytz.timezone('UTC')
    target_tz = pytz.timezone(timezone)

    # Localize the UTC start and end times, then convert them to the target timezone
    end_time = utc_tz.localize(end_time_utc).astimezone(target_tz)
    start_time = utc_tz.localize(start_time_utc).astimezone(target_tz)

    # Return the converted start and end times
    return start_time, end_time

def process_wave_data(waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, waveDirection):
    # Convert wave height (waveHs) from meters to feet
    waveHs_ft = waveHs * 3.28084

    # Convert wave energy density from square meters to square feet
    waveEnergyDensity_ft2 = waveEnergyDensity * 10.764

    # Get the most recent wave direction
    waveDirection_true = waveDirection[0]

    # Return the processed wave data: wave height in feet, wave energy density in square feet, and the most recent wave direction
    return waveHs_ft, waveEnergyDensity_ft2, waveDirection_true

def convert_time_array_timezone(time_array, target_timezone):
    target_tz = pytz.timezone(target_timezone)
    converted_time_array = []
    
    for time in time_array:
        if time.tzinfo is None:
            time = pytz.UTC.localize(time)
        converted_time_array.append(time.astimezone(target_tz).replace(tzinfo=None))
    
    return np.array(converted_time_array, dtype='datetime64')


def main():
    # Set the buoy station ID and time range for data download
    stn = '045'
    end_time_utc = datetime.utcnow()
    start_time_utc = end_time_utc - timedelta(hours=24)

    # Convert start and end times to the Pacific timezone
    start_time, end_time = convert_times_to_timezone(start_time_utc, end_time_utc, 'US/Pacific')

    # Download wave data
    waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, waveDirection, time_array = download_wave_data(stn, start_time, end_time)

    # Process downloaded wave data
    waveHs_ft, waveEnergyDensity_ft2, waveDirection_true = process_wave_data(waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, waveDirection)

    # Calculate the significant wave height
    swh_calculated_ft = calculate_swh(waveEnergyDensity_ft2, waveFrequency)

    # Calculate the M1 ordinate
    m1_ordinate = calculate_m1(waveEnergyDensity_ft2, waveFrequency)

    # Calculate the mean period and length based on the M1 ordinate
    mean_period_m1 = 1/m1_ordinate
    length_m1 = 5.12 * (mean_period_m1 ** 2)

    # Calculate the M1 ordinate for the desired 3' PM SWH
    mean_period_desired = 3.73
    length_desired = 5.12 * (mean_period_desired ** 2)
    m1_ordinate_desired = 1 / mean_period_desired

    # Remove the swell energy from the energy density
    waveEnergyDensity_ft2_no_swell = remove_swell_energy(waveEnergyDensity_ft2, waveFrequency)

    # Calculate the swell energy density
    swell_energy_density = waveEnergyDensity_ft2 - waveEnergyDensity_ft2_no_swell

    # Calculate the modal period of sea energy density for the most recent row
    freq_range_min = 0.025
    freq_range_max = np.min([0.580, np.max(waveFrequency)])
    idx_closest_to_min_freq = np.argmin(np.abs(waveFrequency - freq_range_min))
    idx_closest_to_max_freq = np.argmin(np.abs(waveFrequency - freq_range_max))

    sea_energy_density = waveEnergyDensity_ft2_no_swell[-1, idx_closest_to_min_freq:idx_closest_to_max_freq+1]
    max_energy_density_idx = np.argmax(sea_energy_density)
    modal_period = 1 / waveFrequency[idx_closest_to_min_freq + max_energy_density_idx]

    # Calculate the normalized significant wave height using only the sea energy density
    normalized_swh_ft = calculate_normalized_swh(waveEnergyDensity_ft2, waveFrequency, m1_ordinate, m1_ordinate_desired)

    # Calculate the M0 for the original and modified energy density
    m0_original = calculate_m0(waveEnergyDensity_ft2, waveFrequency)
    m0_modified = calculate_m0(waveEnergyDensity_ft2_no_swell, waveFrequency)

    # Convert time_array to PST
    pst_time_array = convert_time_array_timezone(time_array, 'US/Pacific')

    # Call the function to plot the data
    plot_swh(pst_time_array, waveHs_ft, swh_calculated_ft, normalized_swh_ft)

    # Print the results
    print('Significant Wave Height (ft)', np.round(waveHs_ft[-1], 2))
    print('Estimated Pierson-Moskowitz Wave Height (ft)', np.round(swh_calculated_ft[-1], 2))
    print('Normalized Sea Significant Wave Height (ft)', np.round(normalized_swh_ft[-1], 2))

    # Plot the energy density spectrum
    plot_energy_density_spectrum(waveEnergyDensity_ft2, waveEnergyDensity_ft2_no_swell, swell_energy_density, waveFrequency, -1)

# Call the main function when the script is run
if __name__ == '__main__':
    main()
