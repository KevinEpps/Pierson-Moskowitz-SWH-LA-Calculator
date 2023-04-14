import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz

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
    time_array = netCDF4.num2date(time_var[start_index:end_index+1], time_var.units).tolist()

    nc.close()

    return waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, time_array

# Calculate the effective significant wave height for each time step
def calculate_effective_wave_height(waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz):
    g = 9.81
    Hsea_list = []
    HsPM_list = []
    Cmeas_list = []
    CPM_list = []
    
    # Loop through each time step in the input data and calculate Hsea, HsPM, Cmeas, and CPM for that time step
    for i in range(len(waveHs)):
        waveHs_i = waveHs[i]
        waveTp_i = waveTp[i]
        waveEnergyDensity_i = waveEnergyDensity[i, :]
        waveTz_i = waveTz[i]

        # Calculate the centroid frequency of the measured spectrum (Cmeas)
        Cmeas = np.sum(waveFrequency * waveEnergyDensity_i) / np.sum(waveEnergyDensity_i)
        Cmeas_list.append(Cmeas)

        # Calculate the Pierson-Moskowitz significant wave height (HsPM) and use it to calculate B and fp for the PM spectrum
        HsPM = 0.21 * (g * waveTp_i / (2 * np.pi)) ** (1/2)
        HsPM_list.append(HsPM)
        B = (5/4) * (HsPM ** 2) / (g ** 2) #Testing to see if it works
        fp = 1.0 / (1.56 * waveTz_i)

        # Calculate the centroid frequency of the Pierson-Moskowitz spectrum (CPM)
        CPM = g / (2 * np.pi) * np.sqrt(B / fp) #Testing to see if it works
        #CPM = 3.11 * HsPM / fp
        CPM_list.append(CPM)
        
        # Calculate the effective significant wave height for the current time step
        Hsea = waveHs_i * Cmeas / CPM
        Hsea_list.append(Hsea)

    return np.array(Hsea_list), np.array(HsPM_list), np.array(Cmeas_list), np.array(CPM_list)

def plot_wave_height(time_array, Hsea_feet, stn):

    # Create a new figure and axes objects
    fig, ax = plt.subplots()

    # Plot the Hsea_feet values against the corresponding time_array values
    ax.plot(time_array, Hsea_feet, label='Length-Adjusted Wave Height')

    # Set the title and axis labels for the plot
    plt.title(f'Length-Adjusted Wave Height for Buoy {stn}')
    plt.xlabel('Time (PST)')
    plt.ylabel('Wave Height (ft)')

    # Add a grid to the plot
    plt.grid()

    # Set the date format for the x-axis labels
    date_format = '%m/%d/%Y %H:%M'
    date_formatter = plt.matplotlib.dates.DateFormatter(date_format)
    ax.xaxis.set_major_formatter(date_formatter)

    # Rotate the x-axis labels for better readability
    fig.autofmt_xdate()

    # Add a legend to the plot and display it
    plt.legend()
    plt.show()

#Main function
def main():

    # Set the buoy station ID and time range for data download
    stn = '045'
    end_time_utc = datetime.utcnow()
    start_time_utc = end_time_utc - timedelta(hours=12)

    # Set the timezones for the start and end times
    utc_tz = pytz.timezone('UTC')
    pst_tz = pytz.timezone('US/Pacific')

    # Convert the start and end times to the Pacific timezone
    end_time = utc_tz.localize(end_time_utc).astimezone(pst_tz)
    start_time = utc_tz.localize(start_time_utc).astimezone(pst_tz)

    # Download the wave data for the specified time range
    waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz, time_array = download_wave_data(stn, start_time, end_time)

    # Calculate the effective wave height and Pierson-Moskowitz significant wave height
    Hsea, HsPM, Cmeas, CPM = calculate_effective_wave_height(waveHs, waveTp, waveFrequency, waveEnergyDensity, waveTz)

    # Convert the wave height units from meters to feet
    feet_per_meter = 3.28084
    waveHs_feet = waveHs * feet_per_meter
    Cmeas_feet = Cmeas * feet_per_meter
    CPM_feet = CPM * feet_per_meter
    Hsea_feet = Hsea * feet_per_meter
    HsPM_feet = HsPM * feet_per_meter

    # Print out the measured, Pierson-Moskowitz, and effective wave heights
    print(f"Measured significant wave height: {waveHs_feet} ft")
    print(f"Pierson-Moskowitz significant wave height: {HsPM_feet} ft")
    print(f"Centroid of the measured wave spectrum: {Cmeas_feet} ft")
    print(f"Centroid of Pierson-Moskowitz spectrum: {CPM_feet} ft")
    print(f"Effective significant wave height: {Hsea_feet} ft")

    #time_array = netCDF4.num2date(time_var[start_index:end_index+1], time_var.units).tolist() #not needed when called and stored in def download
    time_array = [datetime.fromisoformat(str(t)).astimezone(pst_tz) for t in time_array]
    time_array = [t - timedelta(hours=6) for t in time_array]

    plot_wave_height(time_array, Hsea_feet, stn)

if __name__ == "__main__":
    main()
