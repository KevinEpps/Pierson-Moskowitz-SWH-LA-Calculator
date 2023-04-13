import netCDF4
import numpy as np
from datetime import datetime, timedelta

def fetch_data(stn):
    data_url = f'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/{stn}p1_rt.nc'
    nc = netCDF4.Dataset(data_url)
    return nc

def to_pst(time_var, time_val):
    utc_time = netCDF4.num2date(time_val, time_var.units)
    pst_time = utc_time - timedelta(hours=7)  # Convert to PST
    return pst_time

def main():
    stations = ['045', '264']
    conversion_factor = 3.28084

    for stn in stations:
        nc = fetch_data(stn)

        # Fetch variables
        wave_time_var = nc.variables['waveTime']
        wave_hs_var = nc.variables['waveHs']
        wave_tp_var = nc.variables['waveTp']
        wave_frequency_var = nc.variables['waveFrequency']
        wave_bandwidth_var = nc.variables['waveBandwidth']

        if 'wavePeakPSD' in nc.variables:
            wave_energy_density_var = nc.variables['wavePeakPSD']
        else:
            wave_energy_density_var = nc.variables['waveEnergyDensity']

        # Get the most recent time event
        latest_time_idx = -1
        latest_time_val = wave_time_var[latest_time_idx]
        latest_time_pst = to_pst(wave_time_var, latest_time_val)

        # Get the corresponding values for the latest time event
        latest_wave_hs = wave_hs_var[latest_time_idx] * conversion_factor
        latest_wave_tp = wave_tp_var[latest_time_idx]
        latest_wave_frequency = wave_frequency_var[:]
        latest_wave_bandwidth = wave_bandwidth_var[:]
        if len(wave_energy_density_var.shape) == 2:
          latest_wave_energy_density = wave_energy_density_var[latest_time_idx, :]
        else:
          latest_wave_energy_density = wave_energy_density_var[:]

        # Calculate Pierson-Moskowitz Length Adjusted Significant wave height
        # Comment out the line below when using the most recent set of values
        mean_wave_energy_density = np.mean(latest_wave_energy_density, axis=0)

        # Uncomment the line below when using the most recent set of values
        latest_wave_energy_density_values = latest_wave_energy_density.reshape(-1, 1)
        latest_wave_bandwidth_reshaped = latest_wave_bandwidth.reshape(-1, 1)

        latest_wave_bandwidth_reshaped = latest_wave_bandwidth.reshape(-1)

        # Comment out the line below when using the most recent set of values
        pm_length_adjusted_hs = np.trapz(mean_wave_energy_density * latest_wave_bandwidth, latest_wave_frequency) * 16 * np.pi**2 / (latest_wave_tp**4)

        # Uncomment the line below when using the most recent set of values
        pm_length_adjusted_hs_latest = np.trapz(latest_wave_energy_density_values * latest_wave_bandwidth_reshaped, latest_wave_frequency) * 16 * np.pi**2 / (latest_wave_tp**4)

        # Comment out the line below when using the most recent set of values
        pm_length_adjusted_hs = (pm_length_adjusted_hs * conversion_factor) ** 0.5

        # Uncomment the line below when using the most recent set of values
        pm_length_adjusted_hs_latest = (pm_length_adjusted_hs_latest * conversion_factor) ** 0.5

        print(f"Buoy {stn}:")
        print(f"  Latest Time (PST): {latest_time_pst}")
        print(f"  Significant Wave Height (feet): {latest_wave_hs:.2f}")
        
        # Comment out the line below when using the most recent set of values
        print(f"  Pierson-Moskowitz Length Adjusted Significant Wave Height (feet): {float(pm_length_adjusted_hs):.2f}")

        # Uncomment the line below when using the most recent set of values
        print(f"  Pierson-Moskowitz Length Adjusted Significant Wave Height (latest values, feet): {pm_length_adjusted_hs_latest.compressed()[0]:.2f}")



if __name__ == "__main__":
    main()
