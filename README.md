# CDIP-Buoy-Processor
Python Script to calculate length adjusted significant wave heights using Pierson-Moskowitz Spectrum and centroid offsets.

This program was developed by Kevin Epps from the Engineering Department at the Amphibious Vehicle Test Branch (AVTB).  It's purpose is to download wave buoy data from the Coastal Data Information Program (CDIP) and display it in a way that is useful for test planning and analysis.  CDIP is run by the Scripps Institution of Oceanography (SIO) at the University of California San Diego (UCSD).

https://cdip.ucsd.edu/


# Length-Adjusted Wave Height Calculator
This script calculates and plots the length-adjusted wave height for a CDIP buoy.

# Requirements
To run this script you need:

Python 3.x
NumPy
netCDF4
matplotlib
pytz

# Usage
To run the script, execute the following command:
CDIP PM LA SWH Calculator_v1.0.py

The script currently uses the two main buoys for AVTB in realtime (most recent 12 hour period):

'045':  Oceanside Harbor
'264':  Red Beach

Future versions may prompt the user for the buoy number and the time window to download data for.

# Output
The script will output the following information:

Measured significant wave height, 
Pierson-Moskowitz significant wave height, 
Centroid of the measured wave spectrum, 
Centroid of Pierson-Moskowitz spectrum, 
Effective significant wave height, 
It will also display a plot of the length-adjusted wave height over time.

# Contributing
Contributions to this script are welcome! If you find a bug or have an idea for a new feature, please open an issue or submit a pull request.
