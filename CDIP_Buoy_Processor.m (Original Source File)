%% CDIP Buoy Processor

% This program was developed by Mike Slivka from the Instrumentation and
% Metrology Laboratory at the Amphibious Vehicle Test Branch (AVTB).  It's
% purpose is to download wave buoy data from the Coastal Data Information
% Program (CDIP) and display it in a way that is useful for test planning
% and analysis.  CDIP is run by the Scripps Institution of Oceanography
% (SIO) at the University of California San Diego (UCSD).
% 
% https://cdip.ucsd.edu/

    if isdeployed % Stand-alone mode.
        [status, result] = system('path');
        currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
    else % MATLAB mode.
        currentDir = pwd;
        addpath(genpath('./')); %Add the code library to the path
        addpath(genpath('.\')); %Add the code library to the path
    end


    
%% Get User Input
% This will open the dialog to get information from the user on what they
% want.

    dialog_cdip % Get input from the user
    
    programstart = clock; % This will be used at the end to caluclate the total process time



%% Make sure previous figures aren't open to interfere
% This is only here because there were a few times this was run with a figure
% already opened and it caused some interference issues...

    close all

    
    
%% Variables

    use_local_time = false; % True uses the "utc_offset", false uses UTC time for plot. 

    Magnetic_Declination = 11.53199; % From NOAA for zip code 92054 (Camp Pendleton) as of 2019-Aug-14
    
    upload_via_ftp = false; % Future use 
       
    
    
%%  Verify Output Folder
% This verifies the selected output folder is valid if you have selected to
% plot the data.  If it is not valid, it will prompt for you to select one
% rather than just error out.

    if create_plots == true % No need to save if you aren't plotting anything
        if 7 == exist(output_folder,'dir') % It exists, but is it a valid location?
            % There is nothing to do if it already exists...
        else
            disp('The workspace variable "output_folder" is not valid.')
            clear output_folder % Clears the invalid "output_folder".
            disp('Select the Output Folder location...')
            output_folder = uigetdir('Choose Output Folder');      
        end
    end
    


%% CONNECT TO THREDDS SERVER AND OPEN NETCDF FILE
% This address is created using variables at the top of this script.
% 43 is Pendleton, 43f is the forecast for Pendleton
% 45 is Oceanside, 45f is the forecast for Oceanside
% 71 is Harvest, 71f is the forecast for Harvest
% D1064 is the Transect buoy for AVTB Beach
% B1418 is the Transect buoy for the test site at Wall Beach at VAFB

    tic
    disp(' ') % Adds a space in the command window 
    disp('Making a connection to the data on the CDIP Thredds server...')   
    
    address1 = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/',buoy_number,'p1_rt.nc'];
    address2 = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_validation/BP',buoy_number,'_forecast.nc'];
    address3 = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/',buoy_number,'p1/',buoy_number,'p1_historic.nc'];
    
    toc

%     Transect buoys for future use...
%     dsD1064 = ncdataset('http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/D1064_nowcast.nc');
%     dsB1418 = ncdataset('http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/B1418_nowcast.nc');

    

%% CONVERT START/END DATES TO MATLAB SERIAL NUMBERS

    startser = datenum(startdate);
    endser = datenum(enddate);

    

%% CALL 'TIME' VARIABLE 

    disp(' ') % Adds a space in the command window
    disp('Comparing times of available data versus what has been requested and creating "Time" variables...')
    tic

    time_available_ser = double(ncread(address1,'waveTime'));
    time_available_ser = unixtime2mat(time_available_ser); % Converts from UNIX time to MATLAB serial time
    time_available = datetime(time_available_ser,'ConvertFrom','datenum'); % Converted to date/time
	time_available_Table = table(time_available,'VariableNames',{'utc_time'});

    if plot_forecast == true
        timevar_f = double(ncread(address2,'waveTime'));
        timeconvert_f = unixtime2mat(timevar_f);
        timedt_f = datetime(timeconvert_f,'ConvertFrom','datenum'); % Convert Matlab serial units 'timecut' section to datetime objects
        timedt_f_Table = table(timedt_f,'VariableNames',{'utc_time'});
    end

    if plot_archive == true
        timevar_a = double(ncread(address3,'waveTime'));
        timeconvert_a = unixtime2mat(timevar_a);
        timedt_a = datetime(timeconvert_a,'ConvertFrom','datenum'); % Convert Matlab serial units 'timecut' section to datetime objects
        timedt_a_Table = table(timedt_a,'VariableNames',{'utc_time'});
    end
    
    toc

    
    
%% ASSIGN NAMES TO VARIABLES TO BE USED IN COMPENDIUM PLOT, AND SPECIFY TIMERANGE
% This pulls out specific variables we want to work with and creates arrays
% for each one.  The "f" is for any forcast data from that buoy as opposed to realtime.

    disp(' ') % Adds a space in the command window
    disp('Downloading data and assigning variables. This can take a moment or two...')
    tic
    
    f = waitbar(0,'Downloading data and assigning variables...');
   
    waitbar(0,f,'Loading Significant Wave Height');
    Significant_Wave_Height_Meters = single(ncread(address1,'waveHs'));
    Significant_Wave_Height_Meters_Table = table(Significant_Wave_Height_Meters,'VariableNames',{'Significant_Wave_Height_Meters'});
    
    waitbar(.10,f,'Loading Significant Wave Height');
    Significant_Wave_Height_Feet = round((Significant_Wave_Height_Meters*3.28084),2);
    Significant_Wave_Height_Feet_Table = table(Significant_Wave_Height_Feet,'VariableNames',{'Significant_Wave_Height_Feet'});
    
    waitbar(.20,f,'Loading Peak Wave Period');
    Peak_Wave_Period_Seconds = single(ncread(address1,'waveTp'));
    Peak_Wave_Period_Seconds_Table = table(Peak_Wave_Period_Seconds,'VariableNames',{'Peak_Wave_Period_Seconds'});
    
    waitbar(.30,f,'Loading Peak Wave Direction');
    Peak_Wave_Direction_Deg_True = single(ncread(address1,'waveDp'));
    Peak_Wave_Direction_Deg_True_Table = table(Peak_Wave_Direction_Deg_True,'VariableNames',{'Peak_Wave_Direction_Deg_True'});
    
    waitbar(.40,f,'Loading Average Wave Period');
    Average_Wave_Period_Seconds = single(ncread(address1,'waveTa'));
    Average_Wave_Period_Seconds_Table = table(Average_Wave_Period_Seconds,'VariableNames',{'Average_Wave_Period_Seconds'});
    
    waitbar(.50,f,'Loading Sea Surface Temperature');
    Sea_Surface_Temp_DegC = single(ncread(address1,'sstSeaSurfaceTemperature'));
    Sea_Surface_Temp_DegC_Table = table(Sea_Surface_Temp_DegC,'VariableNames',{'Sea_Surface_Temp_DegC'});
    
    waitbar(.55,f,'Loading Sea Surface Temperature');
    Sea_Surface_Temp_DegF = round((Sea_Surface_Temp_DegC*1.8+32),1);
    Sea_Surface_Temp_DegF_Table = table(Sea_Surface_Temp_DegF,'VariableNames',{'Sea_Surface_Temp_DegF'});
    
    waitbar(.60,f,'Loading Wave Frequency');
    Wave_Frequency = ncread(address1,'waveFrequency');
    
    waitbar(.70,f,'Loading Wave Energy');
    Wave_Energy_Raw = single(ncread(address1,'waveEnergyDensity'));
    Wave_Energy_Raw = Wave_Energy_Raw';
    Wave_Energy_Raw_Table = table(Wave_Energy_Raw,'VariableNames',{'Wave_Energy_Raw'});
    
    waitbar(.80,f,'Loading Wave Mean Direction');
    Wave_Mean_Direction_Deg_True = single(ncread(address1,'waveMeanDirection'));
    Wave_Mean_Direction_Deg_True = Wave_Mean_Direction_Deg_True';
    Wave_Mean_Direction_Deg_True_Table = table(Wave_Mean_Direction_Deg_True,'VariableNames',{'Wave_Mean_Direction_Deg_True'});

    waitbar(.90,f,'Loading Wave Bandwidth');    
    Wave_Bandwidth = single(ncread(address1,'waveBandwidth'));    

    len_WaveFreq = length(Wave_Frequency);    
    
    waitbar(1.0,f,'Finished Loading Realtime Data');   
    close(f);% Close waitbar
    
    
    if plot_forecast == true

        f = waitbar(0,'Downloading forecast data and assigning variables...');

        waitbar(0,f,'Loading Forecast Significant Wave Height');
        Significant_Wave_Height_Meters_Forecast = single(ncread(address2,'waveHs')); % The "f" implies forcast data as oppsed to realtime.
        Significant_Wave_Height_Meters_Forecast_Table = table(Significant_Wave_Height_Meters_Forecast,'VariableNames',{'Significant_Wave_Height_Meters'});

        Significant_Wave_Height_Feet_Forecast = round((Significant_Wave_Height_Meters_Forecast*3.28084),2);
        Significant_Wave_Height_Feet_Forecast_Table = table(Significant_Wave_Height_Feet_Forecast,'VariableNames',{'Significant_Wave_Height_Feet'});        

        waitbar(.16,f,'Loading Forecast Peak Wave Period');
        Peak_Wave_Period_Seconds_Forecast = single(ncread(address2,'waveTp')); % The "f" implies forcast data as oppsed to realtime.
        Peak_Wave_Period_Seconds_Forecast_Table = table(Peak_Wave_Period_Seconds_Forecast,'VariableNames',{'Peak_Wave_Period_Seconds'});

        waitbar(.32,f,'Loading Forecast Peak Wave Direction');
        Peak_Wave_Direction_Forecast = single(ncread(address2,'waveDp')); % The "f" implies forcast data as oppsed to realtime.
        Peak_Wave_Direction_Forecast_Table = table(Peak_Wave_Direction_Forecast,'VariableNames',{'Peak_Wave_Direction'});

        waitbar(.48,f,'Loading Forecast Average Wave Period');
        Average_Wave_Period_Seconds_Forecast = single(ncread(address2,'waveTa')); % The "f" implies forcast data as oppsed to realtime.
        Average_Wave_Period_Seconds_Forecast_Table = table(Average_Wave_Period_Seconds_Forecast,'VariableNames',{'Average_Wave_Period_Seconds'});

        waitbar(.64,f,'Loading Forecast Wave Energy');      
        Wave_Energy_Raw_Forecast = single(ncread(address2,'waveEnergyDensity'));
        Wave_Energy_Raw_Forecast = Wave_Energy_Raw_Forecast'; % Rotate 90 degrees
        Wave_Energy_Raw_Forecast_Table = table(Wave_Energy_Raw_Forecast,'VariableNames',{'Wave_Energy_Raw'});

        waitbar(.80,f,'Loading Forecast Wave Mean Direction');
        Wave_Mean_Direction_Deg_True_Forecast = single(ncread(address2,'waveMeanDirection'));
        Wave_Mean_Direction_Deg_True_Forecast = Wave_Mean_Direction_Deg_True_Forecast'; % Rotate 90 degrees
        Wave_Mean_Direction_Deg_True_Forecast_Table = table(Wave_Mean_Direction_Deg_True_Forecast,'VariableNames',{'Wave_Mean_Direction_Deg_True'});

        waitbar(1.0,f,'Finished Loading Forecast Data');   
        close(f);% Close waitbar
        
    end
    
    if plot_archive == true

        f = waitbar(0,'Downloading forecast data and assigning variables...');

        waitbar(0.0,f,'Loading Archive Significant Wave Height');
        Significant_Wave_Height_Meters_Archive = single(ncread(address3,'waveHs')); 
        Significant_Wave_Height_Meters_Archive_Table = table(Significant_Wave_Height_Meters_Archive,'VariableNames',{'Significant_Wave_Height_Meters'});

        waitbar(0.11,f,'Loading Archive Significant Wave Height');
        Significant_Wave_Height_Feet_Archive = round((Significant_Wave_Height_Meters_Archive*3.28084),2);
        Significant_Wave_Height_Feet_Archive_Table = table(Significant_Wave_Height_Feet_Archive,'VariableNames',{'Significant_Wave_Height_Feet'});
        
        waitbar(.22,f,'Loading Archive Peak Wave Period');
        Peak_Wave_Period_Seconds_Archive = single(ncread(address3,'waveTp')); 
        Peak_Wave_Period_Seconds_Archive_Table = table(Peak_Wave_Period_Seconds_Archive,'VariableNames',{'Peak_Wave_Period_Seconds'});
        
        waitbar(.33,f,'Loading Archive Peak Wave Direction');
        Peak_Wave_Direction_Archive = single(ncread(address3,'waveDp')); 
        Peak_Wave_Direction_Archive_Table = table(Peak_Wave_Direction_Archive,'VariableNames',{'Peak_Wave_Direction'});
        
        waitbar(.44,f,'Loading Archive Average Wave Period');
        Average_Wave_Period_Seconds_Archive = single(ncread(address3,'waveTa'));
        Average_Wave_Period_Seconds_Archive_Table = table(Average_Wave_Period_Seconds_Archive,'VariableNames',{'Average_Wave_Period_Seconds'});

        waitbar(.55,f,'Loading Archive Wave Energy');
        Wave_Energy_Raw_Archive = single(ncread(address3,'waveEnergyDensity'));
        Wave_Energy_Raw_Archive = Wave_Energy_Raw_Archive'; % Rotate 90 degrees
        Wave_Energy_Raw_Archive_Table = table(Wave_Energy_Raw_Archive,'VariableNames',{'Wave_Energy_Raw'});
        
        waitbar(.66,f,'Loading Archive Wave Mean Direction');
        Wave_Mean_Direction_Deg_True_Archive = single(ncread(address3,'waveMeanDirection'));
        Wave_Mean_Direction_Deg_True_Archive = Wave_Mean_Direction_Deg_True_Archive'; % Rotate 90 degrees
        Wave_Mean_Direction_Deg_True_Archive_Table = table(Wave_Mean_Direction_Deg_True_Archive,'VariableNames',{'Wave_Mean_Direction_Deg_True'});
        
        waitbar(.77,f,'Loading Archive Sea Surface Temperature');
        Sea_Surface_Temp_DegC_Archive = single(ncread(address3,'sstSeaSurfaceTemperature'));
        Sea_Surface_Temp_DegC_Archive_Table = table(Sea_Surface_Temp_DegC_Archive,'VariableNames',{'Sea_Surface_Temp_DegC'});
    
        waitbar(.88,f,'Loading Archive Sea Surface Temperature');
        Sea_Surface_Temp_DegF_Archive = round((Sea_Surface_Temp_DegC_Archive*1.8+32),1);
        Sea_Surface_Temp_DegF_Archive_Table = table(Sea_Surface_Temp_DegF_Archive,'VariableNames',{'Sea_Surface_Temp_DegF'});
        
        waitbar(1.0,f,'Finished Loading Archive Data');   
        close(f);% Close waitbar
        
    end    
    
    toc



%% Create Master Tables

    buoy_data = [time_available_Table Average_Wave_Period_Seconds_Table Peak_Wave_Direction_Deg_True_Table ...
        Peak_Wave_Period_Seconds_Table Significant_Wave_Height_Meters_Table ...
        Significant_Wave_Height_Feet_Table Wave_Energy_Raw_Table Wave_Mean_Direction_Deg_True_Table];
    
    buoy_data_SST = [time_available_Table Sea_Surface_Temp_DegC_Table Sea_Surface_Temp_DegF_Table];

    if plot_forecast == true
        
        buoy_data_forecast = [timedt_f_Table Significant_Wave_Height_Meters_Forecast_Table ...
            Significant_Wave_Height_Feet_Forecast_Table Peak_Wave_Period_Seconds_Forecast_Table ... 
            Peak_Wave_Direction_Forecast_Table Average_Wave_Period_Seconds_Forecast_Table ...
            Wave_Mean_Direction_Deg_True_Forecast_Table Wave_Energy_Raw_Forecast_Table];
    
    end

    if plot_archive == true

        buoy_data_archive = [timedt_a_Table Significant_Wave_Height_Meters_Archive_Table ...
            Significant_Wave_Height_Feet_Archive_Table Peak_Wave_Period_Seconds_Archive_Table ... 
            Peak_Wave_Direction_Archive_Table Average_Wave_Period_Seconds_Archive_Table ...
            Wave_Mean_Direction_Deg_True_Archive_Table Wave_Energy_Raw_Archive_Table];
       
    end
    
    
    
%% FIND INDEX NUMBERS OF CLOSEST TIME ARRAY VALUES TO START/END DATES

    diffstart = abs(time_available_ser-startser); % Compute the difference between the 'startdate' serial number and every serial number in 'timevar', to determine which difference is smallest (which number most closely matches 'startser')
    [~,startidx] = min(diffstart); % index of closest object to 'startser', based on the minimum value in 'diffstart' differences array
    % closeststart = timeconvert(startidx);  % value of closest object to 'startser' (optional - not used in code)

    diffend = abs(time_available_ser-endser);
    [~,endidx] = min(diffend); % index of closest object to 'endser', based on the minimum value in 'diffstart' differences array
    % closestend = timeconvert(endidx); % value of closest object to 'endser' (optional - not used in code)
   
    % Create Subset of data
    buoy_data_subset = buoy_data(startidx:endidx,1:end);
    buoy_data_SST_subset = buoy_data_SST(startidx:endidx,1:end);

    

%% Create Table if it doesn't exist

    length_adjusted_table = table();
    Sea_State_table = table();
    Estimated_Surf_Height_table = table();
    Estimated_Surf_Period_table = table();

        
    
%% Start Calculations
    
    disp(' ')
    disp('Running the data through analysis equations...')
    tic
        
    for m = 1:height(buoy_data_subset)

        process_number = m;
        
    %% Pull Out Single Values
    % Some of the variables pulled from the CDIP realtime data is just the most
    % recent study, and some will be all 7 days worth (or however much data you
    % requested).  This section pulls the most recent value for specific
    % variables.  Any variable with the "_Now" suffix is the datapoint that
    % represents the most recent study as opposed to an array with 7 days of data.


        Significant_Wave_Height_Feet_Now = buoy_data_subset.Significant_Wave_Height_Feet(process_number);  
        Significant_Wave_Height_Feet_Now_str = num2str(round(Significant_Wave_Height_Feet_Now,2)); % Create string for later use

        Peak_Wave_Period_Seconds_Now = buoy_data_subset.Peak_Wave_Period_Seconds(process_number);
        Peak_Wave_Period_Seconds_Now_str = num2str(round(Peak_Wave_Period_Seconds_Now,1)); % Create string for later use

        Peak_Wave_Direction_Deg_True_Now = buoy_data_subset.Peak_Wave_Direction_Deg_True(process_number);
        Peak_Wave_Direction_Deg_True_Now_str = num2str(round(Peak_Wave_Direction_Deg_True_Now,1)); % Create string for later use

        Average_Wave_Period_Seconds_Now = buoy_data_subset.Average_Wave_Period_Seconds(process_number);

        Sea_Surface_Temp_DegF_Now = buoy_data_SST_subset.Sea_Surface_Temp_DegF(process_number);    
        Sea_Surface_Temp_DegF_Now_str = num2str(round(Sea_Surface_Temp_DegF_Now,1));

        Wave_Time_UTC_Now = buoy_data_subset.utc_time(process_number);
        Wave_Time_Local_Now = Wave_Time_UTC_Now - hours(utc_offset);
        Wave_Time_Local_Now_str = datestr(Wave_Time_Local_Now);% Converts date to string for later use



    %% Wave Energy
    % This part is fuzzy.  When you compare the raw SPT file from the buoy,
    % then look at the "waveEnergyDensity" variable from the CDIP server,
    % this is what we did to make it look like what the Length PM formula
    % is expecting.  This makes the PM come out to what we think it's supposed
    % to be.

    % The exact description for that value from CDIP is:
    % "waveEnergyDensity - Band Energy Density (m^2 second).  2D variable
    % based on (waveTime, waveFrequency).  Min: 0."

    % The bottom row of the 2D array is the most recent values.  This cuts
    % that from the 2D array and puts it in it's own 1D array, rotated 90 degrees
    % so it can be used appropriately.

        Wave_Energy_Cut = Wave_Energy_Raw; % Keeps the original Wave Energy as "Raw" for comparison purposes

        Wave_Energy = buoy_data_subset.Wave_Energy_Raw(process_number,1:end);

        Wave_Energy = Wave_Energy'; % Rotate the array 90 degrees, it needs to be 64x1, not 1x64.



    %% All the equations below were taken from the DatawellV08-yyyy-mm-dd-hhmm.xls file.  
    % These equations were built into the excel template and used for data
    % analysis. Some of the equations are still being looked at for more
    % clarification as to why they are performed that way.



    % Create an array for the Bandwidth
        bandwidth = zeros(len_WaveFreq,1);

    % Calculate a bandwidth using the formula below.  The formula uses the
    % midpoint formula to find the bandwidth
        for x = 2:len_WaveFreq-1
            bandwidth(x,1) = (Wave_Frequency(x,1)-Wave_Frequency((x-1),1))/2+(Wave_Frequency((x+1),1)-Wave_Frequency(x,1))/2;
        end

    % The first bandwidth point is equal to the second bandwidth point
        bandwidth(1,1) = bandwidth(2,1);

    % The last bandwidth point is equal the second to last bandwidth point.
        bandwidth(len_WaveFreq,1) = bandwidth(len_WaveFreq-1,1);



    %% Calculate the Smoothed direction from the mean direction. 
    % This is done by a linear regression.

        % This pulls only the most recent line of wave mean direction data from
        % the 7 day download.  It then rotates the array for further use.
        % These two lines were added to make the data from the CDIP website
        % look like the data the Excel formula was expecting to see.
    
        Wave_Mean_Direction_Deg_True = buoy_data_subset.Wave_Mean_Direction_Deg_True(process_number,1:end);
        Wave_Mean_Direction_Deg_True = Wave_Mean_Direction_Deg_True';

        Smoothed_Direction = zeros(len_WaveFreq,1);
        Smoothed_Direction(2,1) = Wave_Mean_Direction_Deg_True(1,1)/4+Wave_Mean_Direction_Deg_True(2,1)/2+Wave_Mean_Direction_Deg_True(3,1)/4;
            for x = 3:len_WaveFreq-2
                Smoothed_Direction(x,1) = Wave_Mean_Direction_Deg_True(x-2,1)/9+2*Wave_Mean_Direction_Deg_True(x-1,1)/9+3*Wave_Mean_Direction_Deg_True(x,1)/9+2*Wave_Mean_Direction_Deg_True(x+1,1)/9+Wave_Mean_Direction_Deg_True(x+2,1)/9;
            end
        Smoothed_Direction(len_WaveFreq-1,1) = Smoothed_Direction(len_WaveFreq-2,1);
        Smoothed_Direction(len_WaveFreq,1) = Smoothed_Direction(len_WaveFreq-2,1);

        % For the first data point find the smoothed direction point by using
        % interpolation
        Smoothed_Direction(1,1) =(Smoothed_Direction(3,1)-Smoothed_Direction(2,1))*(Wave_Frequency(1,1)-Wave_Frequency(2,1))/(Wave_Frequency(3,1)-Wave_Frequency(2,1))+Smoothed_Direction(2,1);



    %% Wave Period
    % Calculate the Wave Period. The wave period is equal to 1/f

        Wave_Period = Wave_Frequency.^(-1);



    %% Wave Length    
    % Wave length = 5.12*Wave Period^2

        Wave_Length = 5.12*Wave_Period.^2;



    %% Total Energy Density    
    % This is done by multiplying the raw energy data by a factor of 10.76

        Total_Energy_Density = Wave_Energy*10.76391042;

    % Calculate the Smoothed Total Energy Density.  Use the formula below to
    % calculate the Smoothed Total Energy Density from the Total Energy Density

        Smoothed_Total_Energy_Density = zeros(len_WaveFreq,1);
            for x = 2:len_WaveFreq-1
                Smoothed_Total_Energy_Density(x,1) = Total_Energy_Density(x-1,1)/4+Total_Energy_Density(x,1)/2+Total_Energy_Density(x+1,1)/4;
            end
        Smoothed_Total_Energy_Density(len_WaveFreq,1) = Smoothed_Total_Energy_Density(len_WaveFreq-1,1);


    % Calculate the M0 Ordinate using the smooth total energy density. This is
    % done by multiplying the bandwidth and the smoothed total energy density
    % together

        Smoothed_Total_Energy_Density_M0 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Smoothed_Total_Energy_Density_M0(x,1) = bandwidth(x,1)*Smoothed_Total_Energy_Density(x,1);
            end


    % Calculate the Smoothed Total Energy Density M1 Ordinate by multiplying the Smoothed 
    % M0 Ordinate by the frequency

        Smoothed_Total_Energy_M1 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Smoothed_Total_Energy_M1(x,1) = Wave_Frequency(x,1)*Smoothed_Total_Energy_Density_M0(x,1);
            end


    % Calculate the smoothed slope energy density which is the difference
    % between the smoothed total energy density 

        Smoothed_Slope_Energy_Density = zeros(len_WaveFreq,1);
            for x = 2:len_WaveFreq
                Smoothed_Slope_Energy_Density(x,1) = Smoothed_Total_Energy_Density(x,1)-Smoothed_Total_Energy_Density(x-1,1);
            end


    %Create the swell seperation Marker
    % Create a zero matrix for the swell seperation marker

        Swell_Seperation_Marker = zeros(len_WaveFreq,1);


    % Create another zero matrix to find the minimum point in the Smoothed
    % Total Energy Density

        Smoothed_Total_Energy_Density_min_mat = zeros(6,1);

        y = 1; % y represents the row number in the smoothed total energy density and will be used in the for loop below


    % Find the minimum smoothed total energy density from the 11th - 16th
    % points.  This was how this was performed in the datawell buoy document
    % and it is still being investigated as to why.

        for x = 11:16
            Smoothed_Total_Energy_Density_min_mat(y,1) = Smoothed_Total_Energy_Density(x,1);
            y = y+1;
        end

        Smoothed_Total_Energy_Density_min = min(Smoothed_Total_Energy_Density_min_mat);


    % The following for loop is to be used to find the swell seperation marker
    % at the minimum value for the smoothed total energy density.

        for x = 1:len_WaveFreq
            if Smoothed_Total_Energy_Density_min == Smoothed_Total_Energy_Density(x,1)
                Swell_Seperation_Marker(x,1) = 2;
                f_sted = Wave_Frequency(x,1);
                Swell_Seperation_Marker_2 = x;
            end
        end


    % Calculate the Swell Seperation Index

        Swell_Seperation_Index = zeros(len_WaveFreq,1);


    % The Swell seperation index is calculated by adding the swell seperation
    % marker to the previous swell seperation index value multiplied by a
    % factor of 1.01

        for x =2:len_WaveFreq
            Swell_Seperation_Index(x,1) = Swell_Seperation_Marker(x,1)+Swell_Seperation_Index(x-1,1)*1.01;
        end



    %% Calculate the Sea Energy Density

        Sea_Energy_Density = zeros(len_WaveFreq,1);
        Swell_Seperation_Marker_sum = Swell_Seperation_Marker(x,1);


    % The Sea Energy density equation is calculated using the frequency values
    % and the smoothed total energy density before the swell seperation marker
    % hits 2.0. 

        for x = 1:Swell_Seperation_Marker_2
            Sea_Energy_Density(x,1) = Smoothed_Total_Energy_Density_min*((Wave_Frequency(x,1)-Wave_Frequency(1,1))/(f_sted-Wave_Frequency(1,1)))^8;
        end


    % After the swell seperation marker hits 2.0, the sea energy density is
    % equl to the smoothed total energy density

        for x = (Swell_Seperation_Marker_2+1):len_WaveFreq
            Sea_Energy_Density(x,1) = Smoothed_Total_Energy_Density(x,1);
        end


    % Calculate the M0 Ordinate using the Sea Energy Density and the bandwidth

        Sea_Energy_Density_M0 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Sea_Energy_Density_M0(x,1) = bandwidth(x,1)*Sea_Energy_Density(x,1);
            end


    % Calculate the M1 Ordinate using the M0 ordinate and the frequency

        Sea_Energy_Density_M1 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Sea_Energy_Density_M1(x,1) = Wave_Frequency(x,1)*Sea_Energy_Density_M0(x,1);
            end


    % Calculate the Swell Energy Density by taking the difference of the sea
    % energy density and the smoothed total energy density

        Swell_Energy_Density = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                if   Smoothed_Total_Energy_Density(x,1)-Sea_Energy_Density(x,1)<0
                     Swell_Energy_Density(x,1) =0;
                else
                    Swell_Energy_Density(x,1) = Smoothed_Total_Energy_Density(x,1)-Sea_Energy_Density(x,1);
                end
            end


    % Calculate the M0 Ordinate using the Swell Energy Density and the
    % bandwidth

        Swell_Energy_Density_M0 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Swell_Energy_Density_M0(x,1) = bandwidth(x,1)*Swell_Energy_Density(x,1);
            end


    % Calculate the M1 Ordinate using the swell energy density M0 ordinate and
    % the frequency

        Swell_Energy_Density_M1 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                Swell_Energy_Density_M1(x,1) = Wave_Frequency(x,1)*Swell_Energy_Density_M0(x,1);
            end



    %% Calculate the Bretschneider & Pierson Moskowitz Energy Densities

    % Calculate the frequency in rad/sec

        f_rad = Wave_Frequency*2*pi;


    % Convert bandwidth to rad/sec

        band_rad = bandwidth*2*pi;


    % Calculate the H1/3 value, this value was not defined in the excel sheet

        H = 4*sum(Swell_Energy_Density_M0).^(.5);


    % Calculate the t-modal (wave period @ max swell_ED)
    % Find the max for the Swell Energy Density

        Swell_max = max(Swell_Energy_Density);
            for x = 1:len_WaveFreq
                if Swell_Energy_Density(x,1)==Swell_max
                   T_modal_BS = Wave_Period(x,1);
                end
            end


    % Calculate the W-modal for BretSchneider using the wave period

        W_modal_BS = 2*pi/T_modal_BS;


    % Convert W-modal from rps to Hz
    % W_modBS = W_modBS/(2*pi);


    % Calculate the coefficients for the Bretschneider equation

        a_coef_BS = W_modal_BS^4*H^2*1.25/4;
        b_coef_BS = -1.25*W_modal_BS^4;


    % Calculate the BrettSchneider Energy Density using the BrettSchneider
    % Equation

        BS = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                BS(x,1) = a_coef_BS/(f_rad(x,1)+.0001)^5*exp(b_coef_BS/(f_rad(x,1)+.0001)^4);
            end


    % Calculate the M0 Ordinate using the Bretschneider Energy Density

        BS_M0 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                BS_M0(x,1) = bandwidth(x,1)*BS(x,1);
            end


    % Calculate the M1 Ordinate using the smooth total energy density

        BS_M1 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                BS_M1(x,1) = Wave_Frequency(x,1)*BS_M0(x,1);
            end


    % Calculate the H1/3 value.  The equation uses the sum of all the Sea
    % Energy Density M0 Ordinate values 

        H_PM = 4*sum(Sea_Energy_Density_M0).^(.5);


    % Calculate the coefficients for the Pierson_Moskowitz equation

        a_coef_PM = .0081*32.2^2; %32.2 is gravity in ft/s^2
        b_coef_PM = -.032*(32.2/H_PM)^2; %32.2 is gravity in ft/s^2


    % Calculate the Pierson_Moskowitz Energy Density using the Pierson
    % Moskowitz equation

        PM = zeros(len_WaveFreq,1);
        for x = 1:len_WaveFreq
             PM(x,1) = a_coef_PM/(f_rad(x,1)+.0001)^5*exp(b_coef_PM/(f_rad(x,1)+.0001)^4);
        end


    % Calculate the M0 Ordinate using the Pierson-Moskowitz Energy Density

        PM_M0 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                PM_M0(x,1) = bandwidth(x,1)*PM(x,1);
            end


    % Calculate the M1 Ordinate using the Pierson-Moskowitz energy density

        PM_M1 = zeros(len_WaveFreq,1);
            for x = 1:len_WaveFreq
                PM_M1(x,1) = Wave_Frequency(x,1)*PM_M0(x,1);
            end


    % Find the Target 3' Sig Pierson-Moskowitz

        a_coef = .0081*32.2^2; % a_coef equation
        b_coef = -.032*(32.2/3)^2; % b_coef equation
        sigPM = zeros(len_WaveFreq,1);


    % Calculate the Target 3' Sig Pierson Moskowitz using the 3' Sig equation

        for x = 1:len_WaveFreq
            sigPM(x,1) = a_coef/(f_rad(x,1)+.0001)^5*exp(b_coef/(f_rad(x,1)+.0001)^4);
        end


    % Calculate the data for swell height

        Swell_Height = 2*(sum(Swell_Energy_Density_M0)^(.5))/.707;


    % Calculate the Swell Length

        Swell_Length = sum(Swell_Energy_Density_M0)/sum(Swell_Energy_Density_M1); % Divide the sums of the swell energy density M0 ordinate by the M1 Ordinate

        Swell_Mean_Period = 1/Swell_Length; % Get the swell mean peirod
        Swell_Mean_Period = 1/Swell_Mean_Period;
        Swell_Length = 5.12*Swell_Mean_Period^2; % Use the equation to get the swell length


    % Calculate the Swell Direction by finding the mean direction at the max
    % swell energy density point

        Swell_max = max(Swell_Energy_Density);
            for x = 1:len_WaveFreq
                if Swell_Energy_Density(x,1)==Swell_max
                   Swell_Direct = Wave_Mean_Direction_Deg_True(x,1);
                end
            end


    % Calculate the Swell Period
    Swell_max = max(Swell_Energy_Density);
        for x = 1:len_WaveFreq
            if Swell_Energy_Density(x,1)==Swell_max
               Swell_Period = Wave_Period(x,1);
            end
        end


    % Calculate the Swell 1/10th height

        Swell_110_Height = 5.1*sum(Swell_Energy_Density_M0)^(.5);


    % Calculate the swell 1/100 height 

        Swell_1100_Height = 6.5*sum(Swell_Energy_Density_M0)^(.5);


    % Calculate the Estimated Surf Height

        Estimated_Surf_Height = Swell_Height/.645;
        Estimated_Surf_Height_str = num2str(round(Estimated_Surf_Height,1));% Converts to string


    % Calculate the estimated surf period

        Estimated_Surf_Period = Swell_Period*.585;
        Estimated_Surf_Period_str = num2str(round(Estimated_Surf_Period,1));% Converts to string


    % Caluclate the Sea Sig Height

        Sea_Significant_Height = ((sum(Sea_Energy_Density_M0))^.5) * 4;


    % Calculate the Average length

        Average_Length = sum(Sea_Energy_Density_M1)/sum(Sea_Energy_Density_M0);
        Average_Length = 1/Average_Length;
        Average_Length = 5.12*Average_Length^2;


    % Calculate the sea modal time

        Sea_Energy_Density_max = max(Sea_Energy_Density);
            for x = 1:len_WaveFreq
                if Sea_Energy_Density(x,1)==Sea_Energy_Density_max
                   Sea_Modal_Period = Wave_Period(x,1);
                end
            end


    % Calculate the Sea Direction  

        Sea_Energy_Density_max = max(Sea_Energy_Density);
        for x = 1:len_WaveFreq
            if Sea_Energy_Density(x,1)==Sea_Energy_Density_max
               Sea_Direction_True = Wave_Mean_Direction_Deg_True(x,1);
            end
        end

        Sea_Direction_True_str = num2str(round(Sea_Direction_True,1));


    % Add the magnetic declination

        Sea_Direction_Magnetic = Sea_Direction_True + Magnetic_Declination;
        Sea_Direction_Magnetic_str = num2str(round(Sea_Direction_Magnetic,1));


    % Calculate the Sea Mean Period

        Sea_Mean_Period = sum(Sea_Energy_Density_M0);
        Sea_Mean_Period = sum(Sea_Energy_Density_M1)/Sea_Mean_Period;
        Sea_Mean_Period = 1/Sea_Mean_Period;


    % Calculate the Sea 1/10 height.

        Sea_110_Height = 5.1*sum(Sea_Energy_Density_M0)^(.5);


    % Calculate the sea 1/100 height 

        Sea_1100_Height = 6.5*sum(Sea_Energy_Density_M0)^(.5);


    % Calculate the total sig height using the Smoothed total energy density M0
    % ordinate

        Calculated_Total_Significant_Height = 4*sum(Smoothed_Total_Energy_Density_M0)^(.5);


    % Find the Wave Length Warning

        if Sea_Modal_Period < 7
            Wave_Length_Warning = 'Synchronous Waves';
        else
            Wave_Length_Warning = 'none';
        end


    % Find the Heading Warning

        if abs(Sea_Direction_True-Swell_Direct)<40
            Heading_Warning = 'Swell & Sea Same Dir';
        else
            Heading_Warning = 'none';
        end


    % Calculate the length adjusted PM Height

        Length_Adjusted_PM_Height = sum(PM_M1)/sum(PM_M0); 
        Length_Adjusted_PM_Height = 1/Length_Adjusted_PM_Height;
        Length_Adjusted_PM_Height = 5.12*Length_Adjusted_PM_Height^2; % These factors were taken from the excel sheet DataWellV08
        Length_Adjusted_PM_Height = (Length_Adjusted_PM_Height/Average_Length)^(.5);

            if Length_Adjusted_PM_Height > 1
                Length_Adjusted_PM_Height = 1;
            end

            if strncmpi(Wave_Length_Warning,'none',4) == 1 % if the wave length warning is none then the following equation is to be used to calculate the length adjusted PM height
    %             Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Calculated_Total_Significant_Height;
                Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Sea_Significant_Height;
            else
    %             Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Calculated_Total_Significant_Height * 1.05;
                Length_Adjusted_PM_Height = Length_Adjusted_PM_Height * Sea_Significant_Height * 1.05;
            end


    % If the heading warning is none then the length adjust pm height is
    % multiplied by 1.025

        if strncmpi(Heading_Warning,'none',4) == 0
            Length_Adjusted_PM_Height = Length_Adjusted_PM_Height*1.025;
        end


    % Convert Length Adjusted PM to string for later use    

        Length_Adjusted_PM_Height_str = num2str(round(Length_Adjusted_PM_Height,2));% Converts to string
        Length_Adjusted_PM_Height_string = string(Length_Adjusted_PM_Height_str);

    % Find the Sea State Level
    % The sea state level is dependent on the length adjusted PM Height

        if (.1 < Length_Adjusted_PM_Height) &&  (Length_Adjusted_PM_Height < 1)
            Sea_State_Level = 'Low ';
        elseif (1 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 1.8)
            Sea_State_Level = 'Mid ';
        elseif (1.8 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 2.1)
            Sea_State_Level = 'High ';
        elseif (2.1 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 3)
            Sea_State_Level = 'Low ';
        elseif (3 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 3.7)
            Sea_State_Level = 'Mid ';
        elseif (3.7 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 4.5)
            Sea_State_Level = 'High ';
        elseif (4.5 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 10)
            Sea_State_Level = 'Survival Range ';
        else
            Sea_State_Level = 'Do Not Operate ';
        end


    % Find the Sea State number
    % The sea state number is dependent on the significant total height and the
    % the length adjusted PM height

        if Calculated_Total_Significant_Height < .1
            Sea_State_Level_Number = '0';
        elseif Length_Adjusted_PM_Height < .7 
            Sea_State_Level_Number = '1';
        elseif (.7 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 2.1)
            Sea_State_Level_Number = '2';
        elseif (2.1 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 4.5)
            Sea_State_Level_Number = '3';
        elseif (4.5 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 8.4)
            Sea_State_Level_Number = '4';
        elseif (8.4 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 13.8)
            Sea_State_Level_Number = '5';
        elseif (13.8 < Length_Adjusted_PM_Height) && (Length_Adjusted_PM_Height < 21.1)
            Sea_State_Level_Number = '6';      
        else
            Sea_State_Level_Number = 'Greater than 6';
        end

        Sea_State = [Sea_State_Level Sea_State_Level_Number];
        Sea_State_short = string(Sea_State);


        length_adjusted_table.Length_Adjusted_PM_Height(m) = round(Length_Adjusted_PM_Height,1);
        Sea_State_table.Sea_State(m) = Sea_State_short;
        Estimated_Surf_Height_table.Estimated_Surf_Height(m) = round(Estimated_Surf_Height,1);
        Estimated_Surf_Period_table.Estimated_Surf_Period(m) = round(Estimated_Surf_Period,1);
        
    end

    calculations_table = [length_adjusted_table Sea_State_table Estimated_Surf_Height_table Estimated_Surf_Period_table];
    
    buoy_data_subset = [buoy_data_subset calculations_table]; % Add to the subset table
    
    toc % Output time up to this point in seconds
    
    
    
%% SET UP FIGURE FOR PLOTTING
% This will figure out what the max screensize is on the computer running
% this.  It will make the figure take up the entire screen instead of just
% making a standard size (very small) plot that has to be maximized to
% read.  Also, this will make sure when it saves a PNG file of the plot,
% it's as big as it can be.

    if create_plots == true
        
        disp(' ') % Adds a space in the command window 
        disp('Creating plots...')
        tic

        screen_size = get(0,'ScreenSize');
        screen_size_original = screen_size;
        screen_size(1) = 10; % Place on left side of screen
        screen_size(2) = (screen_size_original(4) - 840); % Place at top of screen
        screen_size(3) = 630; % Window width
        screen_size(4) = 750; % Window height
        
        fig = figure;
        fig.Units = 'pixels';
        fig.PaperUnits = 'inches';
        fig.PaperType = 'usletter';
        fig.PaperOrientation = 'portrait';
        fig.Position = [screen_size(1) screen_size(2) screen_size(3) screen_size(4)]; % Left Bottom Width Height
        set(0,'defaultfigurecolor',[1 1 1]);
        
        
        % MAKE SUBPLOTS AND ADD DATA
        if use_local_time == true
            xvals = buoy_data_subset.utc_time - hours(utc_offset);
        else
            xvals = buoy_data_subset.utc_time;
        end
        
        
        if plot_one == 1
            yvals1 = buoy_data_subset.Significant_Wave_Height_Feet; % define y-variable for Subplot 1
        else
            yvals1 = buoy_data_subset.Length_Adjusted_PM_Height;
        end
        
        if plot_two == 1
            yvals2 = buoy_data_subset.Peak_Wave_Period_Seconds; % define y-variable for Subplot 2
        else
            yvals2 = buoy_data_subset.Average_Wave_Period_Seconds; % define y-variable for Subplot 2
        end

        if plot_three == 1
            yvals3 = buoy_data_subset.Peak_Wave_Direction_Deg_True; % define y-variable for Subplot 3
        else
            yvals3 = buoy_data_subset.Length_Adjusted_PM_Height; % define y-variable for Subplot 3
        end
            
        % FIRST SUBPLOT
        subplot(12,1,[5,6]);
        plot(xvals,yvals1,'DatetimeTickFormat','dd','color','b')
        hold on;
       % plot(xvals,buoy_data_subset.Length_Adjusted_PM_Height,'color','g')
        grid on
        grid minor
        
        if plot_one == 1
            title('{\color{blue}Significant Wave Height, ft}')
        else
            title('{\color{blue}Length Adjusted PM Height, ft}')
        end
        
        ylabel('Hs, ft')
        hold on;

        if plot_forecast == true
            plot(buoy_data_forecast.utc_time,buoy_data_forecast.Significant_Wave_Height_Feet,'--r')
            if plot_one == 1
                title('{\color{blue}Significant Wave Height, ft}   -   {\color{red}Forcast Significant Wave Height, ft}')
            else
                title('{\color{blue}Length Adjusted PM Height, ft}   -   {\color{red}Forcast Significant Wave Height, ft}')
            end  
        end

        % Add upper/lower marker lines on wave height plot
        if plot_boundary == true 
            
            t1 = length(xvals);
            x1 = linspace(1,t1,t1);
            y1a = ones(1,t1) * min_boundary_line;
            y2a = ones(1,t1) * max_boundary_line;
            
            hold on;
            plot(xvals,y1a,'--k');
            hold on;
            plot(xvals,y2a,'--k');            
            
            if plot_forecast == true
                t2 = length(buoy_data_forecast.utc_time);
                y1b = ones(1,t2) * min_boundary_line;
                y2b = ones(1,t2) * max_boundary_line;
                hold on;
                plot(timedt_f,y1b,'--k');
                hold on;
                plot(timedt_f,y2b,'--k');
            end
        end

        % SECOND SUBPLOT
        subplot(12,1,[8,9])
        plot(xvals,yvals2,'DatetimeTickFormat','dd','color','b')
        hold on;
        grid on
        grid minor
        if plot_two == 1
            title('{\color{blue}Peak Wave Period (sec)}')
        else
            title('{\color{blue}Average Wave Period (sec)}')
        end
        
        ylabel('Tp, sec')
        hold on;

        if plot_forecast == true
            plot(buoy_data_forecast.utc_time,buoy_data_forecast.Peak_Wave_Period_Seconds,'--r')
            if plot_two == 1
                title('{\color{blue}Peak Wave Period (sec)}   -   {\color{red}Forcast Peak Wave Period (sec)}')
            else
                title('{\color{blue}Average Wave Period (sec)}   -   {\color{red}Forcast Average Wave Period (sec)}')
            end  
        end 

        % THIRD SUBPLOT
        subplot(12,1,[11,12])
        if plot_three == 1
            plot(xvals,yvals3,'o','DatetimeTickFormat','dd',... % Use 'plot' function (to be able to use datetime array) but set markers as scatter points
                'MarkerSize',2,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
            title('{\color{blue}Peak Wave Direction (deg true)}')
            ylim([0 360])
            ylabel('Dp, deg')
        else
            plot(xvals,yvals3,'DatetimeTickFormat','dd','color','b')
            title('{\color{blue}Length Adjusted PM Height, ft}')
            ylabel('Hs, ft')
        end
        
        hold on;
        grid on

        xlabel('Day')
        hold on;

        if plot_forecast == true
            
            if plot_three == 1
                plot(buoy_data_forecast.utc_time,Peak_Wave_Direction_Forecast,'--r')
                title('{\color{blue}Peak Wave Direction (deg true)}   -   {\color{red}Forcast Peak Wave Direction (deg true)}')
            else
                title('{\color{blue}Length Adjusted PM Height, ft}')
            end  
            
        end

        if plot_forecast == true
            buoytitle = [buoy_name,' Buoy #',buoy_number,' + Forecast'];
        else
            buoytitle = [buoy_name,' Buoy #',buoy_number];
        end
        
        suptitle(buoytitle) % Set main plot title

        
        % Create Table on Plot
        % I broke the values out this way so they could be rearranged or
        % changed easier than it was previously.
        row1name = 'Local Time of Study';
        row1 = Wave_Time_Local_Now_str;
        
        row2name = 'Significant Wave Height';
        row2 = [Significant_Wave_Height_Feet_Now_str,' feet'];
        
        row3name = 'Length Adjusted PM Height';
        row3 = [Length_Adjusted_PM_Height_str,' feet'];
        
        row4name = 'Sea State';
        row4 = Sea_State;
        
        row5name = 'Peak Wave Direction';
        row5 = [Peak_Wave_Direction_Deg_True_Now_str,' degrees true'];
        
        row6name = 'Peak Wave Period';
        row6 = [Peak_Wave_Period_Seconds_Now_str,' seconds'];
        
        row7name = 'Sea Surface Temperature';
        row7 = [Sea_Surface_Temp_DegF_Now_str,' degress F'];
        
        row8name = 'Estimated Surf Height';
        row8 = [Estimated_Surf_Height_str,' feet'];
        
        row9name = 'Estimated Surf Period';
        row9 = [Estimated_Surf_Period_str,' seconds'];
        
        % Create table from the above values
        table_data = {row1name row1;row2name row2;row3name row3;row4name row4;...
            row5name row5;row6name row6;row7name row7;row8name row8; row9name row9};

        t1 = uitable('Data',table_data);
        t1.Units = 'pixels';
        t1.FontSize = 8;
        t1.FontWeight = 'bold';
        t1.Position = [80 520 1 1]; % Left Bottom Width Height
        t1.ColumnName = ({}); % Blank column names
        t1.RowName = ({});
        t1.ColumnWidth = {244};
        t1.Position(3:4) = t1.Extent(3:4);
        t1.ColumnFormat = {'char' 'char' 'char' 'char' 'char' 'char' 'char' 'char' 'char'};
        t1.RowStriping = 'on'; % Alternating row shading
        t1.Parent = fig;
                
        toc
        
        
        
        %% Pause and Save Plot to Disk
        % This is to give a slight delay to allow proper formatting and display of
        % the plot before actually writing it to disk.
        
        disp(' ') % Adds a space in the command window 
        disp('Saving plot to output folder...')
        tic

        filename = Wave_Time_Local_Now;
        filename = datestr(filename,'yyyy-mm-dd_HH-MM-SS');
        
        picture = (strcat(output_folder,'\',filename,'.png'));

        pause(1); % Gives the figure a moment to open all the way before saving.                 
        saveas(gcf,picture) % Much higher quality than hgexport, better for printing the image
        disp(picture)
        toc
        
        
        
    end % End of Plotting Section
            
    

%% Upload image to server location
% Upload to "lab-data.org" so it can be retrieved via Sharepoint.  

    if upload_via_ftp == true
        tic
        disp('Uploading image to site...');
        
        close(figure(1))
        ftpobj = ftp('ftp.lab-data.org','matlab','rLbBV{Y3t3@&'); % Create FTP object (server, user, password)
        mput(ftpobj,'C:\wavedata\wavedata.png'); % This uploads the selected file via the created ftp object
        close(ftpobj); % This closes the FTP connection
        
        toc
    end



%% Save Table To Disk
% This will create a new table with just the values we want to save to a
% XLS / MAT file.

    if save_table == true

        disp(' ') % Adds a space in the command window 
        disp('Saving plot to output folder...')
        tic

        output_table = buoy_data_subset;
        output_table.Wave_Mean_Direction_Deg_True = [];
        output_table.Wave_Energy_Raw = [];
        output_table.Significant_Wave_Height_Meters = [];
        output_table = movevars(output_table,'Length_Adjusted_PM_Height','After','utc_time');
        output_table = movevars(output_table,'Significant_Wave_Height_Feet','After','Length_Adjusted_PM_Height');
        output_table = movevars(output_table,'Peak_Wave_Period_Seconds','After','Significant_Wave_Height_Feet');

        filename = Wave_Time_Local_Now;
        filename = datestr(filename,'yyyy-mm-dd_HH-MM-SS');
        
        check = any(strcmp('xls_file',save_table_format));
        if check == 1
            savefiletable = strcat(output_folder,'\',filename,'.xls');
            writetable(output_table,savefiletable); % Save table to disk
            message = 'Saved data to disk as a XLS file...';
            disp(message)
            disp(savefiletable)
        else
            savefilemat = strcat(output_folder,'\',filename,'.mat');
            save(savefilemat,'output_table');
            message = 'Saved data to disk as a MAT file...';
            disp(message)
            disp(savefilemat)
        end

        toc

    end



%% Output Messages
% Displays a couple values in the command window for ease of use.
    
    output_message_1 = ['The most recent data is from ',Wave_Time_Local_Now_str,'.'];% Creates message
    output_message_2 = ['The length adjusted PM height is ',Length_Adjusted_PM_Height_str,' feet.'];% Creates message   
    output_message_3 = ['The estimated surf height is ',Estimated_Surf_Height_str,' feet with an estimated period of ',Estimated_Surf_Period_str,' seconds.'];% Creates message
    
    disp(' ') % Adds a space in the command window
    disp(output_message_1) % Displays message
    disp(output_message_2) % Displays message
    disp(output_message_3) % Displays message

    
    
%% Housekeeping / Timing
% This section does any housekeeping required, such as display timing
% messages to the command window, or clean up the workspace if required.

    disp(' ') % Adds a space in the command window 
    programend = etime(clock,programstart); % Estimated time for whole program execution
    programendstr = num2str(programend); % Converts time to string
    processtime=['Total program execution time was ',programendstr,' seconds.']; % Creates message
    disp(processtime) % Displays message
    disp(' ') % Adds a space in the command window 
    

    
    

