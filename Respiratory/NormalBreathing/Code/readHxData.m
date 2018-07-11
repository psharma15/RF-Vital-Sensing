% Read CSV files from Hexoskin - plot or observe the results. 
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [data, dataSampRate, dateVector] = readHxData(dataPath,dataNum)
% function readHxData(dataPath,dataNum)
% IMP: Remember issue with the way Hx data is saved. It saves previous data
% records in the same file as new data. 
% Input:
%   dataPath: Path to data.
%   dataNum: number of data which is arranged in 'folderName'.
% Output:
%   data: data as requested.
%   dataSampRate: Sample rate of the data. It is 0 if data is uneven
%   sampled.
%   dateVector: Absolute time corresponding to each sample in data,
%   arranged in columns representing [yyyy mm dd HH MM SS.FFFF] hours, mins
%   and seconds are in 24 hours time format.

%%
addpath(dataPath);

folderName = {'respiration_abdominal','respiration_thoracic',...
              'ECG_I',...
              'tidal_volume_raw','tidal_volume_adjusted',...
              'breathing_rate',...
              'minute_ventilation_raw','minute_ventilation_adjusted',...
              'inspiration','expiration',...
              'heart_rate','RR_interval','NN_interval',...
             };
sampleRates = [128, 128,    ...
               256,         ...
                 1,   1,    ...
                 1,         ...
                 1,   1,    ...
                 0,   0,    ...
                 1,   0,  0,...  
               ]; % Hz
dataName = folderName{dataNum}; % 'respiration_thoracic','ECG_1'
dataSampRate = sampleRates(dataNum);
fprintf(['Reading Hexoskin data,', dataName,'\n']);

%%
data = csvread([dataName,'.csv'],1,0);
posixTime = data(:,1);
% Converting: number of seconds since 1-Jan-2010 00:00:00
% There is some time offset (4 hours addition - check)
% Leap seconds are difference between atomic clock and UTC. Since 1 Jan
% 2010, 3 leap seconds have been added till March 2018.
convertedTime = datetime(posixTime,'ConvertFrom','epochtime',...
                'Epoch','2010-01-01','TimeZone','America/New_York',...
                'Format','d-MMM-y HH:mm:ss Z')-hours(4)+seconds(3);
% Storing the starting and end time of the data
dateVector = datevec(convertedTime);

% if dataSampRate ~= 0
%     t = 0:1/dataSampRate:((length(data(:,1))/dataSampRate)-(1/dataSampRate));
%     figure
%     plot(t,data(:,2))
%     title(dataName, 'Interpreter', 'none')
%     xlabel('Time (sec)')    
% end

% % Flip inhalation and exhalation
% data = abs(-data);
data = data(:,2);

%% 
rmpath(dataPath);
