%% Read NCS data
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [apap,mic, dateVector, subjectInfo] = readNcsData(dataPath,fileName)
% function readNcsData(dataPath,fileName)
% The NCS data originally is in format [amp_ch1 amp_ch2 ph_ch1 ph_ch2]
% The output as shown is in format [amp_ch1 ph_ch1 amp_ch2 ph_ch2]

%% Reading file
load([dataPath,fileName,'.mat'],'ch12345','subjectInfo');
apap = [ch12345(:,1), ch12345(:,3), ch12345(:,2), ch12345(:,4)];

if size(ch12345,2) == 5
    mic = ch12345(:,5);
else
    mic = [];
end

% Storing data start time [Y M D H MI S] format.
% Check the format in case of any issue. MATLAB's dateVector format saves
% everything as integer type, so no fractional seconds. Here it is
% fractional seconds
dateVector = ([(subjectInfo(4).Value), ...
              (subjectInfo(5).Value), ...
              (subjectInfo(6).Value), ...
              (subjectInfo(7).Value), ...
              (subjectInfo(8).Value), ...
              subjectInfo(9).Value + subjectInfo(10).Value]);
          

