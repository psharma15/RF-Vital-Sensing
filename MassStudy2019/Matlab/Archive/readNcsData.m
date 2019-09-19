%% Read NCS data
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [apap,mic, dateVector, subjectInfo] = readNcsData(dataPath,fileName)
% function readNcsData(dataPath,fileName)

%% Reading file
load([dataPath,fileName,'.mat'],'ch12345','subjectInfo');
apap = [ch12345(:,1), ch12345(:,3), ch12345(:,2), ch12345(:,4)];
mic = ch12345(:,5);

% Storing data start time [Y M D H MI S] format.
% Check the format in case of any issue. MATLAB's dateVector format saves
% everything as integer type, so no fractional seconds. Here it is
% fractional seconds
dateVector = ([(subjectInfo(10).Value), ...
              (subjectInfo(11).Value), ...
              (subjectInfo(12).Value), ...
              (subjectInfo(13).Value), ...
              (subjectInfo(14).Value), ...
              subjectInfo(15).Value + subjectInfo(16).Value]);
          

