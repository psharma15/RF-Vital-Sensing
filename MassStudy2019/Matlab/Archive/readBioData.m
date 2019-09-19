% This code reads BIOPAC data. 
% May 02, 2019
% Pragya Sharma, ps847@cornell.edu


function [data,dateVector] = readBioData(dataPath,fileName)
% function readNcsData(dataPath,fileName)

%% Reading file
load([dataPath,fileName,'.mat'],'data');

% Storing data start time [Y M D H MI S] format
dateVector = bioTime(fileName);

end
          

