%% Read NCS data
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [data, dateVector] = readNcsData(dataPath,fileName)
% function readHxData(dataPath,fileName)
%% Adding data path
% dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
%     'Hexoskin\Data\2'];
% fileName = 'freq2G_1';

%% Reading file
load([dataPath,fileName,'.mat']);
data = ampPh;

% data(:,1) = abs(flipData(1) * data(:,1));
% data(:,2) = abs(flipData(2) * data(:,2));

%% Storing the starting time of the data
formatNcsTime = 'dd-mmm-yyyy HH:MM:SS';
dateVector = datevec(dateTime, formatNcsTime);

