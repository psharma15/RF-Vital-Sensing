%% Read NCS data
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [data, dataSampRate, dateVector] = readNcsData(dataPath,dataNum)
% function readHxData(dataPath,dataNum)
%% Adding data path
% dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
%     'Hexoskin\Data\2'];
% fileName = 'freq2G_1';

%% Reading file
addpath(dataPath);

[fileName, dateTime] = ncsFileInfo(dataNum);

load([fileName,'.mat']);
data = eval(fileName);
dataSampRate = 500; % NCS sampling rate in Hz

% data(:,1) = abs(flipData(1) * data(:,1));
% data(:,2) = abs(flipData(2) * data(:,2));

%% Storing the starting and end time of the data
formatNcsTime = 'dd-mm-yyyy HH:MM:SS.FFF';
dateVector = datevec(dateTime, formatNcsTime);

t = 0:1/dataSampRate:((length(data(:,1))/dataSampRate)-(1/dataSampRate));
% figure
% title(['Data collection (started at: ',num2str(dateVector(4)),':',...
%     num2str(dateVector(5)),':',num2str(uint32(dateVector(6)))]);
% ax(1) = subplot(2,1,1);
% plot(t,data(:,1)); ylabel('NCS Amp')
% ax(2) = subplot(2,1,2);
% plot(t,data(:,2)); ylabel('NCS Ph')
% linkaxes(ax,'x')

%% 
rmpath(dataPath);
