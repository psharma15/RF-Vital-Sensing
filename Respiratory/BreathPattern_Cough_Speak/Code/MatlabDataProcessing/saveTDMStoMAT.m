% Pragya Sharma, 12 July 2018
% Edited, 17 July 2018
% ps847@cornell.edu

function [dateTime,ampPh] = saveTDMStoMAT()
% SAVETDMSTOMAT specifies data path, and file name to be converted from
% TDMS to mat, and saveFile name is the name of the saved mat file
% It uses convertTDMS.m to convert TDMS to mat.

dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
    'BreathPattern_Cough_Speak\Data\Pragya\Jul16\'];
fileName = 'freq2G_1v2h';
saveFileName =  'freq2G_2v2h';
[convertedData,ampPh] = callConvert(dataPath,fileName);
dateTime = convertedData.Data.MeasuredData(5).Property(8).Value;
saveFileName = [dataPath,saveFileName,'.mat'];
save(saveFileName,'ampPh','dateTime');
end

function [convertedData,ampPh] = callConvert(dataPath,fileName)
% CALLCONVERT takes data path and file name as input. It calls 
% dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\BreathPattern_Cough_Speak\Data\Yuna\';
% fileName = 'freq2G_1v1h';
% saveFileName = 'freq2G_1v1h';
filePathName = [dataPath,fileName,'.tdms'];
[convertedData,~,~,~,~] = convertTDMS(false,filePathName);

ampPh = [convertedData.Data.MeasuredData(4).Data, convertedData.Data.MeasuredData(5).Data];
end

