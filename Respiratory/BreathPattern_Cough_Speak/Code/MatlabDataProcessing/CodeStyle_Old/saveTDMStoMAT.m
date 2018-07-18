% Pragya Sharma, 12 July 2018
% ps847@cornell.edu

function [dataTime,ampPh] = saveTDMStoMAT()

dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
    'BreathPattern_Cough_Speak\Data\Pragya\Jul16\Jul16\'];
fileName = 'freq2G_1v2h';
saveFileName =  'freq2G_2v2h';
[convertedData,ampPh] = callConvert(dataPath,fileName);
dataTime = convertedData.Data.MeasuredData(5).Property;
dataTime = dataTime(8).Value;
saveFileName = [dataPath,saveFileName,'.mat'];
save(saveFileName,'ampPh','dataTime');
end

function [convertedData,ampPh] = callConvert(dataPath,fileName)
% SAVETDMSTOMAT takes input data path, and file name to be converted from
% TDMS to mat, and saveFile name is the name of the saved mat file
% It uses convertTDMS.m to convert TDMS to mat.
% dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\BreathPattern_Cough_Speak\Data\Yuna\';
% fileName = 'freq2G_1v1h';
% saveFileName = 'freq2G_1v1h';
filePathName = [dataPath,fileName,'.tdms'];
[convertedData,~,~,~,~]=simpleConvertTDMS(false,filePathName);

ampPh = [convertedData.Data.MeasuredData(4).Data, convertedData.Data.MeasuredData(5).Data];
end

