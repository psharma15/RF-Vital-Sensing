% Pragya Sharma, 12 July 2018
% Edited, 08 March 2019: 2 channel
% ps847@cornell.edu

function [convertedData,ch12345,subjectInfo] = saveTDMStoMAT()
% SAVETDMSTOMAT specifies data path, and file name to be converted from
% TDMS to mat, and saveFile name is the name of the saved mat file
% It uses convertTDMS.m to convert TDMS to mat.

dataPath = 'E:\NCS\HumanStudyMarch2019\Data\Case4\';
fileName = '0427_123254Routine 7';
saveFileName =  fileName;
[convertedData,ch12345,subjectInfo] = callConvert(dataPath,fileName);
saveFileName = [dataPath,saveFileName,'.mat'];
save(saveFileName,'ch12345','subjectInfo');
end

function [convertedData, ch12345, subjectInfo] = callConvert(dataPath,fileName)
% CALLCONVERT takes data path and file name as input. It calls 
% dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\BreathPattern_Cough_Speak\Data\Yuna\';
% fileName = 'freq2G_1v1h';
% saveFileName = 'freq2G_1v1h';
filePathName = [dataPath,fileName,'.tdms'];
[convertedData,~,~,~,~] = convertTDMS(false,filePathName);

% The order is [amp_tx1rx1 amp_tx2rx2 ph_tx1rx1 ph_tx2rx2]
ch12345 = [convertedData.Data.MeasuredData(3).Data, ...
        convertedData.Data.MeasuredData(4).Data,...
        convertedData.Data.MeasuredData(5).Data,...
        convertedData.Data.MeasuredData(6).Data,...
        convertedData.Data.MeasuredData(7).Data];
    
subjectInfo = convertedData.Data.Root.Property; 
end

