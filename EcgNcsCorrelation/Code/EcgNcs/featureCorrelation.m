%% This code is to identify correlated features and estimate similarity  
%  between different features of NCS and ECG. 

%% Reading data
dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Data\Mar03';
fileName = 'freq2G2';
tSynchStabilize = 10;
FsOld = 512;
FsNew = 250;
tDelay = 90e-3; % ECG and NCS delay

[ncsAmpData,ncsPhData,ecgData,t] = readEcgNcs(dataPath,fileName,tSynchStabilize,FsOld,FsNew,tDelay);
nSample = length(ncsAmpData);
if length(ecgData) ~= nSample
    fprintf('Ecg and Ncs have different length. \nStopping...\n');
    return
end

% close all; % Close all figures.

%% Perform Post Processing of NCS to view heartbeat
[ncsAmpFilt,ncsPhFilt,~,~,~] = postProcess(0.8,5,6,[ncsAmpData,ncsPhData],FsNew,[1,1]);
[ecgFilt,~,~,~,~] = postProcess(0.9,80,81,[ecgData, zeros(nSample,1)],FsNew,[1,1]);

% Plotting scaled versions of NCS and ECG
ncsAmp = 0.025e4*ncsAmpFilt;
ncsPh = ncsPhFilt;
ecg = 0.25e-3*ecgFilt;

% Removing first 6 seconds of data as filter output takes time to stabilize
tFiltStabilize = 6;
ncsAmp = ncsAmp(tFiltStabilize*FsNew+1:end);
ncsPh = ncsPh(tFiltStabilize*FsNew+1:end);
ecg = ecg(tFiltStabilize*FsNew+1:end);
t = t(1:end-(tFiltStabilize*FsNew));
nSample = nSample - 6*FsNew;

figure
grid on
yyaxis left
plot(t,ncsAmp); 
ylabel('NCS')
yyaxis right
plot(t,ecg); 
ylabel('ECG');

%% Features of ECG (PQRST)
% Using peakdetect.m, which gives correct Q,R,T points.
peakDetectPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Code\EcgNcs\peakdetect_ECG_v4';
addpath(peakDetectPath);
% Using data at 250Hz. Poorer results for high and low (very poor) rates.
ecgData = ecgData(tFiltStabilize*FsNew+1:end); 

[R_i,~,S_i,~,T_i,~,Q_i,~,heart_rate,buffer_plot] = peakdetect2(ecgData,FsNew,length(ecgData)/FsNew);
rmpath(peakDetectPath);

%% Features of NCS 
