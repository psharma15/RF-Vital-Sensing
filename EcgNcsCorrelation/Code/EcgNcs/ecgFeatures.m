%% This code finds features of ECG: PQRST peaks and their location.


dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Data\Mar03';
addpath(dataPath);
load('freq2G2.mat');
rmpath(dataPath);
ecgData = freq2G2(:,3);

fs = 512;
t = (0:length(ecgData)-1)/fs;

% Using peakdetect.m, which gives correct Q,R,T points.
peakDetectPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Code\EcgNcs\peakdetect_ECG_v4';
addpath(peakDetectPath);
view = length(ecgData)/fs;
[R_i,~,S_i,~,T_i,~,Q_i,~,heart_rate,buffer_plot] = peakdetect2(ecgData,fs,view);
rmpath(peakDetectPath);


% Plot figure
figure
plot(t,ecgData);
hold on
plot(t(R_i),ecgData(R_i),'*');

%% 