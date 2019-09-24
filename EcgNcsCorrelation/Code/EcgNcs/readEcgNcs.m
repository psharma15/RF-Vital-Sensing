%% This code reads ECG-NCS data:

function [ncsAmpTrunc,ncsPhTrunc,ecgDelay,tTrunc] = readEcgNcs(dataPath,fileName,tStabilize,FsOld,FsNew,tDelay)
%% READECGNCS function reads ECG and NCS amp and phase data and generates synchronized data at different sample rate.
% dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Data\21_22Feb2018';
% fileName = 'data5' ('data5.mat' should exist)
% tStabilize = 40 (Time to allow ECG and NCS to synchronize and become stable: seconds)
% FsOld = 512 (Acquired data sampling rate: Hz)
% FsNew = 200 (Processed data sampling rate: Hz)
% tDelay = 90e-3 (Constant time by which ECG lags: seconds)

%%
% Addpath
addpath(dataPath);

load([fileName,'.mat']);
data = eval(fileName);

% Allow data to stabilize
data = data(FsOld*tStabilize+1:end,:);

ncsUnfiltAmp = data(:,1);
ncsUnfiltPh = data(:,2);
ecg = data(:,3);

nSample = length(ncsUnfiltAmp);
idx = 1:nSample;
t = ((idx-1)/FsOld)';
df = FsOld/nSample;
f = -FsOld/2:df:FsOld/2 - df;

%
figure
nFigRow = 3;
nFigCol = 1;
ax1(1) = subplot(nFigRow,nFigCol,1);
plot(t,ncsUnfiltAmp,'k'); grid on;
ax1(2) = subplot(nFigRow,nFigCol,2);
plot(t,ncsUnfiltPh,'r'); grid on;
ax1(3) = subplot(nFigRow,nFigCol,3);
plot(t,ecg,'k'); grid on;

linkaxes(ax1(:),'x');

rmpath(dataPath);

%In test1, we are loosing NCS data as lv is not updated, and also
%because NCS is slow (?). That's why, when data is saved, zeros are
%inserted in NCS if no new data is updated. That needs to be corrected
%and NCS needs to be resampled. Also, ECG needs to be resampled at a
%lower rate. Both needs to be filtered. 
idxZeroCorrect = (ncsUnfiltAmp ~= 0);


% This would introduce little error, but lead to more issues in
% terms of timing and resampling later.
idxZeroCorrect(1) = 1;
idxZeroCorrect(end) = 1;

ncsAmpZeroCorrected = ncsUnfiltAmp(idxZeroCorrect);
ncsPhZeroCorrected = ncsUnfiltPh(idxZeroCorrect);
tZeroCorrected = t(idxZeroCorrect);

figure
nFigRow = 3;
nFigCol = 1;
ax2(1) = subplot(nFigRow,nFigCol,1);
plot(tZeroCorrected,ncsAmpZeroCorrected,'k'); grid on;
ax2(2) = subplot(nFigRow,nFigCol,2);
plot(tZeroCorrected,ncsPhZeroCorrected,'r'); grid on;
ax2(3) = subplot(nFigRow,nFigCol,3);
plot(t,ecg,'k'); grid on;

linkaxes(ax2(:),'x');

[ncsAmpResampled, tResampledNcs] = resample(ncsAmpZeroCorrected,tZeroCorrected,FsNew);
[ncsPhResampled, ~] = resample(ncsPhZeroCorrected,tZeroCorrected,FsNew);
[ecgResampled, tResampledEcg] = resample(ecg,t,FsNew); % This low passes by default

tResampled = tResampledNcs; %or ECG
if tResampledEcg(end) < tResampledNcs(end)
    tResampled = tResampledEcg;
    ncsAmpResampled = ncsAmpResampled(1:length(tResampled));
    ncsPhResampled = ncsPhResampled(1:length(tResampled));
end
if tResampledNcs(end) < tResampledEcg(end)
    tResampled = tResampledNcs;
    ecgResampled = ecgResampled(1:length(tResampled));
end

figure
nFigRow = 3;
nFigCol = 1;
ax3(1) = subplot(nFigRow,nFigCol,1);
plot(tResampled,ncsAmpResampled,'k'); grid on;
ax3(2) = subplot(nFigRow,nFigCol,2);
plot(tResampled,ncsPhResampled,'r'); grid on;
ax3(3) = subplot(nFigRow,nFigCol,3);
plot(tResampled,ecgResampled,'k'); grid on;
linkaxes(ax3(:),'x');

% If delay is known: 90e-3 seconds, constant over multiple data sets and time.
ecgDelay = ecgResampled(floor(tDelay*FsNew+1):end);
ncsAmpTrunc = ncsAmpResampled(1:end-floor(tDelay*FsNew));
ncsPhTrunc = ncsPhResampled(1:end-floor(tDelay*FsNew));
tTrunc = tResampled(1:end-floor(tDelay*FsNew));

figure
nFigRow = 3;
nFigCol = 1;
ax4(1) = subplot(nFigRow,nFigCol,1);
plot(tTrunc,ncsAmpTrunc,'k'); grid on;
ax4(2) = subplot(nFigRow,nFigCol,2);
plot(tTrunc,ncsPhTrunc,'r'); grid on;
ax4(3) = subplot(nFigRow,nFigCol,3);
plot(tTrunc,ecgDelay,'k'); grid on;
linkaxes(ax4(:),'x');

end
