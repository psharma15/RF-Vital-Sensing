% This function estimates HRV features
% Pragya Sharma
% ps847@cornell.edu
% Date created: June 25 2019
% nn: equivalent to rr, 'normal' rr
% sdnn: Standard deviation of all NN intervals
% rmssd: Square root of the mean of the squares of successive differences 
% between adjacent NN intervals

function [hrvNcs,hrvEcg,pkNcs,pkEcg] = hrvFeatureEst(dataNcs,dataEcg,fs,opts)

if size(dataNcs,1) > 1
    [hrvNcs,pkNcs] = hrvFeatureEstNcs(dataNcs,fs(1),opts);
else
    hrvNcs = []; pkNcs = [];
end

if size(dataEcg,1) > 1
    [hrvEcg,pkEcg] = hrvFeatureEstEcg(dataEcg,fs(2),opts);
else
    hrvEcg = []; pkEcg = [];
end

end

function [hrvFeature,pk] = hrvFeatureEstNcs(dataNcs,fs,opts)

if ~isfield(opts,'tWinHR')
    opts.tWinHR = 5;
    fprintf('Default HR estimation window: %3.2f\n',opts.tWinHR);
end
if ~isfield(opts,'delPkCalibTime')
    % Make sure this time interval has no wrong beat
    opts.delPkCalibTime = [0 30];
    fprintf('Heartbeat peak amplitude calibration window: [%2.1f, %2.1f]s. \n',opts.delPkCalibTime);
end
if ~isfield(opts,'pkAmpRelRejThresh')
    % This is minimum peak value of a beat relative to mean peak amplitude,
    % identified from delPkCalibTime, below which it will be rejected.
    opts.pkAmpRelRejThresh = 0.6; 
    fprintf('Heartbeat peak amplitude relative rejection threshold: %f\n',opts.pkAmpRelRejThresh);
end
if ~isfield(opts,'tRRthresh')
    % Minimum and maximum time between two peaks
    % This could be made more statistical, for now it is fixed
    opts.tRRthresh = [0.4,1.2]; % 60./[0.4, 1.2] = [150, 50] BPM
    fprintf('Min-Max RR time: [%2.1f, %2.1f]s\n',opts.tRRthresh);
end
if ~isfield(opts,'harmNum')
    % If NCS is harmonic, select this clearly. Default not harmonic
    opts.harmNum = 1; % 1 is fundamental
end

t = (0:(length(dataNcs)-1))/fs; 
t = t(:);
tPkCalibWindow = zeros(1,2);
tPkCalibWindow(1) = t(opts.delPkCalibTime(1) * fs + 1);
tPkCalibWindow(2) = t(opts.delPkCalibTime(2) * fs + 1);

% -------------------------------------------------------------------------
% Minima and maxima detection on NCS thorax and abdomen data
% -------------------------------------------------------------------------
pk = findMaxMin(dataNcs,fs,opts);

if pk(1).ind(1) == 1
    % Start with a minima
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 0
    % End with a maxima
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end


pkMax = pk(1).idx(pk(1).ind == 1);
pkMin = pk(1).idx(pk(1).ind == 0);

if length(pkMax) ~= length(pkMin)
    fprintf('length(pkMin) ~= length(pkMax) in hrvFeatureEst. Returning.\n');
    hrvFeature = [];
    return;
end

% This is the change in minima to maxima for each beat, assuming size of
% pkMin is same as size of pkMax
delPk =  abs(dataNcs(pkMax,1)-dataNcs(pkMin,1)); % Making it positive always
tDelPkMax = t(pkMax);

% Correction: Ignoring a beat if peak to peak is lower than some percentage
% of the mean peak amplitude.
delPkCalib = delPk((tDelPkMax >= tPkCalibWindow(1)) & (tDelPkMax <= tPkCalibWindow(2)));
meanDelPkCalib = mean(delPkCalib);
% The following beats satisfy the constraint
idxDelPk = delPk >= (meanDelPkCalib .* opts.pkAmpRelRejThresh);
pkMax = pkMax(idxDelPk);
pkMin = pkMin(idxDelPk);

%% RR Interval
% This is time difference between two maximas. This will need correction as
% we may have missed some peaks, and in above peak correction we might have
% ignored some true peaks as well.
rrInterval = t(pkMax(2:end)) - t(pkMax(1:end-1));
% tRR: Associating RR interval time at the second peak
tRR = t(pkMax(2:end)); 

% If harmonic of heartbeat, correct the interval
if opts.harmNum == 2
    rrInterval = rrInterval(2:2:end) + rrInterval(1:2:end-1);
    tRR = tRR(2:2:end);
end

% Correcting: RR interval limited by estimated HR in healthy seater person
idxRRinterval =  (rrInterval <= opts.tRRthresh(2)) & (rrInterval >= opts.tRRthresh(1));
rrInterval = rrInterval(idxRRinterval);
tRR = tRR(idxRRinterval);

% Since other correction uses next rrInterval that could be incorrect,
% editing that
nWinRRcorrect = 10; % Use last 10 RR intervals mean to correct
if length(rrInterval)< nWinRRcorrect*2
    nWinRRcorrect = 1;
    fprintf('Changing nWinRRcorrect.\n');
    
end
for i = 1:nWinRRcorrect
    if (rrInterval(i) < (1/1.3)*mean(rrInterval(i+1:nWinRRcorrect))) || (rrInterval(i) > 1.3*mean(rrInterval(i+1:nWinRRcorrect)))
        rrInterval(i) = mean(rrInterval(i+1:nWinRRcorrect));
    end
end
        

for i = nWinRRcorrect+1:length(rrInterval) - 1
    if (rrInterval(i) > 1.3* mean(rrInterval(i-(nWinRRcorrect):i-1))) || (rrInterval(i) < (1/1.3)* mean(rrInterval(i-(nWinRRcorrect):i-1)))
        rrInterval(i) = rrInterval(i-1); % Use the last interval
    end
end

% Correction: RR interval is not expected to change suddenly with time
% (RR(k) > 0.7?(RR(k ? 1)+RR(k + 1))) -> RR(k) = (RR(k ? 1)+RR(k + 1))/2

for i = 2:length(rrInterval)-1
    if (rrInterval(i) > (0.65*(rrInterval(i-1) + rrInterval(i+1))) || rrInterval(i) < (0.40*(rrInterval(i-1) + rrInterval(i+1))))
        rrInterval(i) = (rrInterval(i-1) + rrInterval(i+1))/2;
    end
end

fprintf('NCS: Mean RR interval: %3.2f s\n',mean(rrInterval));

%% Succesive RR information: In some cases we miss some beats, so the 
% rrInterval is not successive. That information is important for finding
% successive differences - and affects RMSSD, pNN50, etc.
tRRdiff = tRR(2:end) - tRR(1:end-1);
isNextSuccessive = tRRdiff <= opts.tRRthresh(2); % Is next RR interval successive or not

% Saving only successive differences, otherwise 0
sdRR = isNextSuccessive .* (rrInterval(2:end) - rrInterval(1:end-1));
sdRR = sdRR(sdRR~=0);

nn50 = sum(sdRR > 50e-3); % Number of successive differences greater than 50ms

%% Interpolated RR or Inter-Beat Interval (IBI)
% Selecting sampling frequency of 10 Hz.
rrResampleFs = 10;
tRRInterp = tRR(1):1/rrResampleFs:tRR(end);
rrIntervalInterp = interp1(tRR,rrInterval,tRRInterp);

figure
plot(tRR,rrInterval,tRRInterp,rrIntervalInterp);
title('NCS Interpolated RR interval');

%% Frequency analysis: PSD/STFT/wavelet

% Starting with FFT to compare entire baseline vs test durations
n = 2^nextpow2(length(rrIntervalInterp));
yRR = fft(rrIntervalInterp - mean(rrIntervalInterp),n);
f = rrResampleFs*(0:(n/2))/n;
P = abs(yRR).^2/n;

figure
plot(f,P(1:n/2+1)) 
title('Frequency spectrum NCS RR interval')

% LF power: 0.04 - 0.15 Hz, HF: 0.15 - 0.7 Hz
idxLF1 = find(f<=0.04,1,'last');
idxLF2HF1 = find(f>=0.15,1,'first');
idxHF2 = find(f>=0.7,1,'first');

LFpow = trapz(P(idxLF1:idxLF2HF1));
HFpow = trapz(P(idxLF2HF1:idxHF2));
LFHFratio = LFpow/HFpow;

%%
hrvFeature.hr = 60./rrInterval; %BPM
hrvFeature.rrInterval = rrInterval.*1000; %ms
hrvFeature.tRR = tRR; %s
hrvFeature.meanRR = mean(rrInterval).*1000;%ms
hrvFeature.sdnn = std(rrInterval)*1000; % ms: See matlab definition (N-1)
hrvFeature.meanHR = mean(hrvFeature.hr); % BPM
hrvFeature.sdRR = sdRR .* 1000; %ms
hrvFeature.rmssd = sqrt(mean((sdRR).^2))*1000; % ms 
hrvFeature.sdsdRR = std(sdRR)*1000;  %ms - Standard deviation of successive differences between NN intervals
hrvFeature.pNN50 = 100 * nn50/(length(sdRR)+1); % Percentage pairs of adjacent intervals differing by more than 50ms

hrvFeature.rrIntervalInterp = rrIntervalInterp;
hrvFeature.rrResampleFs = rrResampleFs;
hrvFeature.tRRInterp = tRRInterp;
hrvFeature.fftRR = yRR;
hrvFeature.fftFreqXaxis = f;
hrvFeature.fftPowYaxis = P(1:n/2+1);
hrvFeature.nNextPow2 = n;
hrvFeature.LFpow = LFpow;
hrvFeature.HFpow = HFpow;
hrvFeature.LFHFratio = LFHFratio;
end

function [hrvFeature,pk] = hrvFeatureEstEcg(dataEcg,fs,opts)
if ~isfield(opts,'delPkCalibTime')
    % Make sure this time interval has no wrong beat
    opts.delPkCalibTime = [0 10];
    fprintf('Heartbeat peak amplitude calibration window: [%2.1f, %2.1f]s. \n',opts.delPkCalibTime);
end

if ~isfield(opts,'pkAmpRelRejThresh')
    % This is minimum peak value of a beat relative to mean peak amplitude,
    % identified from delPkCalibTime, below which it will be rejected.
    opts.pkAmpRelRejThresh = 0.6; 
    fprintf('Heartbeat peak amplitude relative rejection threshold: %f\n',opts.pkAmpRelRejThresh);
end
if ~isfield(opts,'tRRthresh')
    % Minimum and maximum time between two peaks
    % This could be made more statistical, for now it is fixed
    opts.tRRthresh = [0.4,1.2]; % 60./[0.4, 1.5] = [150, 40] BPM
    fprintf('Min-Max RR time: [%2.1f, %2.1f]s\n',opts.tRRthresh);
end
if ~isfield(opts,'minEcgPkHt')
    opts.minEcgPkHt = 0.5;
    fprintf('Default ECG min peak height %3.2f\n',opts.minEcgPkHt);
end
if ~isfield(opts,'tWinHR')
    opts.tWinHR = 5;
    fprintf('Default HR estimation window: %3.2f\n',opts.tWinHR);
end
if ~isfield(opts,'minEcgPkDist')
    opts.minEcgPkDist = 0.3;
    fprintf('Default ECG min peak dist %3.2f\n',opts.minEcgPkDist);
end

opts.tRRthresh = [0.4,1.2]; % 60./[0.4, 1.5] = [150, 40] BPM
fprintf('Min-Max RR time: [%2.1f, %2.1f]s\n',opts.tRRthresh);

t = (0:(length(dataEcg)-1))/fs;
t = t(:);

% Simple findpeaks for peak detection
[~,pk] = findpeaks(dataEcg,'MinPeakHeight',opts.minEcgPkHt,'MinPeakDistance',opts.minEcgPkDist.*fs);

figure
hold on;
plot(t,dataEcg);
plot(t(pk),dataEcg(pk),'^');

t = (0:(length(dataEcg)-1))/fs; 
t = t(:);
tPkCalibWindow = zeros(1,2);
tPkCalibWindow(1) = t(opts.delPkCalibTime(1) * fs + 1);
tPkCalibWindow(2) = t(opts.delPkCalibTime(2) * fs + 1);

% Feature 1: Mean RR
% This is time difference between two maximas. This will need correction as
% we may have missed some peaks, and in above peak correction we might have
% ignored some true peaks as well.
rrInterval = t(pk(2:end)) - t(pk(1:end-1));
% tRR: Associating RR interval time at the second peak
tRR = t(pk(2:end)); 

% Correcting: RR interval limited by estimated HR in healthy seater person
idxRRinterval =  (rrInterval <= opts.tRRthresh(2)) & (rrInterval >= opts.tRRthresh(1));
rrInterval = rrInterval(idxRRinterval);
tRR = tRR(idxRRinterval);

% Since other correction uses next rrInterval that could be incorrect,
% editing that
nWinRRcorrect = 10; % Use last 10 RR intervals mean to correct

for i = 1:nWinRRcorrect
    if (rrInterval(i) < (1/1.3)*mean(rrInterval(i+1:nWinRRcorrect))) || (rrInterval(i) > 1.3*mean(rrInterval(i+1:nWinRRcorrect)))
        rrInterval(i) = mean(rrInterval(i+1:nWinRRcorrect));
    end
end

for i = nWinRRcorrect+1:length(rrInterval) - 1
    if (rrInterval(i) > 1.3* mean(rrInterval(i-(nWinRRcorrect):i-1))) || (rrInterval(i) < (1/1.3)* mean(rrInterval(i-(nWinRRcorrect):i-1)))
        rrInterval(i) = rrInterval(i-1); % Use the last interval
    end
end

% Correction: RR interval is not expected to change suddenly with time
% (RR(k) > 0.7?(RR(k ? 1)+RR(k + 1))) -> RR(k) = (RR(k ? 1)+RR(k + 1))/2
% (RR(k) < 0.35*(RR(k ? 1)+RR(k + 1))) -> RR(k) = (RR(k ? 1)+RR(k + 1))/2
for i = 2:length(rrInterval)-1
    if (rrInterval(i) > (0.7*(rrInterval(i-1) + rrInterval(i+1))) || rrInterval(i) < (0.35*(rrInterval(i-1) + rrInterval(i+1))))
        rrInterval(i) = (rrInterval(i-1) + rrInterval(i+1))/2;
    end
end

fprintf('ECG: Mean RR interval: %3.2f s\n',mean(rrInterval));

%% Succesive RR information: In some cases we miss some beats, so the 
% rrInterval is not successive. That information is important for finding
% successive differences - and affects RMSSD, pNN50, etc.
tRRdiff = tRR(2:end) - tRR(1:end-1);
isNextSuccessive = tRRdiff <= opts.tRRthresh(2); % Is next RR interval successive or not

% Saving only successive differences, otherwise 0
sdRR = isNextSuccessive .* (rrInterval(2:end) - rrInterval(1:end-1));
sdRR = sdRR(sdRR~=0);
nn50 = sum(sdRR > 50e-3); % Number of successive differences greater than 50ms

%% Interpolated RR or Inter-Beat Interval (IBI)
% Selecting sampling frequency of 10 Hz.
rrResampleFs = 10;
tRRInterp = tRR(1):1/rrResampleFs:tRR(end);
rrIntervalInterp = interp1(tRR,rrInterval,tRRInterp);

figure
plot(tRR,rrInterval,tRRInterp,rrIntervalInterp);
title('ECG Interpolated RR interval');

%% Frequency analysis: PSD/STFT/wavelet

% Starting with FFT to compare entire baseline vs test durations
n = 2^nextpow2(length(rrIntervalInterp));
yRR = fft(rrIntervalInterp - mean(rrIntervalInterp),n);
f = rrResampleFs*(0:(n/2))/n;
P = abs(yRR).^2/n;

figure
plot(f,P(1:n/2+1)) 
title('Frequency spectrum ECG RR interval')

% LF power: 0.04 - 0.15 Hz, HF: 0.15 - 0.7 Hz
idxLF1 = find(f<=0.04,1,'last');
idxLF2HF1 = find(f>=0.15,1,'first');
idxHF2 = find(f>=0.7,1,'first');

LFpow = trapz(P(idxLF1:idxLF2HF1)); % Not normalized, and may depend on time period involved (? some paper)
HFpow = trapz(P(idxLF2HF1:idxHF2)); % Not normalized
LFHFratio = LFpow/HFpow;

%%
hrvFeature.hr = 60./rrInterval; %BPM
hrvFeature.rrInterval = rrInterval.*1000; %ms
hrvFeature.tRR = tRR; %s
hrvFeature.meanRR = mean(rrInterval).*1000;%ms
hrvFeature.sdnn = std(rrInterval)*1000; % ms: See matlab definition (N-1)
hrvFeature.meanHR = mean(hrvFeature.hr); % BPM
hrvFeature.sdRR = sdRR .* 1000; %ms
hrvFeature.rmssd = sqrt(mean((sdRR).^2))*1000; % ms 
hrvFeature.sdsdRR = std(sdRR)*1000;  %ms - Standard deviation of successive differences between NN intervals
hrvFeature.pNN50 = 100 * nn50/(length(sdRR)+1); % Percentage pairs of adjacent intervals differing by more than 50ms

hrvFeature.rrIntervalInterp = rrIntervalInterp;
hrvFeature.rrResampleFs = rrResampleFs;
hrvFeature.tRRInterp = tRRInterp;
hrvFeature.fftRR = yRR;
hrvFeature.fftFreqXaxis = f;
hrvFeature.fftPowYaxis = P(1:n/2+1);
hrvFeature.nNextPow2 = n;
hrvFeature.LFpow = LFpow;
hrvFeature.HFpow = HFpow;
hrvFeature.LFHFratio = LFHFratio;
end    
    

