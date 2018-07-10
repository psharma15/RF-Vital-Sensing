%% Time-Frequency Analysis: Motion Detection

%% From postProcess3.m
% Adding path, if needed
dataPath = 'D:\Research\SummerFall17\CnC\NCS\MotionAffectedData\Exp1\Script1';
addpath(dataPath);

% Reading ecg-ncs data.
load('Tag2Amp 2-41-35 PM.txt') % Tag2Amp 2-41-35 PM.txt Tag2Amp 2-50-33 PM.txt

ncsUnfilt = Tag2Amp_2_41_35_PM(1:1.65e5); % Tag2Amp_2_41_35_PM Tag2Amp_2_50_33_PM
ncsUnfilt = - ncsUnfilt;

clearvars Tag1Amp_2_41_35_PM

fs = 500;
nSample = length(ncsUnfilt);
idx = 1:nSample;
tSample = ((idx-1)/fs)';

rmpath(dataPath);

% Viewing data and frequency content
df = fs/nSample;
f = -fs/2:df:fs/2 - df;
XncsUnfilt = fftshift(fft(ncsUnfilt));

figure
ax1(1) = subplot(2,2,1);
plot(tSample,ncsUnfilt);
title('Unfiltered NCS');
xlabel('Time(sec)');
ylabel('Measured Amplitude');
grid on

ax1(2) = subplot(2,2,3);
plot(f,abs(XncsUnfilt)/nSample);
title('Frequency Spectrum of unfiltered NCS');

% Filtered ncs: Only dc and very high frequency removal
filtNcs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',0.05,0.6,4,5,150,0.1,150,fs);
% M = designmethods(filtNcs);
Hd = design(filtNcs,'kaiserwin'); %kaiserwin is faster than equiripple
% fvtool(Hd)
ncsFiltered = filtfilt(Hd.Numerator,1,ncsUnfilt); % Taking care of group delay

% At least first 30 sec is considered to be at rest. Normalized according
% to that.
nSampPeakDet = find(tSample>=30,1);
ncsPks = findpeaks(ncsFiltered(1:nSampPeakDet),'MinPeakHeight',0.0001,'MinPeakDistance',300); %minpeakheight 1e-3
maxNcs = mean(ncsPks);
ncs = ncsFiltered./(0.8*maxNcs);

ax1(3) = subplot(2,2,2);
plot(tSample,ncs)
xlabel('Time(sec)')
ylabel('Normalized Amplitude')
title('Filtered NCS')
grid on

Xncs = fftshift(fft(ncs));
ax1(4) = subplot(2,2,4);
plot(f,abs(Xncs)/nSample);
title('Frequency Spectrum of filtered NCS');

linkaxes([ax1(1),ax1(3)],'x')
linkaxes([ax1(2),ax1(4)],'x')

%% 
figure
windowLen = 5*fs; % 10 sec is window length, Hamming by default
spectrogram(ncs,windowLen,[],[],fs,'yaxis')
ylim([0.1, 5])