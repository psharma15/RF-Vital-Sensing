% Code snippets for processing

clearvars -except VarName1 VarName2

%% 
% 1. DC Removal
% 2. LP filtering to get breath
% 3. Original DC removed signal - Breath signal (1-2)
% Result: Performs very similar to band pass filtering (as long as order of
% filters are low, otherwise calculation errors, and of course this signal
% has more HF content than band pass filtered, as no LP filter is applied
% at the last stage. 

fs = 500;
ncsUnfiltAmp = VarName1;
nSample = length(ncsUnfiltAmp);
idx = 1:nSample;
t = ((idx-1)/fs)';
df = fs/nSample;
f = -fs/2:df:fs/2 - df;

% Remove dc - low order filter

orderHP = 2;
f3db = 0.08;
filtHP = fdesign.highpass('N,F3db',orderHP,f3db,fs);
hpFilterOptions = designmethods(filtHP);
Hd = design(filtHP,'butter');
ncsHP = filtfilt(Hd.sosMatrix,Hd.ScaleValues,ncsUnfiltAmp);
XncsHP = fftshift(fft(ncsHP));

nRow = 3;
nCol = 2;

figure
ax1(1) = subplot(nRow,nCol,1);
plot(t,ncsHP); grid on;

ax1(2) = subplot(nRow,nCol,2);
plot(f,abs(XncsHP)/nSample); grid on;

% Low pass filtering to extract breath

Fp = 0.7; Fst = 1; 
% Fp = 0.9; Fst = 1.2; 

filtNcs = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,0.1,80,fs);

Hd = design(filtNcs,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
ncsLP = filtfilt(Hd.Numerator,1,ncsHP); % Taking care of group delay
XncsLP = fftshift(fft(ncsLP));

ax1(3) = subplot(nRow,nCol,3);
plot(t,ncsLP); grid on;

ax1(4) = subplot(nRow,nCol,4);
plot(f,abs(XncsLP)/nSample); grid on;

% 
ncsHF = ncsHP-ncsLP;
XncsHF = fftshift(fft(ncsHF));

ax1(5) = subplot(nRow,nCol,5);
plot(t,ncsHF); grid on;

ax1(6) = subplot(nRow,nCol,6);
plot(f,abs(XncsHF)/nSample); grid on;

linkaxes([ax1(1),ax1(3),ax1(5)],'x')
linkaxes([ax1(2),ax1(4),ax1(6)],'x')

%
postProcess1

figure
plot(t,ncsHP,':',t,ncsLP,'-')
hold on
plot(t,ncsHF,'r')
plot(t,ncsFiltered,'g')
grid on

legend('High Pass filtered - DC removal','Low Pass filtered','High - Low','Band Pass filtered')

%% Analyze small section

sampStart = 15*fs;
sampEnd = 20*fs-1;
nSampleTruncated = sampEnd - sampStart + 1;

XTruncated = fftshift(fft(ncsHF(sampStart:sampEnd)));
dfTruncated = fs/nSampleTruncated;
fTruncated = -fs/2:dfTruncated:(fs/2)-dfTruncated;
figure
plot(fTruncated,abs(XTruncated)/nSampleTruncated);

%% Band stop - notch filter
sampStart = 0*fs + 1;
sampEnd = 60*fs;

d = designfilt('bandstopiir','FilterOrder',10,'HalfPowerFrequency1',2.2,...
    'HalfPowerFrequency2',3.6,'DesignMethod','butter','SampleRate',fs);
% fvtool(d,'Fs',fs)
ncsNotch = filtfilt(d,ncsHF(sampStart:sampEnd));

d = designfilt('bandstopiir','FilterOrder',10,'HalfPowerFrequency1',1.4,...
    'HalfPowerFrequency2',1.8,'DesignMethod','butter','SampleRate',fs);
% fvtool(d,'Fs',fs)
ncsNotch2 = filtfilt(d,ncsNotch);

hold on;
plot(t(sampStart:sampEnd),ncsNotch,'k')
hold off

XPk = fftshift(fft(ncsNotch2));
dfPk = fs/length(ncsNotch2);
fPk = -fs/2:dfPk:(fs/2) - dfPk;
figure
plot(fPk,abs(XPk)/length(ncsNotch2));

%% Band pass - peak filter
sampStart = 0*fs + 1;
sampEnd = 60*fs;

peakspec = fdesign.peak('N,F0,BW',2,1,0.6,fs);
peakfilt1 = design(peakspec,'SystemObject',true);
ncsPk1 = filtfilt(peakfilt1.SOSMatrix,peakfilt1.ScaleValues,ncsHF(sampStart:sampEnd));

peakspec = fdesign.peak('N,F0,BW',2,1.3,0.6,fs);
peakfilt2 = design(peakspec,'SystemObject',true);
ncsPk2 = filtfilt(peakfilt2.SOSMatrix,peakfilt2.ScaleValues,ncsHF(sampStart:sampEnd));

peakspec = fdesign.peak('N,F0,BW',2,1.5,0.4,fs);
peakfilt3 = design(peakspec,'SystemObject',true);
ncsPk3 = filtfilt(peakfilt3.SOSMatrix,peakfilt3.ScaleValues,ncsHF(sampStart:sampEnd));

peakspec = fdesign.peak('N,F0,BW',2,2.1,0.5,fs);
peakfilt4 = design(peakspec,'SystemObject',true);
ncsPk4 = filtfilt(peakfilt4.SOSMatrix,peakfilt4.ScaleValues,ncsHF(sampStart:sampEnd));

peakspec = fdesign.peak('N,F0,BW',2,2.7,0.3,fs);
peakfilt5 = design(peakspec,'SystemObject',true);
ncsPk5 = filtfilt(peakfilt5.SOSMatrix,peakfilt5.ScaleValues,ncsHF(sampStart:sampEnd));

ncsPk = 0.5*(3*ncsPk1 + 1*ncsPk2 + 0.9*ncsPk3 + 0.1*ncsPk4 + 0.1*ncsPk5);

hold on;
plot(t(sampStart:sampEnd),ncsPk,'k')
hold off

XPk = fftshift(fft(ncsHF(sampStart:sampEnd)));
% XPk = fftshift(fft(ncsPk));
dfPk = fs/length(ncsPk);
fPk = -fs/2:dfPk:(fs/2) - dfPk;
figure
plot(fPk,abs(XPk)/length(ncsPk));

%% STFT View of Motion Affected Data
figure
% spectrogram(ncsHF,2*fs,0.5*fs,[],fs,'MinThreshold',-130,'yaxis')
% spectrogram(ncsHF,2*fs,0.5*fs,[],fs,'reassign','yaxis')
spectrogram(ncsHF,hamming(2*fs),0.5*fs,[],fs,'reassign','MinThreshold',-130,'yaxis')

% Take window and find fft
% sampStart = 1;
% windowSize = 2*fs;
% nbins = floor(nSample/windowSize) + 1;
% for i = 1:nbins
%     sampEnd = sampStart + windowSize - 1;
%     Xbin = fftshift(fft(ncsHF(sampStart:sampEnd)));
%     dfbin = fs/windowSize;
%     fbin = -fs/2

%% DC Filtered Amplitude and phase weightings (when no breath)
fs = 500;
ncsUnfiltAmp = VarName1;
nSample = length(ncsUnfiltAmp);
idx = 1:nSample;
t = ((idx-1)/fs)';
df = fs/nSample;
f = -fs/2:df:fs/2 - df;

sampStart = 0*fs + 1; %220, 96
% sampEnd = nSample;
sampEnd = 100*fs; % ,107

ncsUnfiltAmpTrunc = -VarName1(sampStart:sampEnd);
ncsUnfiltPhTrunc = VarName2(sampStart:sampEnd);
tTrunc = t(sampStart:sampEnd);

% High Pass filter - for dc or dc and breath removal
orderHP = 20;
f3db = 0.9; % 0.08 for DC filtering only, ~0.8 for DC+breath filtering
filtHP = fdesign.highpass('N,F3db',orderHP,f3db,fs); % Nonlinear phase filter but better filtering characteristics
hpFilterOptions = designmethods(filtHP);
HdHP = design(filtHP,'butter'); % 'butter' with 'N,F3db' specifications
ncsAmpHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltAmpTrunc);
ncsPhHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltPhTrunc);

% HdHP = designfilt('highpassfir','StopbandFrequency',0.001,...
%     'PassbandFrequency',0.05,'PassbandRipple',0.1,'StopbandAttenuation',40,...
%     'DesignMethod','kaiserwin','SampleRate',fs);
% ncsAmpHP = filtfilt(HdHP,ncsUnfiltAmpTrunc);
% ncsPhHP = filtfilt(HdHP,ncsUnfiltPhTrunc);
% fvtool(HdHP,'Fs',fs);

% Low pass filter (No notch filter needed)
Fp = 10; Fst = 11; 
filtLP = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,0.1,80,fs); % Change Ast to 10 if filtfilt creates error
HdLP = design(filtLP,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
% fvtool(HdLP,'Fs',fs);
ncsAmpLP = filtfilt(HdLP.Numerator,1,ncsAmpHP); % Taking care of group delay
ncsPhLP = filtfilt(HdLP.Numerator,1,ncsPhHP); % Taking care of group delay

% orderHP = 20;
% f3db = 10; % 0.08 for DC filtering only, ~0.8 for DC+breath filtering
% filtLP = fdesign.lowpass('N,F3db',orderHP,f3db,fs); % Nonlinear phase filter but better filtering characteristics
% hpFilterOptions = designmethods(filtLP);
% HdLP = design(filtLP,'butter'); % 'butter' with 'N,F3db' specifications
% ncsAmpLP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsAmpHP);
% ncsPhLP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsPhHP);

% % Notch filter @ 60 Hz, Not linear phase
% dNotch = designfilt('bandstopiir','FilterOrder',10,'HalfPowerFrequency1',59,...
%     'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',fs);
% % figure; fvtool(dNotch,'Fs',fs)
% ncsAmpNotch = filtfilt(d,ncsAmpHP);
% ncsPhNotch = filtfilt(d,ncsPhHP);

ncsAmpFilt = ncsAmpLP;
ncsPhFilt = ncsPhLP;

% Frequency spectrum
XncsAmpUnfilt = fftshift(fft(ncsUnfiltAmpTrunc));
XncsPhUnfilt = fftshift(fft(ncsUnfiltPhTrunc));

XncsAmpFilt = fftshift(fft(ncsAmpFilt));
XncsPhFilt = fftshift(fft(ncsPhFilt));
dfTrunc = fs/length(ncsUnfiltAmpTrunc);
fTrunc = -fs/2:dfTrunc:(fs/2) - dfTrunc;

% Assuming heartbeat is the only signal, ratio of weightings
ratioAmpPh = ncsAmpFilt./ncsPhFilt;

nFig = 2;
figure
fontSize = 12;
ax(1) = subplot(nFig,1,1);
yyaxis left
plot(tTrunc,ncsUnfiltAmpTrunc,'markers',24)
xlabel('Time (sec)','FontSize',fontSize)
ylabel('Amplitude','FontSize',fontSize)
yyaxis right
plot(tTrunc,unwrap(ncsUnfiltPhTrunc),'markers',24)
ylabel('Unwrapped Phase','FontSize',fontSize)
grid on;
ax(2) = subplot(nFig,1,2);
yyaxis left
plot(tTrunc,ncsAmpFilt,'markers',24)
xlabel('Time (sec)','FontSize',fontSize)
ylabel('Filtered Amplitude','FontSize',fontSize)
yyaxis right
plot(tTrunc,unwrap(ncsPhFilt),'markers',24)
ylabel('Filtered Phase','FontSize',fontSize)
grid on;
% ax(3) = subplot(nFig,1,3);
% plot(tTrunc,ratioAmpPh,'markers',24)
% xlabel('Time (sec)','FontSize',fontSize)
% ylabel('Ratio of amplitude to phase','FontSize',fontSize)
linkaxes(ax,'x')

figure
ax1(1) = subplot(2,2,1);
plot(fTrunc,abs(XncsAmpUnfilt)/length(ncsAmpHP));
xlabel('Freq (Hz)','FontSize',fontSize)
ylabel('Amp Spectrum','FontSize',fontSize)
ax1(2) = subplot(2,2,2);
plot(fTrunc,abs(XncsPhUnfilt)/length(ncsPhHP));
xlabel('Freq (Hz)','FontSize',fontSize)
ylabel('Ph Spectrum','FontSize',fontSize)
ax1(3) = subplot(2,2,3);
plot(fTrunc,abs(XncsAmpFilt)/length(ncsAmpHP));
xlabel('Freq (Hz)','FontSize',fontSize)
ylabel('Filtered Amp Spectrum','FontSize',fontSize)
ax1(4) = subplot(2,2,4);
plot(fTrunc,abs(XncsPhFilt)/length(ncsPhHP));
xlabel('Freq (Hz)','FontSize',fontSize)
ylabel('Filtered Ph Spectrum','FontSize',fontSize)
linkaxes(ax1,'x')

% Viewing overlapped frequency spectrum
figure
yyaxis left
plot(fTrunc,abs(XncsAmpFilt)/length(ncsAmpHP));
xlabel('Freq (Hz)','FontSize',fontSize)
ylabel('Filtered Amp Spectrum','FontSize',fontSize)
yyaxis right
plot(fTrunc,abs(XncsPhFilt)/length(ncsPhHP));
ylabel('Filtered Ph Spectrum','FontSize',fontSize)
xlim([-10,10])

% Correlation of amplitude and phase measurement (on filtered signal)
[corrAmpPh,lags] = xcorr(ncsAmpFilt,ncsPhFilt);
[~,I] = max(abs(corrAmpPh));
lagDiff = lags(I);
timeDiff = lagDiff/fs; % Time shift between amplitude and phase waveforms.

figure
plot(lags/fs,corrAmpPh,'markers',20)
hold on
plot(timeDiff,max(abs(corrAmpPh)),'.','MarkerSize',30)
title(['Correlation vs lag: Maximum at lag, ',num2str(timeDiff),' sec'],'FontSize',14)
xlabel('Time lag (sec)','FontSize',14)
ylabel('Correlation','FontSize',14)
grid on; axis 'tight'
% ylim([-3.2,4.5])

%% This follows previos section, when breath is not filtered
% Low pass breath signal from amp/ ph 
Fp = 0.6; Fst = 0.9; 
filtNcs = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,0.5,10,fs);
Hd = design(filtNcs,'kaiserwin');
ncsBreath = filtfilt(Hd.Numerator,1,ncsAmpHP); 
XncsBreath = fftshift(fft(ncsBreath));

figure
yyaxis left
plot(tTrunc,ncsBreath,'markers',24)
xlabel('Time (sec)','FontSize',fontSize)
ylabel('Breath from Amplitude','FontSize',fontSize)
yyaxis right
plot(tTrunc,unwrap(ncsPhFilt),'markers',24)
ylabel('Filtered Phase','FontSize',fontSize)
grid on;

%% This part follows previos section, requires information to calculate ratioFactor.
ratioFactor = 14/2.5e-3;
ncsBreathNorm = ratioFactor*ncsBreath;

figure
plot(tTrunc,ncsBreathNorm,tTrunc,ncsPhFilt)

ncsHeart = ncsPhFilt - ncsBreathNorm;
figure
plot(tTrunc,ncsHeart)

%%