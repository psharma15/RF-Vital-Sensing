
function [ncsAmpFilt,ncsPhFilt,ncsUnfiltAmpTrunc,ncsUnfiltPhTrunc,tTrunc] = postProcess2(f3db,fp,fst,dataPath,fileName,signAmp,varargin)
%POSTPROCESS2 function is for post processing NCS signal. It applies
%high pass, low pass filter to get dc-removed, respiratory, or heartbeat 
%signal.
% Example: [ncsAmpFilt,~,~,~,t] = postProcess2(0.9,10,11,dataPath,'b1_1.mat',-1,0,30);
% where, f3db = 0.1, 3dB cutoff frequency of butterworth iir hp filter
% fp, fst = passband and stopband frequency of kaiserwin fir lp filter
% dataPath = Path of folder with data
% fileName = '.mat' data file
% signAmp = This multiplies Amplitude data with +/-1 - take care of
%           polarity reversal.
% varargin{1} = sampStart: Start time(in seconds) of truncated data
% varargin{2} = sampEnd: End time (in seconds) of truncated data

%% Load data
addpath(dataPath);
load(fileName);

%%
fs = 500; % Sampling frequency - constant
ncsUnfiltAmp = VarName1;
nSample = length(ncsUnfiltAmp);
idx = 1:nSample;
t = ((idx-1)/fs)';

%%
sampStart = 0*fs + 1;
sampEnd = nSample;

if nargin > 6
    sampStart = varargin{1}*fs+1;
    if nargin == 8
        sampEnd = varargin{2}*fs+1;
    end
end

%%
ncsUnfiltAmpTrunc = signAmp*VarName1(sampStart:sampEnd);
ncsUnfiltPhTrunc = VarName2(sampStart:sampEnd);
tTrunc = t(sampStart:sampEnd);

% High Pass filter - for dc or dc and breath removal
orderHP = 20;
filtHP = fdesign.highpass('N,F3db',orderHP,f3db,fs); % Nonlinear phase filter but better filtering characteristics
HdHP = design(filtHP,'butter'); % 'butter' with 'N,F3db' specifications
ncsAmpHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltAmpTrunc);
ncsPhHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltPhTrunc);

% Low pass filter 
filtLP = fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fst,0.1,80,fs); % Change Ast to 10 if filtfilt creates error
HdLP = design(filtLP,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
ncsAmpLP = filtfilt(HdLP.Numerator,1,ncsAmpHP); % Taking care of group delay
ncsPhLP = filtfilt(HdLP.Numerator,1,ncsPhHP); % Taking care of group delay

ncsAmpFilt = ncsAmpLP;
ncsPhFilt = ncsPhLP;

% Frequency spectrum
XncsAmpUnfilt = fftshift(fft(ncsUnfiltAmpTrunc));
XncsPhUnfilt = fftshift(fft(ncsUnfiltPhTrunc));

XncsAmpFilt = fftshift(fft(ncsAmpFilt));
XncsPhFilt = fftshift(fft(ncsPhFilt));
dfTrunc = fs/length(ncsUnfiltAmpTrunc);
fTrunc = -fs/2:dfTrunc:(fs/2) - dfTrunc;

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

rmpath(dataPath);
end
