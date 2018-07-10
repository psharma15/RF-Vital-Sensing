
% function [ncs,ncsUnfilt,t] = postProcess1(Fst1,Fp1,Fp2,Fst2,varargin)
%%POSTPROCESS3 function is for post processing NCS signal. It applies
%%bandpass filter in range defined by Fst1,Fp1,Fp2,Fst2. 
% Example: [ncs,~] = postProcess3(0.1,0.6,5,6);
% where, Fst1 = 0.01; Fp1 = 0.6; Fp2 = 5; Fst2 = 6;
% 

% if nargin==4
%     scriptNum = 1;
%     tagNum = 2;
% else
%     scriptNum = varargin{1};
%     tagNum = varargin{2};
% end
% 
% script = {'Script1','Script2','Script3'};
% tagFileName = {'Tag1Amp 2-41-35 PM.txt','Tag2Amp 2-41-35 PM.txt';...
%     'Tag1Amp 2-50-33 PM.txt','Tag2Amp 2-50-33 PM.txt';...
%     'Tag1Amp 3-02-47 PM.txt','Tag2Amp 3-02-47 PM.txt'};

%% Reading ecg-ncs data.
% fileData = load('test.xlsx');
% load 'a9.mat'
% ncsUnfilt = a9;
ncsUnfilt = VarName1;
% ncsUnfilt = -ncsUnfilt;
ncsPh = VarName2;

%%
fs = 500; % This data is sampled at 500 Hz.
nSample = length(ncsUnfilt);
idx = 1:nSample;
t = ((idx-1)/fs)';

%% Plot amp and phase in same figure
figure
yyaxis left
plot(t,ncsUnfilt,'markers',24)
xlabel('Time (sec)','FontSize',16)
ylabel('Amplitude','FontSize',16)

yyaxis right
plot(t,unwrap(ncsPh),'markers',24)
ylabel('Unwrapped Phase','FontSize',16)

% axis('tight'); 
grid on;

%% Viewing data and frequency content
df = fs/nSample;
f = -fs/2:df:fs/2 - df;
XncsUnfilt = fftshift(fft(ncsUnfilt));

nRow = 2;
nCol = 3;

figure
ax(1) = subplot(nRow,nCol,1);
plot(t,ncsUnfilt);
title('Unfiltered NCS');
xlabel('Time(sec)');
ylabel('Measured Amplitude');
grid on

ax(2) = subplot(nRow,nCol,4);
plot(f,abs(XncsUnfilt)/nSample);
title('Frequency Spectrum of unfiltered NCS');

%%
% Fst1 = 0.3; Fp1 = 0.9; Fp2 = 50; Fst2 = 51; % For heart beat
Fst1 = 0.9; Fp1 = 1.2; Fp2 = 50; Fst2 = 51; % For heart beat

% Fp = 10; Fst = 11; % 
% Fst1 = 0.001; Fp1 = 0.3; Fp2 = 5; Fst2 = 6; 
% Fst1 = 0.001; Fp1 = 1; Fp2 = 20; Fst2 = 21; % For eye blinking
% Fst1 = 0.01; Fp1 = 0.6; Fp2 = 20; Fst2 = 21; % 

%% Filtered ncs
filtNcs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,80,0.1,80,fs);
% filtNcs = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,0.1,80,fs);

% M = designmethods(filtNcs);
% Do lineae phase FIR filter. 
Hd = design(filtNcs,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
% fvtool(Hd)
ncsFiltered = filtfilt(Hd.Numerator,1,ncsUnfilt); % Taking care of group delay

% At least first 30 sec is considered to be at rest. Normalized according
% to that.
% nSampPeakDet = find(t>=30,1);
nSampPeakDet = find(t>=30,1);
ncsPks = findpeaks(ncsFiltered(1:nSampPeakDet),'MinPeakHeight',0.0001,'MinPeakDistance',300); %minpeakheight 1e-3
maxNcs = mean(ncsPks);
ncs = ncsFiltered./(0.8*maxNcs);

ax(3) = subplot(nRow,nCol,2);
plot(t,ncs)
xlabel('Time(sec)')
ylabel('Normalized Amplitude')
title('Filtered NCS')
grid on

Xncs = fftshift(fft(ncs));
ax(4) = subplot(nRow,nCol,5);
plot(f,abs(Xncs)/nSample);
title('Frequency Spectrum of filtered NCS');

%%
ax(5) = subplot(nRow,nCol,3);
plot(t,-unwrap(ncsPh))
xlabel('Time(sec)')
ylabel('Phase in degrees')
title('NCS phase')
grid on

linkaxes([ax(1),ax(3),ax(5)],'x')
linkaxes([ax(2),ax(4)],'x')

axis([ax(1),ax(3),ax(5)],'tight')
% end
