% This function finds inhalation and exhalation points in respiration
% waveform extracted from both amplitude and phase of NCS

function [inhalePts,exhalePts] = findInhaleExhale(data,fs)
%% Input/ Output:
% data: [ncs Amp, ncs Ph]
% fs: Sampling frequency of data
% inhalePts, exhalePts: Inhalation and exhalation location in amplitude and
%                       phase saved as structures and can access .Amp, .Ph

%% Finding peaks: both maximum and minimum

% Considering respiration is between 10-30 breaths per minute - normal for
% adult is 12-20 bpm.
minPeakDistance = (60/30)*fs; % 30 is max breath per min rate considered

minPeakProminenceAmp = 4e-4; % This needs to be intuitive.
minPeakProminencePh = 5; 

[~,locInhaleAmp] = findpeaks( data(:,1),'MinPeakDistance',minPeakDistance,...
                   'MinPeakProminence',minPeakProminenceAmp);
[~,locExhaleAmp] = findpeaks(-data(:,1),'MinPeakDistance',minPeakDistance,...
                   'MinPeakProminence',minPeakProminenceAmp);
               
[~,locInhalePh] =  findpeaks( data(:,2),'MinPeakDistance',minPeakDistance,...
                   'MinPeakProminence',minPeakProminencePh);
[~,locExhalePh] =  findpeaks(-data(:,2),'MinPeakDistance',minPeakDistance,...
                   'MinPeakProminence',minPeakProminencePh);

%% Output data
inhalePts.Amp = locInhaleAmp;
inhalePts.Ph = locInhalePh;
exhalePts.Amp = locExhaleAmp;
exhalePts.Ph = locExhalePh;

%%
t = 0:(1/fs):((length(data(:,1))-1)/fs);

figure
ax(1) = subplot(2,1,1);
plot(t,data(:,1)); xlabel('Time (sec)'); ylabel('NCS Amp')
hold on
plot(t(locInhaleAmp),data(locInhaleAmp,1),'*');
plot(t(locExhaleAmp),data(locExhaleAmp,1),'o');
grid on

ax(2) = subplot(2,1,2);
plot(t,data(:,2)); xlabel('Time (sec)'); ylabel('NCS Ph')
hold on 
plot(t(locInhalePh),data(locInhalePh,2),'*');
plot(t(locExhalePh),data(locExhalePh,2),'*');

grid on
linkaxes(ax,'x')

