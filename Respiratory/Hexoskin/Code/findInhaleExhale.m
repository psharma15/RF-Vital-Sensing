% This function finds inhalation (inhalation end) and exhalation 
% (exhalation end) points in respiration waveform extracted from both 
% amplitude and phase of NCS. 
% The peak detection is performed by AMPD (automatic multiscale based peak
% detection). The source is: "An efficient algorithm for automatic peak
% detection in noisy periodic and quasi-periodic signals", F Scholkmann.
% Pragya Sharma, ps847@cornell.edu
% April 15th, 2018.

function [inExAmp,inExPh] = findInhaleExhale(data,fs)
% Input:
%   data: [ncs Amp, ncs Ph]
%   fs: Sampling frequency of data
% Output:
%   inExAmp, inExPh: First coulum gives index in the data, and second gives
%   indicator, 1: Inhalation and 0: exhalation, for amplitude and phase
%   respectively.

%% Finding peaks: both maximum and minimum

% Considering respiration is between 8-60 breaths per minute - normal for
% adult is 8-20 bpm.
freqRange = [4, 60]./60; % In Hz
% peak detection function gives indices corresponding to inhalation and
% exhalation. These indices correspond to location in data, thereby,
% retaining the time information for each inhalation and exhalation.
[locInhaleAmp,locExhaleAmp] = peakDet(data(:,1),freqRange,fs);
[locInhalePh, locExhalePh]  = peakDet(data(:,2),freqRange,fs);

% Somehow repeated inhale exhale points - most likely due to uncorrected
% overlap segments. Can be corrected if nearby. But some issue if exact
% point repeated. So checking that
for i = 2:length(locInhaleAmp)
    if locInhaleAmp(i) == locInhaleAmp(i-1)
        locInhaleAmp(i) = 0;
    end
end
locInhaleAmp = locInhaleAmp(locInhaleAmp ~= 0);

for i = 2:length(locExhaleAmp)
    if locExhaleAmp(i) == locExhaleAmp(i-1)
        locExhaleAmp(i) = 0;
    end
end
locExhaleAmp = locExhaleAmp(locExhaleAmp ~= 0);

for i = 2:length(locInhalePh)
    if locInhalePh(i) == locInhalePh(i-1)
        locInhalePh(i) = 0;
    end
end
locInhalePh = locInhalePh(locInhalePh ~= 0);

for i = 2:length(locExhalePh)
    if locExhalePh(i) == locExhalePh(i-1)
        locExhalePh(i) = 0;
    end
end
locExhalePh = locExhalePh(locExhalePh ~= 0);

%% These are obtained maximas and minimas
t = 0:(1/fs):((length(data(:,1))-1)/fs);

figure
ax(1) = subplot(2,1,1);
plot(t,data(:,1)); xlabel('Time (sec)'); ylabel('NCS Amp')
hold on
plot(t(locInhaleAmp),data(locInhaleAmp,1),'*');
plot(t(locExhaleAmp),data(locExhaleAmp,1),'o');
grid on
title('Inhalation and exhalation points from peak detection algorithm')

ax(2) = subplot(2,1,2);
plot(t,data(:,2)); xlabel('Time (sec)'); ylabel('NCS Ph')
hold on 
plot(t(locInhalePh),data(locInhalePh,2),'*');
plot(t(locExhalePh),data(locExhalePh,2),'o');

grid on
linkaxes(ax,'x')

%% Need to correct Inhale and exhale points

% First combining Inhale and Exhale points 
% Denoting 1: Inhalation, 0: Exhalation
inExAmpIndicator = [ones(length(locInhaleAmp),1); ...
    zeros(length(locExhaleAmp),1)]; % First for amplitude

[inExAmpIdx, sortIdx] = sort([locInhaleAmp(:); locExhaleAmp(:)]);
inExAmpIndicator = inExAmpIndicator(sortIdx);

inExPhIndicator = [ones(length(locInhalePh),1); ...
    zeros(length(locExhalePh),1)]; % For phase
[inExPhIdx, sortIdx] = sort([locInhalePh(:); locExhalePh(:)]);
inExPhIndicator = inExPhIndicator(sortIdx);

% Correcting consecutive maxima or minima. First column gives index, second
% is corresponding indicator.
inExAmp = zeros(length(inExAmpIdx),2); 
inExAmp(1,:) = [inExAmpIdx(1), inExAmpIndicator(1)];
counter = 2;
tempMaxMin = [0,0,0]; % saves [inExAmpIdx, inExAmpIndicator, counter]
dataMaxMin = 1e10;

for iter = 2:length(inExAmpIdx)
    if inExAmpIndicator(iter) ~= inExAmpIndicator(iter-1)
        if tempMaxMin(1) ~= 0
            % Updating max/ min if there were any multiple consecutive
            % max/ min, with a different peak than the one saved for that
            % counter.
            inExAmp(tempMaxMin(3),:) = [tempMaxMin(1), tempMaxMin(2)];
            tempMaxMin = [0,0,0];
        end
        inExAmp(counter,:) = [inExAmpIdx(iter), inExAmpIndicator(iter)];
        counter = counter + 1;
        dataMaxMin = data(inExAmpIdx(iter),1); 
    end
    if (inExAmpIndicator(iter) == inExAmpIndicator(iter-1))
        if dataMaxMin == 1e10
            % This is the max/ min data pt among consecutive max/ min
            dataMaxMin = data(inExAmpIdx(iter-1),1); 
        end
        if ((inExAmpIndicator(iter) == 1) && (data(inExAmpIdx(iter),1) > dataMaxMin)) ...
                || ((inExAmpIndicator(iter) == 0) && (-data(inExAmpIdx(iter),1) > -dataMaxMin))
            dataMaxMin = data(inExAmpIdx(iter),1);
            tempMaxMin = [inExAmpIdx(iter), inExAmpIndicator(iter), ...
                counter-1];
        end
    end
end
inExAmp = inExAmp(1:(counter-1),:);
locInhaleAmpCorrected = inExAmp((inExAmp(:,2) == 1),1);
locExhaleAmpCorrected = inExAmp((inExAmp(:,2) == 0),1);

inExPh = zeros(length(inExPhIdx),2);
inExPh(1,:) = [inExPhIdx(1), inExPhIndicator(1)];
counter = 2;
tempMaxMin = [0,0,0]; % saves [inExIdx, inExIndicator, counter]
dataMaxMin = 1e10;

for iter = 2:length(inExPhIdx)
    if inExPhIndicator(iter) ~= inExPhIndicator(iter-1)
        if tempMaxMin(1) ~= 0
            % Updating max/ min if there were any multiple consecutive
            % max/ min, with a different peak than the one saved for that
            % counter.
            inExPh(tempMaxMin(3),:) = [tempMaxMin(1), tempMaxMin(2)];
            tempMaxMin = [0,0,0];
        end
        inExPh(counter,:) = [inExPhIdx(iter), inExPhIndicator(iter)];
        counter = counter + 1;
        dataMaxMin = data(inExPhIdx(iter),2);
    end
    if (inExPhIndicator(iter) == inExPhIndicator(iter-1))
        if dataMaxMin == 1e10
            % This is the max/ min data pt among consecutive max/ min
            dataMaxMin = data(inExPhIdx(iter-1),2); 
        end
        if ((inExPhIndicator(iter) == 1) && (data(inExPhIdx(iter),2) > dataMaxMin)) ...
                || ((inExPhIndicator(iter) == 0) && (-data(inExPhIdx(iter),2) > -dataMaxMin))
            dataMaxMin = data(inExPhIdx(iter),2);
            tempMaxMin = [inExPhIdx(iter), inExPhIndicator(iter), ...
                counter-1];
        end
    end

end
inExPh = inExPh(1:(counter-1),:);
locInhalePhCorrected = inExPh((inExPh(:,2) == 1),1);
locExhalePhCorrected = inExPh((inExPh(:,2) == 0),1);

%% Plotting corrected locations
figure
title('Corrected Inhalation and Exhalation Indices')
ax(1) = subplot(2,1,1);
plot(t,data(:,1)); xlabel('Time (sec)'); ylabel('NCS Amp')
hold on
plot(t(locInhaleAmpCorrected),data(locInhaleAmpCorrected,1),'*');
plot(t(locExhaleAmpCorrected),data(locExhaleAmpCorrected,1),'o');
grid on

ax(2) = subplot(2,1,2);
plot(t,data(:,2)); xlabel('Time (sec)'); ylabel('NCS Ph')
hold on 
plot(t(locInhalePhCorrected),data(locInhalePhCorrected,2),'*');
plot(t(locExhalePhCorrected),data(locExhalePhCorrected,2),'o');

grid on
linkaxes(ax,'x')
