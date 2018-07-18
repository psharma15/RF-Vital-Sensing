% This program compares respiratory vital sign and its characteristics as
% obtained from Near-Field Coherent Sensing, against Hexoskin smart garment
% data. 
% Terminology used:
% Inspiration: Peak at end of inspiration.
% Expiration: Minima at end of expiration.
% Tidal Volume (TV): Inspiration and expiration amplitude difference for
%                    the last breath.
% Breath Rate (BR): Breath per minute calculated using last 7 breaths.
% TiTt: Fractional 
% Created on: April 24, 2018
% Pragya Sharma, ps847@cornell.edu
% Edited, July 17, 2018: NCS sampling rate, and ncsFileName.

%% ------------------------------------------------------------------------
% Provide input to function.
dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
    'BreathPattern_Cough_Speak\Data\Pragya\Jul16'];
hxFolder = '\user_13412';
% hxDataNum: 1 = resp_abd, 2 = resp_thrx, 3 = ecg, 4 = tidal_vol_raw,
% 5 = tidal_vol_adj, 6 = min_vent_raw, 7 = min_vent_adj
hxDataNumRespAbd = 1; % Reading abdomen respiration
hxDataNumRespTh = 2; % Reading thorax respiration
hxDataNumTV = 4; % Tidal Volume (TV)
hxDataNumBR = 6; % Breath Rate (BR)
hxDataNumInspExp = [9, 10]; % [Beginning of insp, beginning of exp]
hxDataNumEcg = 3; % ECG_I data
hxDataNumHR = 11; % Heart rate
hxDataNumRR = 12; % RR interval 

ncsFileName = '\freq2G_2v2h'; % data at different time instants
ncsFlipData = [1,1]; % -1 to flip
ncsSampRate = 500; % NCS sampling rate in Hz

% Manual time offset is by observation. Use following settings:
% Yuna_Jul12: 18.4
% 
manualTimeOffset = -3.46; % sec: This is by observation 
ncsTstart = 30; % Time is relative to NCS in seconds, keep 150 for stable
dataDuration = 0; % Leave last 30 sec abd data: 10*60-30-ncsTstart

%% ------------------------------------------------------------------------
% Reading NCS-synchronized Hx waveforms
% Reading ECG
[ncsSync,~,hxEcg,tAbsHxEcg,hxEcgSampRate] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumEcg,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart);

% Reading thoracic and abdomen respiratory waveform
[~,~,hxRespTh,tAbsHxResp,hxRespSampRate] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumRespTh,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart);
[~,~,hxRespAbd,~,~] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumRespAbd,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart);

% Reading Hx tidal volume (TV) and breath rate (BR) 
[~,~,hxTV,tAbsHxTV,hxSampRateTV] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumTV,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart);
[~,~,hxBR,tAbsHxBR,hxSampRateBR] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumBR,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart);

% Reading Hx Inspiration and Expiration events
[~,~,hxExp,tAbsHxExp,~] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumInspExp(2),ncsFileName,...
    ncsSampRate,manualTimeOffset,dataDuration,ncsTstart);
[~,~,hxInsp,tAbsHxInsp,~] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumInspExp(1),ncsFileName,...
    ncsSampRate,manualTimeOffset,dataDuration,ncsTstart);

% Reading Hx Heart Rate
[~,~,hxHR,tAbsHxHR,hxSampRateHR] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumHR,ncsFileName,...
    ncsSampRate,manualTimeOffset,dataDuration,ncsTstart);

% Reading RR interval
[~,~,hxRR,tAbsHxRR,~] = ncsHxRawSync(dataPath,hxFolder,hxDataNumRR,...
    ncsFileName, ncsSampRate, manualTimeOffset,dataDuration,ncsTstart);

% Synchronizing timing of Hx data at different time stamps, taking start
% time of highest sampled respiration (abdomen or thorax) as the reference.
tRef = tAbsHxEcg(1,:);
tHeart = etime(tAbsHxEcg, tRef); % Time starts from 0
tResp = etime(tAbsHxResp,tRef); % Time starts relative to tRef
tTV = etime(tAbsHxTV,tRef); % Time starts relative to tRef
tOffsetTV = tTV(1);
tBR = etime(tAbsHxBR,tRef); % Time starts relative to tRef
tOffsetBR = tBR(1); 
tHxInsp = etime(tAbsHxInsp,tRef); % Time starts relative to tRef
tHxExp = etime(tAbsHxExp,tRef); % Time starts relative to tRef
tHR = etime(tAbsHxHR,tRef);
tHxRR = etime(tAbsHxRR,tRef);

%% ------------------------------------------------------------------------
% Process the NCS data to get respiration waveform - both amplitude 
% and phase, after filtering out heartbeat, and get Heartbeat data after
% filtering out respiration
% *********************************************************************** %
% Specify correct amplitude and phase sign, see if you can implement algo
% to detect this as well. 
% *********************************************************************** %
fprintf('Change amplitude and phase sign if needed...\n');
[ncsRespFiltered,~,~,~] = postProcess(0,1,1.2,ncsSync,ncsSampRate,ncsFlipData);
[ncsHeartFiltered,~,~,~] = postProcess(0.9,10,15,ncsSync,ncsSampRate,ncsFlipData);

%% ------------------------------------------------------------------------
% Downsample NCS resp to hx respiration sample rate = 128 Hz
ncsRespSampRate = hxRespSampRate;
ncsResp = resample(ncsRespFiltered,ncsRespSampRate,ncsHighSampRate);

% Downsample NCS heart to hx ECG sample rate = 256 Hz
ncsHeartSampRate = hxEcgSampRate;
ncsHeart = resample(ncsHeartFiltered,ncsHeartSampRate,ncsHighSampRate);

%% ------------------------------------------------------------------------
% Size check: As time-sync is limited by the min sampling rate between NCS 
% and Hexoskin, sometimes size is one-time point off (As observed after
% re-sampling). Ensuring that size remains same by truncating to the
% minimum size.
if (size(ncsResp,1) - size(hxRespAbd,1)) == 1
    ncsResp = ncsResp(1:length(hxRespAbd),:);
end
if (size(hxRespAbd,1) - size(ncsResp,1)) == 1
    ncsResp = [ncsResp; [ncsResp(end,1), ncsResp(end,2)]];
end
if (size(ncsHeart,1) - size(hxEcg,1)) == 1
    ncsHeart = ncsHeart(1:length(hxEcg),:);
end
if (size(hxEcg,1) - size(ncsHeart,1)) == 1
    ncsHeart = [ncsHeart; [ncsHeart(end,1), ncsHeart(end,2)]];
end
%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS breaths with TV
close all;

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
% tHxResp = 0:1/hx
nFig = 2;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(ax1(1),tResp,hxRespTh./max(hxRespTh),':','color','k','LineWidth',2); %./max(hxRespTh)
hold on
plot(ax1(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]); %./max(hxRespAbd)
plotCute1([],'a.u.',ax1(1),[],[],0);
% ylim([0.9967,1.001])
hold off

yyaxis right
plot(ax1(1),tTV,hxTV);
plotCute1('Time (sec)','mL',ax1(1),...
    'Hexoskin Respiration (thorax & abdomen) and Tidal Volume',{'Thorax (a.u.)',...
    'Abdomen (a.u.)','Tidal Volume (mL)'},1);

ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(ax1(2),tResp,ncsResp(:,1)); 
plotCute1([],'a.u.',ax1(2),[],[],0);

yyaxis right
plot(ax1(2),tResp,ncsResp(:,2)); 
plotCute1('Time (sec)','a.u.',ax1(2),...
    'NCS Respiration (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1);
 
linkaxes(ax1,'x')

%% ------------------------------------------------------------------------
% Find 'inhalation end' and 'exhalation end' indices in amplitude and phase
% inExAmp: Ending of inhalation - peak
% inExPh: Ending of exhalation - minima
% Considering respiration is between 8-60 breaths per minute - normal for
% adult is 8-20 bpm.
freqRangeBR = [4, 60]./60; % In Hz
[inExAmp, inExPh] = findInhaleExhale(ncsResp,ncsRespSampRate,freqRangeBR,tResp);

%% ------------------------------------------------------------------------
% *********************************************************************** %
% Calibrate and Estimate tidal volume with same sampling frequency as Hx.
% Remember to provide correct calibration time. 
% *********************************************************************** %
% calibTime = [tHxTV(1), tHxTV(end)];
calibTime = [17, 22]; % 10-50

[tvCoeffAmpPhSum,tvCoeffAmp,tvCoeffPh,ncsUncalibAmpPhTV] = ...
    ncsEstTV(ncsResp,inExAmp,inExPh,hxTV,ncsRespSampRate,hxSampRateTV,tOffsetTV,calibTime);
ncsTVAmpPhSum = tvCoeffAmpPhSum(1).*ncsUncalibAmpPhTV(:,1) + ...
                tvCoeffAmpPhSum(2).*ncsUncalibAmpPhTV(:,2);
ncsTVAmp = tvCoeffAmp.*ncsUncalibAmpPhTV(:,1);
ncsTVPh = tvCoeffPh.*ncsUncalibAmpPhTV(:,2);

figure('Units', 'pixels', ...
    'Position', [100 100 900 500]);
ax2 = gca;
plot(tTV,hxTV,'color',[0 0.9 0],'LineWidth',2);
hold on
plot(tTV,ncsTVAmpPhSum,'color',[0.5 0.2 0.7],'LineStyle',':','LineWidth',2);
plot(tTV,ncsTVAmp,'color',[0.5,0.7,0.7],'LineWidth',2);
plot(tTV,ncsTVPh,'--','color',[0.9 0.5 0.2]);
hold off
grid on

xLabel = 'Time (sec)';
yLabel = 'Tidal Volume (mL)';
plotTitle = 'Estimated TV from Hexoskin and calibrated TV from NCS';
% plotLegend = {'Hx TV','Ncs TV: A*amp + B*ph','NCS TV: C*amp'};

plotLegend = {'Hx TV','Ncs TV: A*amp + B*ph','NCS TV: C*amp','NCS TV: D*ph'};
plotCute1(xLabel,yLabel,ax2,plotTitle,plotLegend,1);
axis(ax2,'tight')

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS breaths with BR

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);

nFig = 2;
ax3(1) = subplot(nFig,1,1);
yyaxis left
plot(ax3(1),tResp,hxRespTh./max(hxRespTh),':','color','k'); % ./max(hxRespTh)
hold on
plot(ax3(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]); % ./max(hxRespAbd)
plotCute1([],'a.u.',ax3(1),[],[],0);
ylim([0.995,1.001])
hold off

yyaxis right
plot(ax3(1),tBR,hxBR);
plotCute1('Time (sec)','BPM',ax3(1),...
    'Hexoskin Respiration (thorax & abdomen) and Breath Rate',{'Thorax (a.u.)',...
    'Abdomen (a.u.)','Breath Rate (BPM)'},1);

ax3(2) = subplot(nFig,1,2);
yyaxis left
plot(ax3(2),tResp,ncsResp(:,1)); 
plotCute1([],'a.u.',ax3(2),[],[],0);

yyaxis right
plot(ax3(2),tResp,ncsResp(:,2)); 
plotCute1('Time (sec)','a.u',ax3(2),...
    'NCS Respiration (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1);
 
linkaxes(ax3,'x')

%% ------------------------------------------------------------------------
% Calculating breath rate from NCS and comparing against Hexoskin. 

[ncsBR,tNcsBR] = ncsEstBR(ncsResp,inExAmp,inExPh,ncsRespSampRate,tBR);

figure('Units', 'pixels', 'Position', [100 100 900 500]);
ax4 = gca;
plot(tBR,hxBR,'LineWidth',2,'color',[0.5,0.2,0.7],'LineStyle',':');
hold on
plot(tNcsBR,ncsBR(:,1),'color',[0 0.9 0],'LineWidth',2);
plot(tNcsBR,ncsBR(:,2),'--','color',[0.9 0.5 0.2],'LineWidth',2);
hold off
grid on

xLabel = 'Time (s)';
yLabel = 'Breath Per Minute (BPM)';
plotTitle = 'Estimated BR from Hexoskin NCS';
plotLegend = {'Hx BR','NCS Amp BR','NCS Ph BR'};
plotCute1(xLabel,yLabel,ax4,plotTitle,plotLegend,1);
axis(ax4,'tight')

%% ------------------------------------------------------------------------
% Calculating fractional inspiratory time (Ti/Tt)
% Results using amplitude and phase waveform for NCS
[ncsAmpTiTt,ncsPhTiTt] = ncsEstTiTt(ncsResp,inExAmp,inExPh,ncsRespSampRate);
hxTiTt = hxEstTiTt(tHxInsp,tHxExp,hxInsp,hxExp);

figure('Units', 'pixels', 'Position', [100 100 900 500]);
ax5 = gca; 
plot(hxTiTt(:,1),hxTiTt(:,2),'LineWidth',2,'color',[0.5,0.2,0.7],'LineStyle',':');
hold on
plot(ncsAmpTiTt(:,1),ncsAmpTiTt(:,2),'color',[0 0.9 0],'LineWidth',2);
plot(ncsPhTiTt(:,1),ncsPhTiTt(:,2),'--','color',[0.9 0.5 0.2],'LineWidth',2);
hold off
grid on

xLabel = 'Time (s)';
yLabel = 'Ti/Tt';
plotTitle = 'Fractional Inspiratory Time (Ti/Tt) from Hexoskin and NCS';
plotLegend = {'Hx','NCS Amp','NCS Ph'};
plotCute1(xLabel,yLabel,ax5,plotTitle,plotLegend,1);
axis(ax5,'tight')

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS heartbeat with HR
% close all;

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);

nFig = 2;
ax6(1) = subplot(nFig,1,1);
yyaxis left
plot(ax6(1),tHeart,hxEcg,':','color','k','LineWidth',2);

yyaxis right
plot(ax6(1),tHR,hxHR);
plotCute1('Time (s)','mL',ax6(1),...
    'Hexoskin Heartbeat',{'Hx ECG (a.u.)',...
    'Heart Rate (BPM)'},1);

ax6(2) = subplot(nFig,1,2);
yyaxis left
plot(ax6(2),tHeart,ncsHeart(:,1)); 
plotCute1([],'a.u.',ax6(2),[],[],0);

yyaxis right
plot(ax6(2),tHeart,ncsHeart(:,2)); 
plotCute1('Time (s)','a.u.',ax6(2),...
    'NCS Heartbeat (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1);
 
linkaxes(ax6,'x')

%% ------------------------------------------------------------------------
% Peak detection in NCS Heartbeat

freqRangeHR =  [0.5,220/60];
[heartAmpMinMax, heartPhMinMax] = findInhaleExhale(ncsHeart,ncsHeartSampRate,freqRangeHR,tHeart);
ncsHR = ncsEstHR(ncsHeart,heartAmpMinMax,heartPhMinMax,ncsHeartSampRate,tHR);

figure('Units', 'pixels', 'Position', [100 100 900 500]);
ax7 = gca;
plot(tHR,hxHR,'LineWidth',2,'color',[0.5,0.2,0.7],'LineStyle',':');
hold on
plot(tHR,ncsHR(:,1),'color',[0 0.9 0],'LineWidth',2);
plot(tHR,ncsHR(:,2),'--','color',[0.9 0.5 0.2],'LineWidth',2);
hold off
xLabel = 'Time (s)';
yLabel = 'Heart Rate (BPM)';
plotTitle = 'Estimated HR from Hexoskin NCS';
plotLegend = {'Hx HR','NCS Amp HR','NCS Ph HR'};
plotCute1(xLabel,yLabel,ax7,plotTitle,plotLegend,1);
axis(ax7,'tight')
