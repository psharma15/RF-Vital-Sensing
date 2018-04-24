% This program compares respiratory vital sign and its characteristics as
% obtained from Near-Field Coherent Sensing, against Hexoskin smart garment
% data. 
% April 24, 2018
% Pragya Sharma, ps847@cornell.edu

%% ------------------------------------------------------------------------
% Provide input to function.
dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
    'Hexoskin\Data\3'];
hxFolder = '\user_13412b';
% hxDataNum: 1 = resp_abd, 2 = resp_thrx, 3 = ecg, 4 = tidal_vol_raw,
% 5 = tidal_vol_adj, 6 = min_vent_raw, 7 = min_vent_adj
hxDataNumRespAbd = 1; % Reading abdomen respiration
hxDataNumRespTh = 2; % Reading thorax respiration
hxDataNumTV = 4; % Tidal Volume (TV)
hxDataNumBR = 6; % Breath Rate (BR)
ncsDataNum = 11; % data at different time instants
% Manual time offset is by observation. 
% 20.5 for '2', -1 for '3a: 1-13', 6 for '3b: 14:20'
manualTimeOffset = -1; % sec: This is by observation 
ncsTstart = 0; % Time is relative to NCS in seconds, keep 140 for stable
dataDuration = 0; % Leave last 30 sec abd data: 10*60-32-ncsTstart

%% ------------------------------------------------------------------------
% Reading NCS-synchronized Hx waveforms
% Reading thoracic and abdomen respiratory waveform
[ncsSync,~,ncsHighSampRate,hxRespTh,tAbsHxResp,hxRespSampRate] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumRespTh,ncsDataNum,...
    manualTimeOffset,dataDuration,ncsTstart);
[~,~,~,hxRespAbd,~,~] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumRespAbd,ncsDataNum,...
    manualTimeOffset,dataDuration,ncsTstart);

% Synchronize NCS and Hx tidal volume (TV) estimation data. 
[~,~,~,hxTV,tAbsHxTV,hxSampRateTV] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumTV,ncsDataNum,...
    manualTimeOffset,dataDuration,ncsTstart);

% Synchronizing NCS and Hx breathing rate (BR) estimation data.
% Expecting ncsSync to remain the same
[~,~,~,hxBR,tAbsHxBR,hxSampRateBR] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNumBR,ncsDataNum,...
    manualTimeOffset,dataDuration,ncsTstart);

% Synchronizing timing of Hx data at different time stamps, taking start
% time of highest sampled respiration (abdomen or thorax) as the reference.
tRef = tAbsHxResp(1,:);
tResp = etime(tAbsHxResp,tRef); % Time starts from 0
tTV = etime(tAbsHxTV,tRef); % Time starts relative to tHxResp
tOffsetTV = tTV(1);
tBR = etime(tAbsHxBR,tRef); % Time starts relative to tHxResp
tOffsetBR = tBR(1);

%% ------------------------------------------------------------------------
% Process the NCS data to get respiration waveform - both amplitude 
% and phase, after filtering out heartbeat.
% *********************************************************************** %
% Specify correct amplitude and phase sign, see if you can implement algo
% to detect this as well. 
% *********************************************************************** %
fprintf('Change amplitude and phase sign if needed ...\n');
ncsFlipData = [1,1]; % -1 to flip
[ncsRespFiltered,~,~,~] = postProcess(0,1,1.4,ncsSync,ncsHighSampRate,ncsFlipData);

%% ------------------------------------------------------------------------
% Downsample NCS to hx respiration sample rate = 128 Hz
ncsSampRate = hxRespSampRate;
ncsResp = resample(ncsRespFiltered,ncsSampRate,ncsHighSampRate);

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS breaths with TV
close all;

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
% tHxResp = 0:1/hx
nFig = 2;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(ax1(1),tResp,hxRespTh./max(hxRespTh),':','color','k','LineWidth',2);
hold on
plot(ax1(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]);
plotCute1([],'a.u.',ax1(1),[],[],0);
ylim([0.9967,1.001])
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
[inExAmp, inExPh] = findInhaleExhale(ncsResp,ncsSampRate);

%% ------------------------------------------------------------------------
% *********************************************************************** %
% Calibrate and Estimate tidal volume with same sampling frequency as Hx.
% Remember to provide correct calibration time. 
% *********************************************************************** %
% calibTime = [tHxTV(1), tHxTV(end)];
calibTime = [40, 100]; % 10-50

[tvCoeffAmpPhSum,tvCoeffAmp,tvCoeffPh,ncsUncalibAmpPhTV] = ...
    ncsEstTV(ncsResp,inExAmp,inExPh,hxTV,ncsSampRate,hxSampRateTV,tOffsetTV,calibTime);
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
plot(ax3(1),tResp,hxRespTh./max(hxRespTh),':','color','k');
hold on
plot(ax3(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]);
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

[ncsBR,tNcsBR] = ncsEstBR(ncsResp,inExAmp,inExPh,ncsSampRate,tBR);

figure('Units', 'pixels', ...
    'Position', [100 100 900 500]);
ax4 = gca;
plot(tBR,hxBR,'LineWidth',2,'color',[0.5,0.2,0.7],'LineWIdth',2,'LineStyle',':');
hold on
plot(tNcsBR,ncsBR(:,1),'color',[0 0.9 0],'LineWidth',2);
plot(tNcsBR,ncsBR(:,2),'--','color',[0.9 0.5 0.2],'LineWidth',2);
hold off
grid on

xLabel = 'Time (sec)';
yLabel = 'Breath Per Minute (BPM)';
plotTitle = 'Estimated BR from Hexoskin NCS';
plotLegend = {'Hx BR','NCS Amp BR','NCS Ph BR'};
plotCute1(xLabel,yLabel,ax4,plotTitle,plotLegend,1);
axis(ax4,'tight')

%% ------------------------------------------------------------------------
% Calculating fractional inspiratory time (Ti/Tt)
% Resutls using amplitude and phase waveform for NCS
[ampTiTt,phTiTt] = ncsEstTiTt(ncsResp,inExAmp,inExPh,ncsSampRate);

%% ------------------------------------------------------------------------
% Results to report: 
% RMSE (Root mean squared Error):
