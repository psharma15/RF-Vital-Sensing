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
dataPath = 'E:\NCS\Respiratory\NormalBreathing\Data\v2\Feb15_2019';
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

ncsFileName = '\freq1_8G_xiphoid_Neg5cmv2a 0215_1410'; % data at different time instants
ncsSampRate = 500; % NCS sampling rate in Hz
% Manual time offset is by observation. Use following settings:
% Yuna_Jul12: 18.4, Yuna_Jul18: 15.57, Yuna_Jul19: -.455
% Pragya_Jul16: -3.16
% Data4: -3.89
% Pragya_Jan31_2019: -3.812 for first two active ones (in order of time), 
% -46.682 for 1816
% Pragya_Feb15_2019: -1.982

manualTimeOffset = -1.982; % sec: This is by observation 
ncsTstart = 0; % Time is relative to NCS in seconds, keep 150 for stable
dataDuration = 0; % Leave last 30 sec abd data: 10*60-30-ncsTstart

% Plot props
legHorz = 'Horizontal'; 

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
% Specify correct amplitude and phase sign. This is IMPORTANT, mostly when
% finding fractional inspiratory time, where it is important to identify
% inspiration and expiration separately, and cannot be interchanged. The
% inspiration == minima points and expiration == maxima points.
% See if you can implement algo to detect this as well. 
% *********************************************************************** %
ncsFlipData = [1,1]; % -1 to flip

fprintf('Change amplitude and phase sign by %d and % d respectively...\n',ncsFlipData(1),ncsFlipData(2));

[ncsRespFiltered,~,~,~] = postProcess(0.0,5,10,ncsSync,ncsSampRate,ncsFlipData);
ncsRespDetrend = detrend(ncsRespFiltered); % Linear detrending

figure
yyaxis left
plot(ncsRespDetrend(:,1));
yyaxis right
plot(ncsRespDetrend(:,2));

ncsRespProcessed = ncsRespDetrend;

[ncsHeartFiltered,~,~,~] = postProcess(0.9,10,12,ncsSync,ncsSampRate,ncsFlipData);
ncsHeartProcessed = ncsHeartFiltered;

%% ------------------------------------------------------------------------
% Downsample NCS resp to hx respiration sample rate = 128 Hz
ncsRespSampRate = hxRespSampRate;
ncsResp = resample(ncsRespProcessed,ncsRespSampRate,ncsSampRate);

% Downsample NCS heart to hx ECG sample rate = 256 Hz
ncsHeartSampRate = hxEcgSampRate;
ncsHeart = resample(ncsHeartProcessed,ncsHeartSampRate,ncsSampRate);

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

% This is not right graph to plot in publication, because hxRespTh and 
% hxRespAbd are plotted on same y-axis, but are not comparable, as these
% are normalized to *different* values.
close all;

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
% tHxResp = 0:1/hx
nFig = 2;
ax0(1) = subplot(nFig,1,1);
yyaxis left
plot(ax0(1),tResp,hxRespTh./max(hxRespTh),':','color','k','LineWidth',2); %./max(hxRespTh)
hold on
plot(ax0(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]); %./max(hxRespAbd)
plotCute1('Time (s)','a.u.',ax0(1),[],[],0);
% ylim([0.9967,1.001])
hold off

yyaxis right
plot(ax0(1),tTV,hxTV);
plotCute1('Time (s)','mL',ax0(1),...
    'Hexoskin Respiration (thorax & abdomen) and Tidal Volume',{'Thorax (a.u.)',...
    'Abdomen (a.u.)','Tidal Volume (mL)'},1);

ax0(2) = subplot(nFig,1,2);
yyaxis left
plot(ax0(2),tResp,ncsResp(:,1)); 
plotCute1('Time (s)','a.u.',ax0(2),[],[],0);

yyaxis right
plot(ax0(2),tResp,ncsResp(:,2)); 
plotCute1('Time (s)','a.u.',ax0(2),...
    'NCS Respiration (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1);
 
linkaxes(ax0,'x')

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS breaths

figure('Units', 'pixels', ...
    'Position', [100 100 530 360]);
nFig = 2;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(ax1(1),tResp,hxRespTh,'color',[0.6,0.2,0.6],'LineWidth',1); 
plotCute1('Time (s)','Thorax (mL)',ax1(1),[],[],0);
yyaxis right
plot(ax1(1),tResp,hxRespAbd,'-','color',[0.3,0.75,0.93]); 
plotCute1('Time (s)','Abdomen (mL)',ax1(1),...
    'Hexoskin Respiration (thorax & abdomen)',{'Thorax (a.u.)',...
    'Abdomen (a.u.)'},1,legHorz);

ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(ax1(2),tResp,ncsResp(:,1)); 
plotCute1('Time (s)','a.u.',ax1(2),[],[],0);

yyaxis right
plot(ax1(2),tResp,ncsResp(:,2)); 
plotCute1('Time (s)','a.u.',ax1(2),...
    'NCS Respiration (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1,legHorz);
 
linkaxes(ax1,'x')

%% ------------------------------------------------------------------------
% Find 'inhalation end' and 'exhalation end' indices in amplitude and phase
% inExAmp: Ending of inhalation - peak
% inExPh: Ending of exhalation - minima
% Considering respiration is between 8-60 breaths per minute - normal for
% adult is 8-20 bpm.
freqRangeBR = [6, 60]./60; % In Hz
[inExAmp, inExPh] = findInhaleExhale(ncsResp,ncsRespSampRate,freqRangeBR,tResp);
[inExAmpHx, inExPhHx] = findInhaleExhale([hxRespTh,hxRespAbd],hxRespSampRate,freqRangeBR,tResp);

%% ------------------------------------------------------------------------
% *********************************************************************** %
% Calibrate and Estimate tidal volume with same sampling frequency as Hx.
% Remember to provide correct calibration time. 
% *********************************************************************** %
% calibTime = [tHxTV(1), tHxTV(end)];
calibTime = [200, 350]; 

fitFunc = 'quad';
[tvCoeffAmpPhSum,tvCoeffAmp,tvCoeffPh,ncsUncalibAmpPhTV] = ...
    ncsEstTV(ncsResp,inExAmp,inExPh,hxTV,ncsRespSampRate,hxSampRateTV,fitFunc,tOffsetTV,calibTime);

% [tvCoeffAmpPhSum,tvCoeffAmp,tvCoeffPh,ncsUncalibAmpPhTV] = ...
%     ncsEstTV2(ncsResp,inExAmp,inExPh,hxTV,ncsRespSampRate,hxSampRateTV,...
%     tOffsetTV,calibTime,[40,80],0.2);
% 
ncsTVAmpPhSum = tvCoeffAmpPhSum(1).*ncsUncalibAmpPhTV(:,1) + ...
                tvCoeffAmpPhSum(2).*ncsUncalibAmpPhTV(:,2);
            
switch(fitFunc)
    case 'linear'
        ncsTVAmp = tvCoeffAmp.*ncsUncalibAmpPhTV(:,1);
    case 'quad'
        ncsTVAmp = tvCoeffAmp(1).*(ncsUncalibAmpPhTV(:,1).^2)+...
                   tvCoeffAmp(2).*(ncsUncalibAmpPhTV(:,1))+tvCoeffAmp(3);
    otherwise
        fprintf('Enter correct fitFunc value.\n');
end
ncsTVPh = tvCoeffPh.*ncsUncalibAmpPhTV(:,2);

figure('Units', 'pixels', ...
    'Position', [100 100 670 240]);
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
plotCute1(xLabel,yLabel,ax2,plotTitle,plotLegend,1,legHorz);
axis(ax2,'tight')

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS breaths with BR
% This is not right graph to plot in publication, because hxRespTh and 
% hxRespAbd are plotted on same y-axis, but are not comparable, as these
% are normalized to *different* values.
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
% Take care of which algo to follow and also both hexoskin waveforms may
% give DIFFERENT resutls
[ncsBR,tNcsBR] = ncsEstBR2(ncsResp,inExAmp,inExPh,ncsRespSampRate,tBR,[40 50],0.01);
[hxBR2,tHxBR2] = ncsEstBR2([hxRespTh,hxRespAbd],inExAmpHx,inExPhHx,ncsRespSampRate,tBR,[40 50],0.01);
% [ncsBR,tNcsBR] = ncsEstBR(ncsResp,inExAmp,inExPh,ncsRespSampRate,tBR);
% [hxBR2,tHxBR2] = ncsEstBR([hxRespTh,hxRespAbd],inExAmpHx,inExPhHx,ncsRespSampRate,tBR);
% ****** Updating Hexosking BR calculation.

figure('Units', 'pixels', 'Position', [100 100 550 200]);
ax4 = gca;
plot(tHxBR2,hxBR2(:,2),'LineWidth',2,'color',[0.5,0.2,0.7],'LineStyle',':');
hold on
plot(tNcsBR,ncsBR(:,1),'color',[0 0.9 0],'LineWidth',2);
plot(tNcsBR,ncsBR(:,2),'--','color',[0.9 0.5 0.2],'LineWidth',2);
hold off
grid on

xLabel = 'Time (s)';
yLabel = 'Breath Per Minute (BPM)';
plotTitle = 'Estimated BR from Hexoskin NCS';
plotLegend = {'Hx BR','NCS Amp BR','NCS Ph BR'};
plotCute1(xLabel,yLabel,ax4,plotTitle,plotLegend,1,legHorz);
axis(ax4,'tight')

%% ------------------------------------------------------------------------
% Calculating fractional inspiratory time (Ti/Tt)
% Results using amplitude and phase waveform for NCS
[ncsAmpTiTt,ncsPhTiTt] = ncsEstTiTt(ncsResp,inExAmp,inExPh,ncsRespSampRate);
hxTiTt = hxEstTiTt(tHxInsp,tHxExp,hxInsp,hxExp);

figure('Units', 'pixels', 'Position', [100 100 550 200]);
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
plotCute1(xLabel,yLabel,ax5,plotTitle,plotLegend,1,legHorz);
axis(ax5,'tight')

%% ------------------------------------------------------------------------
% Plotting Hexoskin and NCS heartbeat with HR
% close all;

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);

nFig = 2;
ax6(1) = subplot(nFig,1,1);
yyaxis left
plot(ax6(1),tHeart,hxEcg,':','color',[0.3 1 0.5],'LineWidth',2);

yyaxis right
plot(ax6(1),tHR,hxHR);
plotCute1('Time (s)','BPM',ax6(1),...
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

ifWavelet = 1;
if ifWavelet == 1
    % Going for wavelet decomposition when heart waveform is poor - this would
    % lose HRV but can improve HR estimation
    wName = 'db10';
    recCoeff = 7; % Reconstructed coefficient
    [ncsHeartAmpWavelet,~,~] = wavedecrec1(ncsHeart(:,1),tHeart,wName,recCoeff,0);
    [ncsHeartPhWavelet,~,~] = wavedecrec1(ncsHeart(:,2),tHeart,wName,recCoeff,0);
    ncsHeartProcessed = [ncsHeartAmpWavelet,ncsHeartPhWavelet];
else
    ncsHeartProcessed = ncsHeart;
end

freqRangeHR =  [0.5,220/60];
[heartAmpMinMax, heartPhMinMax] = findInhaleExhale(ncsHeartProcessed,ncsHeartSampRate,freqRangeHR,tHeart);
ncsHR = ncsEstHR(ncsHeartProcessed,heartAmpMinMax,heartPhMinMax,ncsHeartSampRate,tHR);

figure('Units', 'pixels', 'Position', [100 100 530 180]);
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
plotCute1(xLabel,yLabel,ax7,plotTitle,plotLegend,1,legHorz);
axis(ax7,'tight')

% plotExport(['hr',ncsFileName(2:end)])