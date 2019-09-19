% Reading NCS data using our USRP 2 channel - one near heart, and one on
% the abdomen. Performing some basic analysis of that, including:
% Heartbeat processing, heart rate, HRV etc.
% Respiration waveform processing, breath rate, lung volume, TiTt etc.
% Possible motion information etc.
% Pragya Sharma, ps847@cornell.edu, 08 March 2019

%% ------------------------------------------------------------------------
% Provide input to function.
dataPath = 'E:\NCS\HumanStudyMarch2019\Data\Case1\';
ncsFile = '0414_164811Routine 1'; % data at different time instants
bioFile = 'Case1_14April_biopac';
fsNcsHigh = 50e3; % NCS sampling rate in Hz
fsBioHigh = 2e3;
tStartEndOff = [5,0]; % Start and end offset wrt NCS data in seconds
tManualOff = 0; % Manual offset between NCS and BIOPAC in seconds
ifCalib = 1; % This data is used to calibrate
ncsCalibFile = '0414_162409Calib1';

%% ------------------------------------------------------------------------
% Reading synchronized data
[ncsHS,bioHS,micHS] = ncsBioSync(dataPath,ncsFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tStartEndOff(1),tStartEndOff(2),...
                                 tManualOff);

tNcsHS = (0:(length(ncsHS)-1))/fsNcsHigh;
tBioHS = (0:(length(bioHS)-1))/fsBioHigh;

figure
nFig = 3;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(tNcsHS,-ncsHS(:,3));xlabel('Time (s)'); ylabel('NCS Amp');grid on;
yyaxis right
plot(tNcsHS,-ncsHS(:,4)); ylim([-200,200]); xlabel('Time (s)'); ylabel('NCS Ph');
ax1(2) = subplot(nFig,1,2);
plot(tBioHS,bioHS(:,2)); hold on;
plot(tBioHS,bioHS(:,3)); grid on;
xlabel('Time (s)'); ylabel('BIOPAC resp (mV)');
legend('Thorax','Abdomen'); hold off;
ax1(3) = subplot(nFig,1,3);
plot(tBioHS,bioHS(:,4));grid on;
xlabel('Time (s)'); ylabel('BIOPAC airflow (L/s)');
linkaxes(ax1,'x');

%% ------------------------------------------------------------------------
% if ifCalib == 1
    
%% ------------------------------------------------------------------------
% Process the NCS data to get respiration waveform - both amplitude 
% and phase, after filtering out heartbeat, and get Heartbeat data after
% filtering out respiration
% *********************************************************************** %

% Perform optional downsampling before processing
ncsSampRate = 500;
ncsData = resample(ncsDataHS,ncsSampRate,fsNcsHigh);

tNcsHS = 0:1/fsNcsHigh:((length(ncsDataHS)-1)/fsNcsHigh);
t = 0:1/ncsSampRate:((length(ncsData)-1)/ncsSampRate);
figure
plot(tNcsHS,ncsDataHS(:,3))
hold on
plot(t,ncsData(:,3))
legend('50k Samp/sec','500 Samp/sec')

% Optional data ordering, somehow data at end is first.
% ncsData = ncsData(end:-1:1,:);
tLim = [0.2,t(end) - 0.2];
ncsFlipData = [1,1]; % -1 to flip
fprintf('Change Ch 1 amplitude and phase sign by [%d,%d]...\n',ncsFlipData(1),ncsFlipData(2));
[ncsCh1resp,~,~,tTrunc] = postProcess(0,0.9,2.5,ncsData(:,1:2),ncsSampRate,ncsFlipData,tLim(1),tLim(2));
[ncsCh1heart,~,~,~] = postProcess(0.7,20,26,ncsData(:,1:2),ncsSampRate,ncsFlipData,tLim(1),tLim(2));
ncsFlipData = [1,1]; % -1 to flip
fprintf('Change Ch 2 amplitude and phase sign by [%d,%d]...\n',ncsFlipData(1),ncsFlipData(2));
[ncsCh2resp,~,~,~] = postProcess(0.0,2,5,ncsData(:,3:4),ncsSampRate,ncsFlipData,tLim(1),tLim(2));
[ncsCh2heart,~,~,~] = postProcess(0.5,15,20,ncsData(:,3:4),ncsSampRate,ncsFlipData,tLim(1),tLim(2));

%% ------------------------------------------------------------------------
% Plotting NCS Ch1 Ch2: heart and respiration

figure('Units', 'pixels', ...
    'Position', [100 100 530 360]);
nFig = [2, 2];
ax1(1) = subplot(nFig(1),nFig(2),1);
yyaxis left
plot(ax1(1),tTrunc,ncsCh1heart(:,1),'color',[0.6,0.2,0.6],'LineWidth',1); 
plotCute1('Time (s)','Ch1 Amp(a.u.)',ax1(1),[],[],0);
% yyaxis right
% plot(ax1(1),tTrunc,ncsCh1heart(:,2),'-','color',[0.3,0.75,0.93]); 
% plotCute1('Time (s)','Ch1 Ph(a.u.)',ax1(1),'Heart',{'Amp','Ph'},1);

ax1(2) = subplot(nFig(1),nFig(2),2);
yyaxis left
plot(ax1(2),tTrunc,ncsCh2resp(:,1)); 
plotCute1('Time (s)','Ch2 Amp(a.u.)',ax1(2),'Respiration',[],0);
% yyaxis right
% plot(ax1(2),tTrunc,ncsCh2resp(:,2)); 
% plotCute1('Time (s)','Ch2 Ph(degrees)',ax1(2),'Respiration',{'Amp','Ph'},1);

ax1(3) = subplot(nFig(1),nFig(2),3);
yyaxis left
plot(ax1(3),tTrunc,ncsCh2heart(:,1)); 
plotCute1('Time (s)','Ch2 Amp(a.u.)',ax1(3),'Heart',[],0);
% yyaxis right
% plot(ax1(3),tTrunc,ncsCh2heart(:,2)); 
% plotCute1('Time (s)','Ch2 Ph(degrees)',ax1(3),'Heart',{'Amp','Ph'},1);

ax1(4) = subplot(nFig(1),nFig(2),4);
yyaxis left
plot(ax1(4),tTrunc,ncsCh1resp(:,1)); 
plotCute1('Time (s)','Ch1 Amp(a.u.)',ax1(4),'Respiration',[],0);
% yyaxis right
% plot(ax1(4),tTrunc,ncsCh1resp(:,2)); 
% plotCute1('Time (s)','Ch1 Ph(degrees)',ax1(4),'Respiration',{'Amp','Ph'},1);

linkaxes(ax1,'x')

%% ------------------------------------------------------------------------
% Find 'inhalation end' and 'exhalation end' indices in amplitude and phase
% inExAmp: Ending of inhalation - peak
% inExPh: Ending of exhalation - minima
% Considering respiration is between 8-60 breaths per minute - normal for
% adult is 8-20 bpm.
% For the way findInhaleExhale function is written (to discard consecutive
% max and min), inhalation is supposed to be maxima (==1) and exhalation is
% supposed to be minima (==0).

freqRangeBR = [8, 80]./60; % In Hz
[inExAmpCh2, inExPhCh2] = findInhaleExhale(ncsCh2resp,ncsSampRate,freqRangeBR,tTrunc);

%% Detecting average peak-to-peak signal for a particular duration
inhaleAmpIdx = inExAmpCh2(inExAmpCh2(:,2)==0,1); 
exhaleAmpIdx = inExAmpCh2(inExAmpCh2(:,2)==1,1);
if length(inhaleAmpIdx) > length(exhaleAmpIdx)
    inhaleAmpIdx = inhaleAmpIdx(1:end-1);
elseif length(exhaleAmpIdx) > length(inhaleAmpIdx)
    exhaleAmpIdx = exhaleAmpIdx(1:end-1);
end
ncsCh2PkPk = ncsCh2resp(exhaleAmpIdx,1) - ncsCh2resp(inhaleAmpIdx,1);

meanPkPk = mean(ncsCh2PkPk(2:6));
fprintf('Mean Pk-Pk respiration during normal breathing is %f. \n',meanPkPk);


