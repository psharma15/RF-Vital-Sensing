% This file reads the excel for stress test reaction time (RT) and
% synchronizes it with NCS file name. There is option to do volume
% calibration. Other processing includes respiration rate and heart rate
% estimation. Then HRV features are estimated.

% Pragya Sharma, ps847@cornell.edu
% Created on: 14 June 2019

% function [] = ncsReactionTimeProcess(dataPath,rtFile,ncsData,bioData,fs,...
%               tStartOff,tEndOff, tManualOff)
% -------------------------------------------------------------------------
% Inputs:
% dataPath: Path for NCS and BIOPAC data
% rtFile: Name of the reaction time stress/ attention test file
% ncsData: NCS data (pre-synced with Biopac)
% bioData: Biopac data (pre-synced with NCS)
% fs: Sampling frequency of NCS in Hz. Assuming same for NCS, Biopac.
% tStartOff: Start time offset selected during NCS/Biopac synchonization.
% tEndOff: End time offset selected during NCS/Biopac synchonization.
% tManualOFf: To manually offset the reaction time wrt NCS and Biopac.
% Since we are manually synchonizing with the Beep sound, there may be an
% offser.
% -------------------------------------------------------------------------
% Outputs:
% ncsTrial, bioTrial: NCS and biopac data synchonized to trial time.
% rtTrial: Reaction time data during trial arranged in column format as:
% [reaction_time_ms,  
% ncsTest, bioTest: NCS and biopac data synchronized to 
% -------------------------------------------------------------------------

% Provide input to function.
caseNum = 7;
dataPath = ['C:\Research\NCS\HumanStudyData\Case',num2str(caseNum),'\'];
ncsCalibFile = '0828_132142Calib4';
ncsRelaxFile = '0828_133802Routine2a'; % NCS file when listening to music
ncsAttnTestFile = '0828_134555Routine2b'; % NCS data during clock test
bioFile = 'bio_case8_2019-08-28T13_30_45';
bioCalibFile = 'bio_case8_2019-08-28T12_17_27';
rtFile = 'attn_case8_0828';
fsNcsHigh = 50e3; % NCS sampling rate in Hz
fsBioHigh = 2e3;

tOffNcsStEndRelax = [5,2]; % Start and end offset wrt NCS data in seconds. 
tOffNcsStEndAttn = [5,10]; 
tOffNcsStEndCalib = [3,3]; 

tOffNcsThAbd = [0,0,0]; %[Relax, Attn, Calib] Manual offset between NCS Th/Abd in seconds. +ve means Abd is shifted to right.
tOffNcsBio = [0,0,0]; % [Relax, Attn, Calib] Manual offset between NCS and BIOPAC in seconds. Biopac is offset wrt NCS.

signNcs = [-1 -1 -1 -1]; signNcsCalib = [-1 -1 -1 -1]; % sign[ampTh phTh ampAbd phAbd]

%% ------------------------------------------------------------------------
% Define conditions for optional processing here
ifVolCalib = 0; % If some associated calibration data is needed
%ifLoadCalibCoeff = 0; % If pre-saved calibration coefficient are present
ifDownSamp = [1 1]; % If [ncs biopac] data is to be downsampled
fsDS = [500, 500]; % [Ncs,Bio] downsampling frequencies
unwrapPh = [1 1];
unwrapPhCalib = unwrapPh;

%% ------------------------------------------------------------------------
% Reading synchronized data
[ncsRelax,bioRelax,~,~,tRelax,~,fig(2)] = ncsBioSync(dataPath,ncsRelaxFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tOffNcsStEndRelax,...
                                 tOffNcsBio(1),ifDownSamp,fsDS,signNcs,unwrapPh);

[ncsAttn,bioAttn,~,fs,tAttn,~,fig(1)] = ncsBioSync(dataPath,ncsAttnTestFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tOffNcsStEndAttn,...
                                 tOffNcsBio(2),ifDownSamp,fsDS,signNcs,unwrapPh);
                             
          
[ncsCalib,bioCalib,~,~,tCalib,~,fig(3)] = ...
    ncsBioSync(dataPath,ncsCalibFile,bioCalibFile,fsNcsHigh,...
               fsBioHigh,tOffNcsStEndCalib,...
               tOffNcsBio(3),ifDownSamp,fsDS,signNcsCalib,unwrapPhCalib);

%% Optional: Using cross correlation to estimate time shift, in case
% manually is difficult.
tCorr = [10, 50]; 
nStart = tCorr(1)*fs(1); % Assuming same fs for both ncs and biopac
nEnd = tCorr(2)*fs(1); 
[r, lags] = xcorr(bioAttn(nStart:nEnd,3),ncsAttn(nStart:nEnd,3));
figure; plot(lags,r)
[rMax,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDevCalib = lagsMax/fs(1);
fprintf('Suggested NCS calibration time offset is %f\n',tDevCalib);   
 
%%
opts1.filtType = 'Hp';
opts1.f3db = 0.1; 
bioCalib(:,2:3) = filterLpHp(bioCalib(:,2:3),fs(2),opts1);
bioRelax(:,2:3) = filterLpHp(bioRelax(:,2:3),fs(2),opts1);
bioAttn(:,2:3) = filterLpHp(bioAttn(:,2:3),fs(2),opts1);

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.1; opts2.fpLP = 0.9; opts2.fstLP = 1.5;

ncsCalib(:,3) = filterLpHp(ncsCalib(:,3),fs(1),opts2);
ncsCalib(:,1) = filterLpHp(ncsCalib(:,1),fs(1),opts2);
ncsCalib(:,2) = filterLpHp(ncsCalib(:,2),fs(1),opts2);
ncsCalib(:,4) = filterLpHp(ncsCalib(:,4),fs(1),opts2);
ncsRelaxRespAbd = filterLpHp(ncsRelax(:,3),fs(1),opts2);
ncsRelaxRespTh = filterLpHp(ncsRelax(:,1),fs(1),opts2);
ncsAttnRespAbd = filterLpHp(ncsAttn(:,3),fs(1),opts2);
ncsAttnRespTh = filterLpHp(ncsAttn(:,1),fs(1),opts2);

%% ------------------------------------------------------------------------
% Calibration and volume estimation phase
if ifVolCalib == 1
    % If calibration is selected, we calculate lung volume or tidal volume
    % from BIOPAC airflow data.

    % Baseline estimation: Mean deviation from zero when no airflow.
    % Using the routine's data to estimate that.
    tMeanDev = [20, 80]; % Time window in sec to estimate baseline
    figure
    plot(tAttn(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2)),bioAttn(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2),4))
    meanDev = mean(bioAttn(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2),4));    
    
    % TV estimated at the expiration end for each breath cycle using airflow. 
    % Cycle starts from an inspiration and ends with expiration.
    [tvAirflow,volAirflow,fig(4)] = bioAirflowTV(bioCalib(:,4),fs(2),meanDev);  

    % TV calibration of Biopac chest belts
    opts1.calibType = 'vol'; opts1.fitEqn = 'BiasedLinear'; 
    opts1.tWin = 2; opts1.minInterceptDist = 0.1;

    % For volume calibration, signals need to be high pass filtered to
    % remove baseline: both biopac belts, and volume airflow
    volAirflow = filterLpHp(volAirflow,fs(2),opts1);
    
    % Start and stop calibration times: Performing calibration during
    % normal breathing only: ONLY FOR calibType='vol'
    opts1.tCalib = [2,12]; 
    [beltCalibCoeff,vBeltCalib] = bioBeltVolCalib(bioCalib(:,2:3),fs(2),tvAirflow,volAirflow,opts1); 
    
    
    opts2.tCalib = opts1.tCalib; 
    opts2.calibType = 'vol';     opts2.fitEqn = 'AbdLinear'; 
    opts2.tWin = 4;              opts2.minInterceptDist = 0.2;
    [ncsCalibCoeff,vNcsCalib] = ncsVolCalib(ncsCalib(:,3),fs(1),tvAirflow,volAirflow,opts2);
    
    % Optional: Using cross correlation to estimate time shift, in case
    % manually is difficult.
    tCorr = [2, 12]; 
    nStart = tCorr(1)*fs(1); % Assuming same fs for both ncs and biopac
    nEnd = tCorr(2)*fs(1); 
    [r, lags] = xcorr(bioCalib(nStart:nEnd,3),ncsCalib(nStart:nEnd,3));
    figure; plot(lags,r)
    [rMax,rMaxIdx] = max(abs(r));
    lagsMax = lags(rMaxIdx);
    tDevCalib = lagsMax/fs(1);
    fprintf('Suggested NCS calibration time offset is %f\n',tDevCalib);
    
    
    % Calibration: Plot Biopac and NCS volumes compared to airflow volume
    fig(5) = figure('Position',[400 200 800 600]);
    nFig = 3;
    ax1(1) = subplot(nFig,1,1);
    plot(tCalib,bioCalib(:,2))
    hold on; plot(tCalib,bioCalib(:,3));
    plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
    ax1(2) = subplot(nFig,1,2);
    plot(tCalib,ncsCalib(:,1));
    hold on; plot(tCalib,ncsCalib(:,3));
    plotCute1('Time (s)','NCS (V)',ax1(2),[],{'NCS Th','NCS Abd'},1);
    ax1(3) = subplot(nFig,1,3);
    plot(tCalib,volAirflow,'color',[79, 209, 98]/256,'LineWidth',2); hold on;
    plot(tCalib,vBeltCalib,'color',[15, 36, 193]/256); 
    plot(tCalib,vNcsCalib,'color',[229, 42, 25]/256)
    leg = {'Airflow','Belt','NCS'};
    plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1,'Horizontal');
    linkaxes(ax1,'x')
    
    % ------------------------------------------------------------------------
    % Use generated calibration coefficients to estimate Biopac and NCS Vol/TV
    opts1.tWin = 4; opts1.minInterceptDist = 0.2;
    vBeltAttn = bioBeltVol(bioAttn(:,2:3),beltCalibCoeff,fs(2),opts1);

    opts2.tWin = 4;
    vNcsAttn = ncsVol(ncsAttnRespAbd,ncsCalibCoeff,fs(1),opts2);

    vBeltRelax = bioBeltVol(bioRelax(:,2:3),beltCalibCoeff,fs(2),opts1);

    opts2.tWin = 4;
    vNcsRelax = ncsVol(ncsRelaxRespAbd,ncsCalibCoeff,fs(1),opts2);

    % Plot calibrated Biopac and NCS volumes (both using calib coefficients)
    fig(6) = figure('Position',[400 200 600 600]);
    nFig = 3;
    ax1(1) = subplot(nFig,1,1);
    plot(tAttn,bioAttn(:,2))
    hold on; plot(tAttn,bioAttn(:,3));
    plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
    ax1(2) = subplot(nFig,1,2);
    plot(tAttn,ncsAttnRespTh);
    hold on; plot(tAttn,ncsAttnRespAbd);
    plotCute1('Time (s)','NCS (V)',ax1(2),[],{'NCS Th','NCS Abd'},1);
    ax1(3) = subplot(nFig,1,3);
    plot(tAttn,vBeltAttn); hold on;
    plot(tAttn,vNcsAttn)
    leg = {'V Belt','V NCS'};
    plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1);
    linkaxes(ax1,'x')

    % Plot calibrated Biopac and NCS volumes (both using calib coefficients)
    fig(7) = figure('Position',[400 200 600 600]);
    nFig = 3;
    ax1(1) = subplot(nFig,1,1);
    plot(tRelax,bioRelax(:,2))
    hold on; plot(tRelax,bioRelax(:,3));
    plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
    ax1(2) = subplot(nFig,1,2);
    plot(tRelax,ncsRelaxRespTh);
    hold on; plot(tRelax,ncsRelaxRespAbd);
    plotCute1('Time (s)','NCS (V)',ax1(2),[],{'NCS Th','NCS Abd'},1);
    ax1(3) = subplot(nFig,1,3);
    plot(tRelax,vBeltRelax); hold on;
    plot(tRelax,vNcsRelax)
    leg = {'V Belt','V NCS'};
    plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1);
    linkaxes(ax1,'x')
end

%% ------------------------------------------------------------------------
% Now looking at the reaction time data
optsRT = detectImportOptions([dataPath,rtFile,'.xlsx']);
rtData = readtable([dataPath,rtFile,'.xlsx'],optsRT);

tManualRToff = 0; % 0 sec offset
rtTime = rtData{:,1}/1000; % Converting ms to s
% Correcting time offset introduced above, with 15s offset due to 'click to
% start' waiting time. As that is manual, additional manual offset is done.
rtTime = rtTime - tOffNcsStEndAttn(1) + 15 + tManualRToff;
rtDataNew = [rtData,array2table(rtTime(:),'VariableNames',{'Time'})];

% Cases when it jumped and was correctly detected.
idxJumpDet = rtData{:,6} == 1;
% Cases when jump was missed.
idxJumpMiss = rtData{:,8} == 1;
% Cases when no jump, wrongly pressed space bar.
idxNoJumpWrong = rtData{:,7} == 1;

% NCS filtered to show heartbeat
opts2.filtType = 'LpHp';
opts2.f3db = 0.7; opts2.fpLP = 5; opts2.fstLP = 7;
ncsHeartTh = filterLpHp(ncsAttn(:,1),fs(1),opts2);
ncsHeartAbd = filterLpHp(ncsAttn(:,3),fs(2),opts2);

% NCS filtered to show respiration
ncsAttnRespTh = filterLpHp(ncsAttn(:,1),fs(1),opts2);

%% Respiration and BR estimation: Attention test
opts2.filtType = 'Hp'; 
opts2.f3db = 0.1; opts2.fpLP = 0.9; opts2.fstLP = 1.5;
ncsAttnRespAbd = filterLpHp(ncsAttn(:,3),fs(1),opts2);
ncsAttnRespAbdPh = filterLpHp(ncsAttn(:,4),fs(1),opts2);

opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
[brNcsAttn,ncsAttnRespPk] = brEst(ncsAttnRespAbd,fs(1),opts3);
pkMaxNcsAttn = ncsAttnRespPk(1).idx(ncsAttnRespPk(1).ind == 1);

[brBioAttnTh,bioAttnRespPkTh] = brEst(bioAttn(:,2),fs(1),opts3);
[brBioAttnAbd,bioAttnRespPkAbd] = brEst(bioAttn(:,3),fs(1),opts3);
pkMaxBioAttnTh = bioAttnRespPkTh(1).idx(bioAttnRespPkTh(1).ind == 1);
pkMaxBioAttnAbd = bioAttnRespPkAbd(1).idx(bioAttnRespPkAbd(1).ind == 1);

%% Plot original figure: The resolution of ncs, @500 Sps is 2ms, while
% reaction time changes on the order of 1ms. 
fig(8) = figure('Position',[100,100,600,600]);
nFig = 4;
ax1(1) = subplot(nFig,1,1);
yyaxis left
%plot(t,ncs(:,1)); hold on;
plot(tAttn,ncsAttnRespAbd);hold on
plot(tAttn(pkMaxNcsAttn),ncsAttnRespAbd(pkMaxNcsAttn),'^');
plotCute1('Time (s)','NCS (V)',ax1(1),[],[],1);

yyaxis right
plot(tAttn,ncsAttnRespAbdPh);
leg = {'Resp Amp','Peak','Resp Ph'};
plotCute1('Time (s)','NCS (deg)',ax1(1),[],leg,1);

ax1(2) = subplot(nFig,1,2);
plot(tAttn,bioAttn(:,2)); hold on;
plot(tAttn(pkMaxBioAttnTh),bioAttn(pkMaxBioAttnTh,2),'^');
plot(tAttn,bioAttn(:,3));
plot(tAttn(pkMaxBioAttnAbd),bioAttn(pkMaxBioAttnAbd,3),'^');

leg = {'Resp Th','Th Peak','Resp Abd','Abd Peak'};
plotCute1('Time (s)','Biopac (V)',ax1(2),[],leg,1);

ax1(3) = subplot(nFig,1,3);
plot(tAttn,brNcsAttn); hold on
plot(tAttn,brBioAttnTh); plot(tAttn,brBioAttnAbd);

leg = {'BR NCS','BR Bio Th','BR Abd Th'};
plotCute1('Time (s)','Breath Rate (BPM)',ax1(3),[],leg,1);

ax1(4) = subplot(nFig,1,4);
stem(rtTime(idxJumpDet),rtData{idxJumpDet,2},'go'); hold on
stem(rtTime(idxJumpMiss),rtData{idxJumpMiss,2},'r+');
stem(rtTime(idxNoJumpWrong),rtData{idxNoJumpWrong,2},'rs');
ylim([0,1500])
leg = {'True Jump','Missed Jump','Wrong'};
plotCute1('Time (s)','Reaction Time (ms)',ax1(4),[],leg,1);


linkaxes(ax1,'x');

%% Heartbeat and HR estimation
opts2.filtType = 'LpHp'; 
opts2.orderHP = 5; opts2.Ast = 20;
opts2.f3db = 0.7; opts2.fpLP = 1.5; opts2.fstLP = 1.8;
ncsAttnHrTh = filterLpHp(ncsAttn(:,1),fs(1),opts2);
ncsAttnHrThPh = filterLpHp(ncsAttn(:,2),fs(1),opts2);

opts3.tWinHR = 4;
opts3.tWin = 0.5;
opts3.minInterceptDist = 0.05;
[hrNcsAttn, ncsAttnHrPeak] = hrEst(ncsAttnHrThPh,fs(1),opts3);
pkMaxNcsAttn = ncsAttnHrPeak(1).idx(ncsAttnHrPeak(1).ind == 1);

[hrBioAttn, bioAttnHrPeak] = ecgHR(bioAttn(:,1),fs(2),opts3);

%% Plot original figure: The resolution of ncs, @500 Sps is 2ms, while
% reaction time changes on the order of 1ms. 
fig(9) = figure('Position',[200,200,600,600]);
nFig = 4;
ax1(1) = subplot(nFig,1,1);
%yyaxis left
%plot(t,ncs(:,1)); hold on;
plot(tAttn,ncsAttnHrThPh); hold on;
plot(tAttn(pkMaxNcsAttn),ncsAttnHrThPh(pkMaxNcsAttn),'k^');
leg = {'Heart','Peak'};
plotCute1('Time (s)','NCS (V)',ax1(1),[],leg,1);

ax1(2) = subplot(nFig,1,2);
plot(tAttn,bioAttn(:,1)); hold on;
plot(tAttn(bioAttnHrPeak),bioAttn(bioAttnHrPeak,1),'^');
leg = {'ECG','ECG Peak'};
plotCute1('Time (s)','Biopac ECG(V)',ax1(2),[],leg,1);

ax1(3) = subplot(nFig,1,3);
plot(tAttn,hrNcsAttn);
hold on
plot(tAttn,hrBioAttn);
leg = {'NCS','ECG'};
plotCute1('Time (s)','HR (BPM)',ax1(3),[],leg,1);

ax1(4) = subplot(nFig,1,4);
stem(rtTime(idxJumpDet),rtData{idxJumpDet,2},'go'); hold on
stem(rtTime(idxJumpMiss),rtData{idxJumpMiss,2},'r+');
stem(rtTime(idxNoJumpWrong),rtData{idxNoJumpWrong,2},'rs');
ylim([0,1500])
leg = {'True Jump','Missed Jump','Wrong'};
plotCute1('Time (s)','Reaction Time (ms)',ax1(4),[],leg,1);

linkaxes(ax1,'x');

%% Respiration and BR estimation: Relax test
opts2.filtType = 'Hp'; 
opts2.f3db = 0.1; opts2.fpLP = 0.9; opts2.fstLP = 1.5;
ncsRelaxRespAbd = filterLpHp(ncsRelax(:,3),fs(1),opts2);
ncsRelaxRespAbdPh = filterLpHp(ncsRelax(:,4),fs(1),opts2);

opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
[brNcsRelax,ncsRelaxRespPk] = brEst(ncsRelaxRespAbd,fs(1),opts3);
pkMaxNcsRelax = ncsRelaxRespPk(1).idx(ncsRelaxRespPk(1).ind == 1);

[brBioRelaxTh,bioRelaxRespPkTh] = brEst(bioRelax(:,2),fs(1),opts3);
[brBioRelaxAbd,bioRelaxRespPkAbd] = brEst(bioRelax(:,3),fs(1),opts3);
pkMaxBioRelaxTh = bioRelaxRespPkTh(1).idx(bioRelaxRespPkTh(1).ind == 1);
pkMaxBioRelaxAbd = bioRelaxRespPkAbd(1).idx(bioRelaxRespPkAbd(1).ind == 1);

%% Heartbeat and HR estimation

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.7; opts2.fpLP = 1.8; opts2.fstLP = 2.0;
ncsRelaxHrTh = filterLpHp(ncsRelax(:,1),fs(1),opts2);
ncsRelaxHrThPh = filterLpHp(ncsRelax(:,2),fs(1),opts2);

opts3.tWinHR = 4;
opts3.tWin = 0.8;
opts3.minInterceptDist = 0.2;
[hrNcsRelax, ncsRelaxHrPeak] = hrEst(ncsRelaxHrTh,fs(1),opts3);
pkMaxNcsRelax = ncsRelaxHrPeak(1).idx(ncsRelaxHrPeak(1).ind == 1);

[hrBioRelax, bioRelaxHrPeak] = ecgHR(bioRelax(:,1),fs(2),opts3);
 
%%
close all

%% 
% -------------------------------------------------------------------------
% Estimating HRV features for comparing relaxation and stress test
% -------------------------------------------------------------------------
opts.tWinHR = 3;
opts.tWin = 0.5;
opts.minInterceptDist = 0.1;

opts.minEcgPkHt = 0.1;

n = 2^nextpow2(length(ncsRelax(:,2)));
yRR = fft(ncsRelax(:,2),n);
f = fs(1)*(0:(n/2))/n;
P = abs(yRR)/n;

figure
plot(f,P(1:n/2+1)); xlim([0.5, 3]) 
title('Frequency spectrum NCS relax')

n = 2^nextpow2(length(ncsAttn(:,1)));
yRR = fft(ncsAttn(:,1),n);
f = fs(1)*(0:(n/2))/n;
P = abs(yRR)/n;

figure
plot(f,P(1:n/2+1));  xlim([0.5, 3])  
title('Frequency spectrum NCS attn')

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.6; opts2.fpLP = 1.5; opts2.fstLP = 1.9;
% Heartbeat signal
ncsRelaxHrTh = filterLpHp(ncsRelax(:,1),fs(1),opts2);
ncsRelaxHrThPh = filterLpHp(ncsRelax(:,2),fs(1),opts2);
ncsRelaxHrAbd = filterLpHp(ncsRelax(:,3),fs(1),opts2);
ncsRelaxHrAbdPh = filterLpHp(ncsRelax(:,4),fs(1),opts2);

ncsAttnHrTh = filterLpHp(ncsAttn(:,1),fs(1),opts2);
ncsAttnHrThPh = filterLpHp(ncsAttn(:,2),fs(1),opts2);
ncsAttnHrAbd = filterLpHp(ncsAttn(:,3),fs(1),opts2);
ncsAttnHrAbdPh = filterLpHp(ncsAttn(:,4),fs(1),opts2);


% Filtering if T wave is very high and noisy, thus wrongly detected as R
% peak.
opts2.filtType = 'LpHp';
opts2.f3db = 4;
opts2.fpLP = 20; opts2.fstLP = 25;
ecgFiltRelax = filterLpHp(bioRelax(:,1),fs(2),opts2);
ecgFiltAttn = filterLpHp(bioAttn(:,1),fs(2),opts2);

opts.pkAmpRelRejThresh = 0.4; 
opts.tRRthresh = [0.4,1.2]; 

[hrvNcsRelax,hrvBioRelax,~,~] = hrvFeatureEst(ncsRelaxHrTh,ecgFiltRelax,fs,opts);
[hrvNcsAttn,hrvBioAttn,~,~] = hrvFeatureEst(ncsAttnHrTh,ecgFiltAttn,fs,opts);

% close all

%%
% During relaxed stage
figure('Position',[300,100,600,550])
nFig = 3;
ax1(1) = subplot(nFig,1,1);
plot(hrvNcsRelax.tRR,hrvNcsRelax.rrInterval);hold on
plot(hrvBioRelax.tRR,hrvBioRelax.rrInterval);
titl = ['Relaxation: Normal RR interval time, meanNCS: ',num2str(mean(hrvNcsRelax.rrInterval),3),...
        ' meanECG: ',num2str(mean(hrvBioRelax.rrInterval),3)];
leg = {'NCS','ECG'};
plotCute1('Time (s)','RR time (ms)',ax1(1),titl,leg,1);

ax1(2) = subplot(nFig,1,2);
plot(hrvNcsRelax.tRR,hrvNcsRelax.hr);
hold on;
plot(hrvBioRelax.tRR,hrvBioRelax.hr);
leg = {'NCS','ECG'};
titl = 'Heart Rate';
plotCute1('Time (s)','HR (BPM)',ax1(2),titl,leg,1,'Horizontal');

ax1(3) = subplot(nFig,1,3);
plot(hrvNcsRelax.fftFreqXaxis,hrvNcsRelax.fftPowYaxis); hold on
plot(hrvBioRelax.fftFreqXaxis,hrvBioRelax.fftPowYaxis);
leg = {'NCS','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax1(3),titl,leg,1,'Horizontal');

linkaxes(ax1(1:2),'x');

% During attention test

figure('Position',[100,100,600,650])
nFig = 3;
ax2(1) = subplot(nFig,1,1);
plot(hrvNcsAttn.tRR,hrvNcsAttn.rrInterval);hold on
plot(hrvBioAttn.tRR,hrvBioAttn.rrInterval);
titl = ['Attention: Normal RR interval time, meanNCS: ',num2str(mean(hrvNcsAttn.rrInterval),3),...
        ' meanECG: ',num2str(mean(hrvBioAttn.rrInterval),3)];
leg = {'NCS','ECG'};
plotCute1('Time (s)','RR time (ms)',ax2(1),titl,leg,1);

ax2(2) = subplot(nFig,1,2);
plot(hrvNcsAttn.tRR,hrvNcsAttn.hr);
hold on;
plot(hrvBioAttn.tRR,hrvBioAttn.hr);
leg = {'NCS','ECG'};
titl = 'Heart Rate';
plotCute1('Time (s)','HR (BPM)',ax2(2),titl,leg,1,'Horizontal');


ax2(3) = subplot(nFig,1,3);
plot(hrvNcsRelax.fftFreqXaxis,hrvNcsRelax.fftPowYaxis); hold on
plot(hrvBioRelax.fftFreqXaxis,hrvBioRelax.fftPowYaxis);
leg = {'NCS','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax2(3),titl,leg,1,'Horizontal');

% ax2(4) = subplot(nFig,1,4);
% stem(rtTime(idxJumpDet),rtData{idxJumpDet,2},'go'); hold on
% stem(rtTime(idxJumpMiss),rtData{idxJumpMiss,2},'r+');
% stem(rtTime(idxNoJumpWrong),rtData{idxNoJumpWrong,2},'rs');
% ylim([0,1500])
% leg = {'True Jump','Missed Jump','Wrong'};
% plotCute1('Time (s)','Reaction Time (ms)',ax2(4),[],leg,1);
% 
linkaxes(ax2([1:2]),'x');
%%
close all
%% Heartbeat second harmonic
opts2.f3db = 2.1; opts2.fpLP = 3.1; opts2.fstLP = 3.3;
ncsRelaxHrTh1Harm = filterLpHp(ncsRelax(:,1),fs(1),opts2);
ncsRelaxHrThPh1Harm = filterLpHp(ncsRelax(:,2),fs(1),opts2);
ncsAttnHrTh1Harm = filterLpHp(ncsAttn(:,1),fs(1),opts2);
ncsAttnHrThPh1Harm = filterLpHp(ncsAttn(:,2),fs(1),opts2);

opts.harmNum = 2;
% opts.tWin = 2;
% opts.minInterceptDist = 0.05; opts.pkAmpRelRejThresh = 0.6; 
[hrvNcsRelax2Harm,hrvBioRelax,~,~] = hrvFeatureEst(ncsRelaxHrTh1Harm,ecgFiltRelax,fs,opts);
[hrvNcsAttn2Harm,hrvBioAttn,~,~] = hrvFeatureEst(ncsAttnHrTh1Harm,ecgFiltAttn,fs,opts);

% hrvNcsRelax2HarmMat = [hrvNcsRelax2Harm.meanRR; hrvNcsRelax2Harm.sdnn; ...
%     hrvNcsRelax2Harm.meanHR; hrvNcsRelax2Harm.rmssd; ...
%     hrvNcsRelax2Harm.sdsdRR; hrvNcsRelax2Harm.pNN50; hrvNcsRelax2Harm.LFpow;...
%     hrvNcsRelax2Harm.HFpow; hrvNcsRelax2Harm.LFHFratio];
% 

%% Plot second harmonic HRV features
% During relaxed stage
figure('Position',[300,100,600,550])
nFig = 3;
ax1(1) = subplot(nFig,1,1);
hold on
plot(hrvNcsRelax2Harm.tRR,hrvNcsRelax2Harm.rrInterval);
plot(hrvBioRelax.tRR,hrvBioRelax.rrInterval);
titl = ['Normal RR interval time'];
leg = {'NCS','ECG'};
plotCute1('Time (s)','RR time (ms)',ax1(1),titl,leg,1);

ax1(2) = subplot(nFig,1,2);
hold on;
plot(hrvNcsRelax2Harm.tRR,hrvNcsRelax2Harm.hr);
plot(hrvBioRelax.tRR,hrvBioRelax.hr);
leg = {'NCS','ECG'};
titl = 'Heart Rate';
plotCute1('Time (s)','HR (BPM)',ax1(2),titl,leg,1,'Horizontal');

ax1(3) = subplot(nFig,1,3);
hold on
plot(hrvNcsRelax2Harm.fftFreqXaxis,hrvNcsRelax2Harm.fftPowYaxis);
plot(hrvBioRelax.fftFreqXaxis,hrvBioRelax.fftPowYaxis);
leg = {'NCS','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax1(3),titl,leg,1,'Horizontal');

linkaxes(ax1(1:2),'x');

% During attention test

figure('Position',[100,100,600,650])
nFig = 3;
ax2(1) = subplot(nFig,1,1);
hold on
plot(hrvNcsAttn2Harm.tRR,hrvNcsAttn2Harm.rrInterval);
plot(hrvBioAttn.tRR,hrvBioAttn.rrInterval);
titl = ['Normal RR interval time'];
leg = {'NCS','ECG'};
plotCute1('Time (s)','RR time (ms)',ax2(1),titl,leg,1);

ax2(2) = subplot(nFig,1,2);
hold on;
plot(hrvNcsAttn2Harm.tRR,hrvNcsAttn2Harm.hr);
plot(hrvBioAttn.tRR,hrvBioAttn.hr);
leg = {'NCS','ECG'};
titl = 'Heart Rate';
plotCute1('Time (s)','HR (BPM)',ax2(2),titl,leg,1,'Horizontal');


ax2(3) = subplot(nFig,1,3);
plot(hrvNcsAttn2Harm.fftFreqXaxis,hrvNcsAttn2Harm.fftPowYaxis); hold on
plot(hrvBioRelax.fftFreqXaxis,hrvBioRelax.fftPowYaxis);
leg = {'NCS','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax2(3),titl,leg,1,'Horizontal');

% ax2(4) = subplot(nFig,1,4);
% stem(rtTime(idxJumpDet),rtData{idxJumpDet,2},'go'); hold on
% stem(rtTime(idxJumpMiss),rtData{idxJumpMiss,2},'r+');
% stem(rtTime(idxNoJumpWrong),rtData{idxNoJumpWrong,2},'rs');
% ylim([0,1500])
% leg = {'True Jump','Missed Jump','Wrong'};
% plotCute1('Time (s)','Reaction Time (ms)',ax2(4),[],leg,1);
% 
linkaxes(ax2([1:2]),'x');
%%
close all;

%% Plot both fundamental and harmonic overlayed
% During relaxed stage
figure('Position',[300,100,600,700])
nFig = 3;
ax1(1) = subplot(nFig,1,1);
plot(hrvNcsRelax.tRR,hrvNcsRelax.rrInterval,'-','color',[1,0.7,0.4]);hold on
plot(hrvNcsRelax2Harm.tRR,hrvNcsRelax2Harm.rrInterval,'-','color',[0.2,0.6,1]);
plot(hrvBioRelax.tRR,hrvBioRelax.rrInterval,'-','color',[0.5,0.5,0.5]);
titl = ['Relaxation Exercise'];
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
plotCute1([],'RR time (ms)',ax1(1),titl,leg,1,'Horizontal');

ax1(2) = subplot(nFig,1,2);
plot(hrvNcsRelax.tRR,hrvNcsRelax.hr,'color',[1,0.7,0.4]); hold on;
plot(hrvNcsRelax2Harm.tRR,hrvNcsRelax2Harm.hr,'color',[0.2,0.6,1]);
plot(hrvBioRelax.tRR,hrvBioRelax.hr,'color',[0.5,0.5,0.5]);
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
titl = 'Heart Rate';
plotCute1('Time (s)','HR (BPM)',ax1(2),[],[],0);

ax1(3) = subplot(nFig,1,3);
plot(hrvNcsRelax.fftFreqXaxis,hrvNcsRelax.fftPowYaxis,'color',[1,0.7,0.4]); hold on
plot(hrvNcsRelax2Harm.fftFreqXaxis,hrvNcsRelax2Harm.fftPowYaxis,'color',[0.2,0.6,1]);
plot(hrvBioRelax.fftFreqXaxis,hrvBioRelax.fftPowYaxis,'color',[0.5,0.5,0.5]);
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax1(3),titl,[],0);

linkaxes(ax1(1:2),'x');

% During attention test

figure('Position',[100,100,600,700])
nFig = 4;
ax2(1) = subplot(nFig,1,1);
plot(hrvNcsAttn.tRR,hrvNcsAttn.rrInterval,'-','color',[1,0.7,0.4]);hold on
plot(hrvNcsAttn2Harm.tRR,hrvNcsAttn2Harm.rrInterval,'-','color',[0.2,0.6,1]);
plot(hrvBioAttn.tRR,hrvBioAttn.rrInterval,'-','color',[0.5,0.5,0.5]);
titl = ['Attention Test'];
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
plotCute1([],'RR time (ms)',ax2(1),titl,leg,1,'Horizontal');

ax2(2) = subplot(nFig,1,2);
plot(hrvNcsAttn.tRR,hrvNcsAttn.hr,'color',[1,0.7,0.4]);hold on;
plot(hrvNcsAttn2Harm.tRR,hrvNcsAttn2Harm.hr,'color',[0.2,0.6,1]);
plot(hrvBioAttn.tRR,hrvBioAttn.hr,'color',[0.5,0.5,0.5]);
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
titl = 'Heart Rate';
plotCute1([],'HR (BPM)',ax2(2),[],[],0);


ax2(4) = subplot(nFig,1,4);
plot(hrvNcsAttn.fftFreqXaxis,hrvNcsAttn.fftPowYaxis,'color',[1,0.7,0.4]);hold on
plot(hrvNcsAttn2Harm.fftFreqXaxis,hrvNcsAttn2Harm.fftPowYaxis,'color',[0.2,0.6,1]); 
plot(hrvBioAttn.fftFreqXaxis,hrvBioAttn.fftPowYaxis,'color',[0.5,0.5,0.5]);
leg = {'NCS','NCS 2^{nd} Harm','ECG'};
titl = 'Frequency Analysis';
xlim([0.04,0.7])
plotCute1('Frequency (Hz)','Power (a.u.)',ax2(4),titl,[],0,'Horizontal');

ax2(3) = subplot(nFig,1,3);
stem(rtTime(idxJumpDet),rtData{idxJumpDet,2},'o','color',[0.7,1,0.4]); hold on
stem(rtTime(idxJumpMiss),rtData{idxJumpMiss,2},'+','color',[1,0.4,0.4]);
stem(rtTime(idxNoJumpWrong),rtData{idxNoJumpWrong,2},'s','color',[1,0.4,0.4]);
ylim([0,1500])
leg = {'True Jump','Missed Jump','Wrong'};
plotCute1('Time (s)','Reaction Time (ms)',ax2(3),[],leg,1,'Horizontal');

linkaxes(ax2([1:3]),'x');

%% Saving the results
hrvSaveFormatNcsRelax = [hrvNcsRelax.meanRR, hrvNcsRelax.sdnn, hrvNcsRelax.meanHR,...
    hrvNcsRelax.rmssd, hrvNcsRelax.sdsdRR, hrvNcsRelax.pNN50,hrvNcsRelax.LFpow,...
    hrvNcsRelax.HFpow,hrvNcsRelax.LFHFratio];
hrvSaveFormatBioRelax = [hrvBioRelax.meanRR, hrvBioRelax.sdnn, hrvBioRelax.meanHR,...
    hrvBioRelax.rmssd, hrvBioRelax.sdsdRR, hrvBioRelax.pNN50,hrvBioRelax.LFpow,...
    hrvBioRelax.HFpow,hrvBioRelax.LFHFratio];
hrvSaveFormatNcsAttn = [hrvNcsAttn.meanRR, hrvNcsAttn.sdnn, hrvNcsAttn.meanHR,...
    hrvNcsAttn.rmssd, hrvNcsAttn.sdsdRR, hrvNcsAttn.pNN50,hrvNcsAttn.LFpow,...
    hrvNcsAttn.HFpow,hrvNcsAttn.LFHFratio];
hrvSaveFormatBioAttn = [hrvBioAttn.meanRR, hrvBioAttn.sdnn, hrvBioAttn.meanHR,...
    hrvBioAttn.rmssd, hrvBioAttn.sdsdRR, hrvBioAttn.pNN50,hrvBioAttn.LFpow,...
    hrvBioAttn.HFpow,hrvBioAttn.LFHFratio];
hrvSaveFormatNcsRelax2Harm = [hrvNcsRelax2Harm.meanRR, hrvNcsRelax2Harm.sdnn, hrvNcsRelax2Harm.meanHR,...
    hrvNcsRelax2Harm.rmssd, hrvNcsRelax2Harm.sdsdRR, hrvNcsRelax2Harm.pNN50,hrvNcsRelax2Harm.LFpow,...
    hrvNcsRelax2Harm.HFpow,hrvNcsRelax2Harm.LFHFratio];
hrvSaveFormatNcsAttn2Harm = [hrvNcsAttn2Harm.meanRR, hrvNcsAttn2Harm.sdnn, hrvNcsAttn2Harm.meanHR,...
    hrvNcsAttn2Harm.rmssd, hrvNcsAttn2Harm.sdsdRR, hrvNcsAttn2Harm.pNN50,hrvNcsAttn2Harm.LFpow,...
    hrvNcsAttn2Harm.HFpow,hrvNcsAttn2Harm.LFHFratio];
hrvSaveFormat = [hrvSaveFormatBioRelax;hrvSaveFormatBioAttn;...
                 hrvSaveFormatNcsRelax;hrvSaveFormatNcsAttn;...
                 hrvSaveFormatNcsRelax2Harm;hrvSaveFormatNcsAttn2Harm];

            
         
%savefig(fig, [dataPath,'Analysis\',ncsAttnTestFile,'_RT.fig']);
% save([dataPath,'Analysis\',ncsAttnTestFile,'_RT.mat'], ...
%     'ncsRelaxHrTh','ncsAttnHrThPh',...
%     'hrvNcsRelax','hrvNcsAttn','hrvBioRelax','hrvBioAttn',...
%     'ncsAttn','bioAttn','fs',...
%     'ncsRelax','bioRelax',...
%     'tAttn','tCalib','tRelax',...
%     'opts1','opts2','tStartEndOff','tManualOff',...
%     'signNcs','signNcsCalib','tManualOffCalib',...
%     'ncsCalib','bioCalib','meanDev',...
%     'ncsAttnRespAbd','ncsAttnRespTh','ncsRelaxRespAbd','ncsRelaxRespTh',... 
%     'hrNcsAttn','ncsAttnHrPeak','brNcsAttn','ncsAttnRespPk',...
%     'hrBioAttn','brBioAttnAbd','brBioAttnTh',...
%     'hrNcsRelax','ncsRelaxHrPeak','brNcsRelax','ncsRelaxRespPk',...
%     'hrBioRelax','brBioRelaxAbd','brBioRelaxTh',...
%     'rtData','tManualRToff','rtTime',...
%     'tvAirflow','volAirflow','beltCalibCoeff','vBeltCalib','ncsCalibCoeff','vNcsCalib',...
%     'tCalibStartEndOff',...
%     'vBeltAttn','vNcsAttn','vBeltRelax','vNcsRelax','ncsCalibFile','tCalibStartEndOff');
% 

% save([dataPath,'Analysis\',ncsAttnTestFile,'_AttnTest.mat'], ...
%     'ncsRelaxHrTh','ncsAttnHrThPh',...
%     'hrvNcsRelax','hrvNcsAttn','hrvBioRelax','hrvBioAttn',...
%     'ncsAttn','bioAttn','fs',...
%     'ncsRelax','bioRelax',...
%     'tAttn','tCalib','tRelax',...
%     'opts1','opts2','tStartEndOff','tManualOff',...
%     'signNcs','signNcsCalib','tManualOffCalib',...
%     'rtData','rtDataNew','tManualRToff','rtTime',...
%     'tvAirflow','volAirflow','beltCalibCoeff','vBeltCalib','ncsCalibCoeff','vNcsCalib',...
%     'tCalibStartEndOff',...
%     'ncsCalibFile','tCalibStartEndOff');

save([dataPath,'Analysis\case',num2str(caseNum),'_',ncsFile(12:end),'AttnTest.mat']);
