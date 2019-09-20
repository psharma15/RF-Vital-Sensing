% Reading NCS data using our USRP 2 channel - one near heart, and one on
% the abdomen. Performing some basic analysis of that, including:
% Heartbeat processing, heart rate, HRV etc.
% Respiration waveform processing, breath rate, lung volume, TiTt etc.
% Possible motion information etc.
% Pragya Sharma, ps847@cornell.edu, 08 March 2019

%% ------------------------------------------------------------------------
% Provide input to function.
dataPath = 'C:\Research\NCS\HumanStudyData\Case7\';
ncsFile = '0828_104608Routine2a'; % data at different time instants
%ncsReactionTimeFile = '0613_210551Routine2b_ReactionTime';
bioFile = 'bio_case7_2019-08-28T10_40_21';
bioCalibFile = 'bio_case7_2019-08-28T09_50_43';
fsNcsHigh = 50e3; % NCS sampling rate in Hz
fsBioHigh = 2e3;
tStartEndOff = [5,1]; % Start and end offset wrt NCS data in seconds
tManualOff = 0; % Manual offset between NCS and BIOPAC in seconds
ncsCalibFile = '0828_103057Calib1';
tCalibStartEndOff = [0.5,0.5]; % Start and end offset wrt NCS data in seconds
fsDS = [500, 500]; % [Ncs,Bio] downsampling frequencies
signNcs = [-1 -1 1 -1]; % sign[ampTh phTh ampAbd phAbd]
tManualOffCalib = tManualOff;
signNcsCalib = [-1 -1 1 1];
unwrapPh = [1 1];
unwrapPhCalib = unwrapPh;

%% ------------------------------------------------------------------------
% Define conditions for optional processing here
ifCalib = 1; % If some associated calibration data is needed
ifLoadCalibCoeff = 0; % If pre-saved calibration coefficient are present
ifDownSamp = [1 1]; % If [ncs biopac] data is to be downsampled

%% ------------------------------------------------------------------------
% Reading synchronized data
[ncs,bio,~,~,tNcs,tBio,fig(1)] = ncsBioSync(dataPath,ncsFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tStartEndOff(1),tStartEndOff(2),...
                                 tManualOff,ifDownSamp,fsDS,signNcs,unwrapPh);
                             
[ncsCalib,bioCalib,~,fs,tNcsCalib,tBioCalib,fig(2)] = ...
    ncsBioSync(dataPath,ncsCalibFile,bioCalibFile,fsNcsHigh,...
               fsBioHigh,tCalibStartEndOff(1),tCalibStartEndOff(2),...
               tManualOffCalib,ifDownSamp,fsDS,signNcsCalib,unwrapPhCalib);

                             
%% Optional: Using cross correlation to estimate time shift, in case
% manually is difficult.
tCorr = [5, 12]; 
nStart = tCorr(1)*fs(1); % Assuming same fs for both ncs and biopac
nEnd = tCorr(2)*fs(1); 
[r, lags] = xcorr(bioCalib(nStart:nEnd,3),ncsCalib(nStart:nEnd,3));
figure; plot(lags,r)
[rMax,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDevCalib = lagsMax/fs(1);
fprintf('Suggested NCS calibration time offset is %f\n',tDevCalib);   
 
%% Filtering
% For volume calibration, signals need to be high pass filtered to
% remove baseline: both biopac belts, and volume airflow
opts1.filtType = 'LpHp';
opts1.f3db = 0.05; opts1.orderHP = 8;
opts1.fpLP = 0.6; opts1.fstLP = 1;
bioCalib(:,2:3) = filterLpHp(bioCalib(:,2:3),fs(2),opts1);
bio(:,2:3) = filterLpHp(bio(:,2:3),fs(2),opts1);

opts2.filtType = 'LpHp'; opts2.orderHP = 8;
opts2.f3db = 0.05; opts2.fpLP = 0.6; opts2.fstLP = 1;
ncsCalib(:,3) = filterLpHp(ncsCalib(:,3),fs(1),opts2); % abd amp
ncsCalib(:,1) = filterLpHp(ncsCalib(:,1),fs(1),opts2); % th amp
ncsCalib(:,4) = filterLpHp(ncsCalib(:,4),fs(1),opts2); % abd ph
ncsCalib(:,2) = filterLpHp(ncsCalib(:,2),fs(1),opts2); % th ph

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.05; opts2.fpLP = .6; opts2.fstLP = 1;
ncsRespAbd = filterLpHp(ncs(:,3),fs(1),opts2);
ncsRespTh = filterLpHp(ncs(:,1),fs(1),opts2);
ncsRespAbdPh = filterLpHp(ncs(:,4),fs(1),opts2);
ncsRespThPh = filterLpHp(ncs(:,2),fs(1),opts2);

%% ------------------------------------------------------------------------
% Calibration phase If calibration is selected, we calculate lung volume or 
% tidal volume from BIOPAC airflow data.

if ifCalib == 1
    rng(5)
    % Baseline estimation: Mean deviation from zero when no airflow.
    % Using the routine's data to estimate that.
    tMeanDev = [20, 80]; % Time window in sec to estimate baseline
    figure
    plot(tBio(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2)),bio(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2),4))
    meanDev = mean(bio(tMeanDev(1)*fs(2)+1:tMeanDev(2)*fs(2),4));    
    
    % TV estimated at the expiration end for each breath cycle using airflow. 
    % Cycle starts from an inspiration and ends with expiration.
    [tvAirflow,volAirflow,fig(3)] = bioAirflowTV(bioCalib(:,4),fs(2),meanDev);  

    % TV calibration of Biopac chest belts
    opts1.calibType = 'vol'; opts1.fitEqn = 'BiasedLinear'; 
    opts1.tWin = 2; opts1.minInterceptDist = 0.1;
    volAirflow = filterLpHp(volAirflow,fs(2),opts1);
    
    % Start and stop calibration times: Performing calibration during
    % normal breathing only: ONLY FOR calibType='vol'
    opts1.tCalib = [8,15]; 
    [beltCalibCoeff,vBeltCalib,gof] = bioBeltVolCalib([bioCalib(:,2),bioCalib(:,3)],fs(2),tvAirflow,volAirflow,opts1); 
        
    opts2.tCalib = opts1.tCalib; 
    opts2.calibType = 'vol';     opts2.fitEqn = 'BiasedLinear'; 
    opts2.tWin = 4;              opts2.minInterceptDist = 0.2;
    [ncsCalibCoeff,vNcsCalib] = ncsVolCalib([ncsCalib(:,2),ncsCalib(:,3)],fs(1),tvAirflow,volAirflow,opts2);
       
    % Calibration phase: Plot Biopac and NCS volumes compared to airflow volume
    fig(4) = figure('Position',[400 200 800 600]);
    nFig = 3;
    ax1(1) = subplot(nFig,1,1);
    plot(tBioCalib,bioCalib(:,2))
    hold on; plot(tBioCalib,bioCalib(:,3));         
    plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
    ax1(2) = subplot(nFig,1,2);
    plot(tNcsCalib,ncsCalib(:,1));
    hold on; plot(tNcsCalib,ncsCalib(:,3));
    plotCute1('Time (s)','NCS',ax1(2),[],{'Th','Abd'},1);
    ax1(3) = subplot(nFig,1,3);
    plot(tBioCalib,volAirflow,'color',[79, 209, 98]/256,'LineWidth',2); hold on;
    plot(tBioCalib,vBeltCalib,'color',[15, 36, 193]/256); 
    plot(tNcsCalib,vNcsCalib,'color',[229, 42, 25]/256)
    leg = {'Airflow','Belt','NCS'};
    plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1,'Horizontal');
    linkaxes(ax1,'x')
end

%% ------------------------------------------------------------------------
% Use generated calibration coefficients to estimate Biopac and NCS Vol/TV
opts1.tWin = 4; opts1.minInterceptDist = 0.2;
vBelt = bioBeltVol(bio(:,2:3),beltCalibCoeff,fs(2),opts1);

opts2.tWin = 4;
vNcs = ncsVol([ncsRespThPh,ncsRespAbd],ncsCalibCoeff,fs(1),opts2);

% Plot calibrated Biopac and NCS volumes (both using calib coefficients)
fig(5) = figure('Position',[400 200 600 600]);
nFig = 3;
ax1(1) = subplot(nFig,1,1);
plot(tBio,bio(:,2))
hold on; plot(tBio,bio(:,3));
plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
ax1(2) = subplot(nFig,1,2);
plot(tNcs,ncsRespTh);
hold on; plot(tNcs,ncsRespAbd);
plotCute1('Time (s)','NCS',ax1(2),[],{'Th','Abd'},1);
ax1(3) = subplot(nFig,1,3);
plot(tNcs,vBelt); hold on;
plot(tNcs,vNcs)
leg = {'V Belt','V NCS'};
plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1);
linkaxes(ax1,'x')

%% ------------------------------------------------------------------------
% Correlation and error in Volume estimates
tCompCalib = [5, 22.5]; % Start and end times for comparison. Ignoring coughing etc.
nCompCalib = tCompCalib.*fs(1) + 1; % Assuming same fs.
rmseVolCalibBeltAirflow = sqrt(mean((vBeltCalib(nCompCalib(1):nCompCalib(2)) - ...
    volAirflow(nCompCalib(1):nCompCalib(2))).^2));
rCorrCalibBeltAirflow = corrcoef(vBeltCalib(nCompCalib(1):nCompCalib(2)),...
    volAirflow(nCompCalib(1):nCompCalib(2)));
fprintf('--------------------- Calibration phase (Unit: L) --------------------- \n');
fprintf('rmse Vol(belt,airlflow) = %3.2f, corr Vol(belt,airflow) = %3.2f. \n',...
    rmseVolCalibBeltAirflow,rCorrCalibBeltAirflow(1,2));

rmseVolCalibNcsAirflow = sqrt(mean((volAirflow(nCompCalib(1):nCompCalib(2)) - ...
    vNcsCalib(nCompCalib(1):nCompCalib(2))).^2));
rCorrCalibNcsAirflow = corrcoef(volAirflow(nCompCalib(1):nCompCalib(2)),...
    vNcsCalib(nCompCalib(1):nCompCalib(2)));
fprintf('rmse Vol(ncs,airlflow) = %3.2f, corr Vol(ncs,airflow) = %3.2f. \n',...
    rmseVolCalibNcsAirflow,rCorrCalibNcsAirflow(1,2));

rmseVolCalibBeltNcs = sqrt(mean((vBeltCalib(nCompCalib(1):nCompCalib(2)) - ...
    vNcsCalib(nCompCalib(1):nCompCalib(2))).^2));
rCorrCalibBeltNcs = corrcoef(vBeltCalib(nCompCalib(1):nCompCalib(2)),...
    vNcsCalib(nCompCalib(1):nCompCalib(2)));
fprintf('rmse Vol(ncs,belt) = %3.2f, corr Vol(ncs,belt) = %3.2f. \n',...
    rmseVolCalibBeltNcs,rCorrCalibBeltNcs(1,2));

tComp = [120, 210]; % [20, 310] Start and end times for comparison. Ignoring coughing etc.
nComp = tComp.*fs(1) + 1; % Assuming same fs.
rmseVol = sqrt(mean((vBelt(nComp(1):nComp(2)) - vNcs(nComp(1):nComp(2))).^2));
rCorr = corrcoef(vBelt(nComp(1):nComp(2)),vNcs(nComp(1):nComp(2)));
fprintf('------------------- Post Calibration (t: [%d,%d]s) ------------------- \n',tComp(1),tComp(2));
fprintf('rmse Vol(ncs,belt) = %3.2f, corr Vol(ncs,belt) = %3.2f. \n',...
    rmseVol,rCorr(1,2));

%% ------------------------------------------------------------------------
savefig(fig, [dataPath,'Analysis\',ncsFile,'_linCalibVol.fig']);
save([dataPath,'Analysis\',ncsFile,'_linCalibVol.mat'], 'ncs','bio','fs',...
    'ncsCalib','bioCalib','meanDev','ncsRespAbd','ncsRespTh','tvAirflow',...
    'volAirflow','beltCalibCoeff','vBeltCalib','ncsCalibCoeff','vNcsCalib',...
    'vBelt','vNcs','tBio','tNcs','tBioCalib','tNcsCalib','opts1','opts2',...
    'tStartEndOff','tManualOff','ncsCalibFile','tCalibStartEndOff',...
    'signNcs','signNcsCalib','tManualOffCalib');


