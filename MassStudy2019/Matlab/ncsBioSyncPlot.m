% Reading NCS data using our USRP 2 channel - one near heart, and one on
% the abdomen. Only for synchronization and visualization.
% Pragya Sharma, ps847@cornell.edu, 08 March 2019

%% ------------------------------------------------------------------------
% Provide input to function.
dataPath = 'C:\Research\NCS\HumanStudyData\Case21\';
ncsFile = '0904_153600calib3'; % data at different time instants
%ncsReactionTimeFile = '0613_210551Routine2b_ReactionTime';
bioFile = 'biopac_2019-09-04T15_17_14';
fsNcsHigh = 50e3; % NCS sampling rate in Hz
fsBioHigh = 2e3;
tStartEndOff = [5,1]; % Start and end offset wrt NCS data in seconds
tManualOffNcsBio = 0.14; % Manual offset between NCS and BIOPAC in seconds
fsDS = [500, 500]; % [Ncs,Bio] downsampling frequencies
signNcs = [-1 -1 1 1]; % sign[ampTh phTh ampAbd phAbd]
unwrapPh = [1 1];

%% ------------------------------------------------------------------------
% Define conditions for optional processing here
ifCalib = 0; % If some associated calibration data is needed
ifLoadCalibCoeff = 0; % If pre-saved calibration coefficient are present
ifDownSamp = [1 1]; % If [ncs biopac] data is to be downsampled

%% ------------------------------------------------------------------------
% Reading synchronized data
[ncs,bio,~,fs,tNcs,tBio,fig(1)] = ncsBioSync(dataPath,ncsFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tStartEndOff(1),tStartEndOff(2),...
                                 tManualOffNcsBio,ifDownSamp,fsDS,signNcs,unwrapPh);
                             
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

%% Plot data you need to visualize
opts2.filtType = 'LpHp'; 
opts2.f3db = 0.05; opts2.fpLP = 1.0; opts2.fstLP = 1.5;
ncsRespAbd = filterLpHp(ncs(:,3),fs(1),opts2);
ncsRespTh = filterLpHp(ncs(:,1),fs(1),opts2);
ncsRespAbdPh = filterLpHp(ncs(:,4),fs(1),opts2);
ncsRespThPh = filterLpHp(ncs(:,2),fs(1),opts2);
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
