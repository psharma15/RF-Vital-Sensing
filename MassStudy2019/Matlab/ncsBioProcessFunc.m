% Reading NCS data using our USRP 2 channel - one near heart, and one on
% the abdomen. Performing some basic analysis of that, including:
% Heartbeat processing, heart rate, HRV etc.
% Respiration waveform processing, breath rate, lung volume.
% Derived from ncsBioProcess.m
% Pragya Sharma, ps847@cornell.edu, 24 Sept 2019
 
function ncsBioProcessFunc(caseNum,ncsCalibFile,ncsFile,bioCalibFile, bioFile)
%%
% caseNum = 10;
% ncsCalibFile = '0830_132313calib3';
% ncsFile = '0830_132942Routine3'; % data at different time instants
% bioFile = 'caseNum_2019-08-30T13_03_36';
% bioCalibFile = 'caseNum_2019-08-30T13_03_36';

% [caseOrder,routineName,fileName,~,~] = case7_30Info();
% i=55;
% caseNum = caseOrder(i);
% ncsCalibFile = fileName{i,1};
% ncsFile = fileName{i,2};
% bioCalibFile = fileName{i,3};
% bioFile = fileName{i,4};
nSec = 1;

% set(gcf,'Visible','off');
% set(0,'DefaultFigureVisible','off');

dataPath = ['C:\Research\NCS\HumanStudyData','\Case',num2str(caseNum),'\'];
fsNcsHigh = 50e3; fsBioHigh = 2e3; % Sampling rate in Hz

tOffNcsStEnd = [5,2]; tOffNcsStEndCalib = [2,2]; % Time offset wrt NCS start and end in sec.
tOffNcsThAbd = [0,0]; % [Routine,calib] Manual time offset between NCS Th/Abd in sec. +ve means Abd is shifted to right.
tOffNcsBio = [0, 0]; % [Routine,calib] Manual time offset b/w NCS (Th) & BIOPAC in sec 

signNcs = [1 1 1 1]; signNcsCalib = [1 1 1 1]; % sign[ampTh phTh ampAbd phAbd]

%% ------------------------------------------------------------------------
% Define conditions for optional processing here
ifCalib = 1; % If some associated calibration data is needed
ifDownSamp = [1 1]; % If [ncs biopac] data is to be downsampled
fsDS = [500, 500]; % [Ncs,Bio] downsampling frequencies: Do NOT make them different now, because the other parts of code are not designed to work with that.
if fsDS(1) ~= fsDS(2)
    warning('Please enter same down-sampling frequencies for NCS and BIOPAC or edit sections of code');
    return;
end
unwrapPh = [1 1]; unwrapPhCalib = unwrapPh;

%% ------------------------------------------------------------------------
% Reading synchronized data
[ncs,bio,~,~,tNcs,~,fig(1)] = ncsBioSync(dataPath,ncsFile,bioFile,fsNcsHigh,...
                                 fsBioHigh,tOffNcsStEnd,tOffNcsThAbd(1),...
                                 tOffNcsBio(1),ifDownSamp,fsDS,signNcs,unwrapPh);
                             
[ncsCalib,bioCalib,~,fs,tNcsCalib,~,fig(2)] = ...
    ncsBioSync(dataPath,ncsCalibFile,bioCalibFile,fsNcsHigh,...
               fsBioHigh,tOffNcsStEndCalib,tOffNcsThAbd(2),...
               tOffNcsBio(2),ifDownSamp,fsDS,signNcsCalib,unwrapPhCalib); 

%% ------------------------------------------------------------------------
% Filtering
% Respiration: high pass filtered to remove baseline and low pass to remove
% heartbeat
opts1.filtType = 'LpHp';
opts1.f3db = 0.05; opts1.orderHP = 5;
opts1.fpLP = 0.8; opts1.fstLP = 1.2;
bioCalib(:,2:3) = filterLpHp(bioCalib(:,2:3),fs(2),opts1);
bio(:,2:3) = filterLpHp(bio(:,2:3),fs(2),opts1);

opts2.filtType = 'LpHp'; opts2.orderHP = 5;
opts2.f3db = 0.05; opts2.fpLP = 0.8; opts2.fstLP = 1.2;
ncsCalibRespThAmp = filterLpHp(ncsCalib(:,1),fs(1),opts2); % th amp
ncsCalibRespThPh = filterLpHp(ncsCalib(:,2),fs(1),opts2); % th ph
ncsCalibRespAbdAmp = filterLpHp(ncsCalib(:,3),fs(1),opts2); % abd amp
ncsCalibRespAbdPh = filterLpHp(ncsCalib(:,4),fs(1),opts2); % abd ph

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.05; opts2.fpLP = .8; opts2.fstLP = 1.2;
ncsRespThAmp = filterLpHp(ncs(:,1),fs(1),opts2);
ncsRespThPh = filterLpHp(ncs(:,2),fs(1),opts2);
ncsRespAbdAmp = filterLpHp(ncs(:,3),fs(1),opts2);
ncsRespAbdPh = filterLpHp(ncs(:,4),fs(1),opts2);

close all;

%% ------------------------------------------------------------------------
% Finding correct sign of the NCS waveforms
% *** *** Make sure this time doesn't have motion interference! *** ***
tCorr = [20,40]; % Correlate NCS respiration with Biopac to determine sign
tCorrCalib = [5,10]; 

fprintf('\nCorrecting NCS sign: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));
fprintf('\nCorrecting NCS calib sign: correlating window [%d,%d]s\n',tCorrCalib(1),tCorrCalib(2)); 

bioCorrTh = bio(tCorr(1)*fs(2):tCorr(2)*fs(2),2); 
bioCorrAbd = bio(tCorr(1)*fs(2):tCorr(2)*fs(2),3);
[r1,p1] = corrcoef(bioCorrTh,ncsRespThAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); % Bio Th, NCS Th
r1 = r1(1,2); p1 = p1(1,2); %#ok<*NASGU>
[r2,p2] = corrcoef(bioCorrTh,ncsRespThPh(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r2 = r2(1,2); p2 = p2(1,2);
[r3,p3] = corrcoef(bioCorrAbd,ncsRespAbdAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); % Bio Abd, NCS Abd
r3 = r3(1,2); p3 = p3(1,2);
[r4,p4] = corrcoef(bioCorrAbd,ncsRespAbdPh(tCorr(1)*fs(1):tCorr(2)*fs(1)));
r4 = r4(1,2); p4 = p4(1,2);

[rC1,pC1] = corrcoef(bioCalib(tCorrCalib(1)*fs(2):tCorrCalib(2)*fs(2),2),ncsCalibRespThAmp(tCorrCalib(1)*fs(1):tCorrCalib(2)*fs(1)));
rC1 = rC1(1,2); pC1 = pC1(1,2);
[rC2,pC2] = corrcoef(bioCalib(tCorrCalib(1)*fs(2):tCorrCalib(2)*fs(2),2),ncsCalibRespThPh(tCorrCalib(1)*fs(1):tCorrCalib(2)*fs(1)));
rC2 = rC2(1,2); pC2 = pC2(1,2);
[rC3,pC3] = corrcoef(bioCalib(tCorrCalib(1)*fs(2):tCorrCalib(2)*fs(2),3),ncsCalibRespAbdAmp(tCorrCalib(1)*fs(1):tCorrCalib(2)*fs(1)));
rC3 = rC3(1,2); pC3 = pC3(1,2);
[rC4,pC4] = corrcoef(bioCalib(tCorrCalib(1)*fs(2):tCorrCalib(2)*fs(2),3),ncsCalibRespAbdPh(tCorrCalib(1)*fs(1):tCorrCalib(2)*fs(1)));
rC4 = rC4(1,2); pC4 = pC4(1,2);

fprintf('Updating NCS signs...  ');
signNcs = sign([r1, r2, r3, r4]);
signNcsCalib = sign([rC1, rC2, rC3, rC4]);
    
% Now that we have the signs, updating the waveforms:
ncs = [signNcs(1).*ncs(:,1), signNcs(2).*ncs(:,2), signNcs(3).*ncs(:,3), signNcs(4).*ncs(:,4)];
ncsCalib = [signNcsCalib(1).*ncsCalib(:,1), signNcsCalib(2).*ncsCalib(:,2), signNcsCalib(3).*ncsCalib(:,3), signNcsCalib(4).*ncsCalib(:,4)];
fprintf('Done. \n');

%% ------------------------------------------------------------------------
% % Next step is to determine which NCS thorax and abdomen waveforms to be
% used for breath volume estimation: Based on better correlation in
% specified time duration (should change to overall?)

% Re calculating - you can change the time
tCorr = [2,120]; % Correlate NCS respiration with Biopac to determine sign
fprintf('\nSelecting NCS Amp/Ph for respiration: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));

[r1,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),2),ncsRespThAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r1 = r1(1,2);
[r2,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),2),ncsRespThPh(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r2 = r2(1,2); 
[r3,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),3),ncsRespAbdAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r3 = r3(1,2); 
[r4,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),3),ncsRespAbdPh(tCorr(1)*fs(1):tCorr(2)*fs(1)));
r4 = r4(1,2); 


[~,idxRespCorr] = max(abs([r1,r2,rC1,rC2])); % Find max out of both routine and calibration waveform
thMinCorr = 0.6;
if ((idxRespCorr == 1) && (abs(rC1) >= thMinCorr)) || ((idxRespCorr == 3) && (abs(r1) >= thMinCorr))
% Strict conditions to use amplitude as respiration, default is phase
    ncsRespTh = ncsRespThAmp;
    ncsCalibRespTh = ncsCalibRespThAmp;
    rRespThBio = r1;
    fprintf('Thorax: Using NCS Amplitude waveform for respiration.\n');
    rRespCalibThBio = rC1;

else
    rRespThBio = r2; rRespCalibThBio = rC2;
    ncsRespTh = ncsRespThPh; 
    ncsCalibRespTh = ncsCalibRespThPh;
    fprintf('Thorax: Using NCS Phase waveform for respiration.\n');
end

[~,idxRespCorr] = max(abs([r3,r4,rC3,rC4])); % Find max out of both routine and calibration waveform
thMinCorr = 0.6;
if ((idxRespCorr == 1) && (abs(rC3) >= thMinCorr)) || ((idxRespCorr == 3) && (abs(r3) >= thMinCorr))
% Strict conditions to use amplitude as respiration, default is phase
    rRespAbdBio = r3; rRespCalibAbdBio = rC3;
    ncsRespAbd = ncsRespAbdAmp;
    ncsCalibRespAbd = ncsCalibRespAbdAmp;
    fprintf('Abdomen: Using NCS Amplitude waveform for respiration.\n');
else
    rRespAbdBio = r4; rRespCalibAbdBio = rC4;
    ncsRespAbd = ncsRespAbdPh; 
    ncsCalibRespAbd = ncsCalibRespAbdPh;
    fprintf('Abdomen: Using NCS Phase waveform for respiration.\n');
end

%% Optional: Using cross correlation to estimate time shift, in case
% manually is difficult. But this is not very accurate
tCorr = [5, 12]; 
fprintf('\nTime sync, Calibration: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));

nStart = tCorr(1)*fs(1); % Assuming same fs for both ncs and biopac
nEnd = tCorr(2)*fs(1); 
maxLag = 1000;
[r, lags] = xcorr(bioCalib(nStart:nEnd,2),ncsCalibRespTh(nStart:nEnd),maxLag);
figure; plot(lags,r)
[~,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDevCalib = lagsMax/fs(1);
fprintf('Suggested NCS Th-Bio Th calibration time offset is %f\n',tDevCalib);   

[r, lags] = xcorr(ncsCalibRespTh(nStart:nEnd),ncsCalibRespAbd(nStart:nEnd),maxLag);
figure; plot(lags,r)
[~,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDevCalibThAbd = lagsMax/fs(1);
fprintf('Suggested NCS Th-Abd calibration time offset is %f\n',tDevCalibThAbd); 

tCorr = [2, 100]; 
fprintf('\nTime sync, Routine: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));

nStart = tCorr(1)*fs(1); % Assuming same fs for both ncs and biopac
nEnd = tCorr(2)*fs(1); 
[r, lags] = xcorr(bio(nStart:nEnd,2),ncsRespTh(nStart:nEnd),maxLag);
figure; plot(lags,r)
[~,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDev = lagsMax/fs(1);
fprintf('Suggested NCS Th-Bio Th time offset is %f\n',tDev);  

[r, lags] = xcorr(ncsRespTh(nStart:nEnd),ncsRespAbd(nStart:nEnd),maxLag);
figure; plot(lags,r)
[~,rMaxIdx] = max(abs(r));
lagsMax = lags(rMaxIdx);
tDevThAbd = lagsMax/fs(1);
fprintf('Suggested NCS Th-Abd time offset is %f\n',tDevThAbd);  

%% ONLY call this block once!!
% Now shift the waveforms to compensate the time difference
% Only shift if time difference is more than a threshold and less than a
% threshold: Ideally minimum should be sampling frequency

fprintf('Performing synchronization based on time-shift estimate...  \n ');
thMinCorr = 0.4;
tOffMinMax = [0.005, 1.0];
% First the calibration waveforms (Bio, NCS Th):
if (abs(tDevCalib)>=tOffMinMax(1)) && (abs(tDevCalib) <= tOffMinMax(2)) && (abs(rRespCalibThBio) > thMinCorr)
    nSampDev = abs(tDevCalib) * fs(1); % same BIO and NCS frequencies
    if tDevCalib > 0
        bioCalib = bioCalib(nSampDev+1:end,:);
        ncsCalib = ncsCalib(1:end-nSampDev,:);
    else
        bioCalib = bioCalib(1:end-nSampDev,:);
        ncsCalib = ncsCalib(nSampDev+1:end,:);
    end
    tNcsCalib = tNcsCalib(1:end-nSampDev);
    fprintf('\nSync: NCS calib Th with Biopac.\n');
else
    tDevCalib = 0; % To not be recorded, since we are not correcting
end

% Calibration waveforms (NCS Th, NCS Abd):
if (abs(tDevCalibThAbd)>=tOffMinMax(1)) && (abs(tDevCalibThAbd) <= tOffMinMax(2)) && (abs(rRespCalibAbdBio) > thMinCorr) && (abs(rRespCalibThBio) > thMinCorr)
    nSampDev = abs(tDevCalibThAbd) * fs(1); % same BIO and NCS frequencies
    temp = ncsCalib; % Temporary variable
    ncsCalib = zeros(size(ncsCalib,1)-nSampDev,4);
    if tDevCalibThAbd > 0
        ncsCalib(:,1:2) = temp(nSampDev+1:end,1:2);
        ncsCalib(:,3:4) = temp(1:end-nSampDev,3:4);
        bioCalib = bioCalib(nSampDev+1:end,:);
    else
        ncsCalib(:,1:2) = temp(1:end-nSampDev,1:2);
        ncsCalib(:,3:4) = temp(nSampDev+1:end,3:4);
        bioCalib = bioCalib(1:end-nSampDev,:);
    end
    tNcsCalib = tNcsCalib(1:end-nSampDev);
    fprintf('\nSync: NCS calib Th with NCS calib Abd.\n');
else
    tDevCalibThAbd = 0; % To not record as a shift
end

tBioCalib = tNcsCalib; % Make sure sampling frequencies are set same

% Waveforms from routine (Bio, NCS Th):
if (abs(tDev)>=tOffMinMax(1)) && (abs(tDev) <= tOffMinMax(2)) && (abs(rRespThBio)> thMinCorr) 
    nSampDev = abs(tDev) * fs(1); % same BIO and NCS frequencies
    if tDev > 0
        bio = bio(nSampDev+1:end,:);
        ncs = ncs(1:end-nSampDev,:);
    else
        bio = bio(1:end-nSampDev,:);
        ncs = ncs(nSampDev+1:end,:);
    end
    tNcs = tNcs(1:end-nSampDev);
    fprintf('Sync: NCS Th with Biopac.\n');
else
    tDev = 0;
end

%  waveforms (NCS Th, NCS Abd):
if (abs(tDevThAbd)>=tOffMinMax(1)) && (abs(tDevThAbd) <= tOffMinMax(2)) && (abs(rRespThBio)>thMinCorr) && (abs(rRespAbdBio)>thMinCorr)
    nSampDev = abs(tDevThAbd) * fs(1); % same BIO and NCS frequencies
    temp = ncs; % Temporary variable
    ncs = zeros(size(ncs,1)-nSampDev,4);
    if tDevThAbd > 0
        ncs(:,1:2) = temp(nSampDev+1:end,1:2);
        ncs(:,3:4) = temp(1:end-nSampDev,3:4);
        bio = bio(nSampDev+1:end,:);
    else
        ncs(:,1:2) = temp(1:end-nSampDev,1:2);
        ncs(:,3:4) = temp(nSampDev+1:end,3:4);
        bio = bio(1:end-nSampDev,:);
    end
    tNcs = tNcs(1:end-nSampDev);
    fprintf('Sync: NCS Th with NCS Abd.\n');
else
    tDevThAbd = 0;
end
   
tBio = tNcs;

fprintf('Done \n');

%%
% We need filtering again for respiration/ should have truncated filtered
% waveforms above.
opts2.filtType = 'LpHp'; opts2.orderHP = 5;
opts2.f3db = 0.05; opts2.fpLP = 0.8; opts2.fstLP = 1.2;
ncsCalibRespThAmp = filterLpHp(ncsCalib(:,1),fs(1),opts2); % th amp
ncsCalibRespThPh = filterLpHp(ncsCalib(:,2),fs(1),opts2); % th ph
ncsCalibRespAbdAmp = filterLpHp(ncsCalib(:,3),fs(1),opts2); % abd amp
ncsCalibRespAbdPh = filterLpHp(ncsCalib(:,4),fs(1),opts2); % abd ph

opts2.filtType = 'LpHp'; 
opts2.f3db = 0.05; opts2.fpLP = .8; opts2.fstLP = 1.2;
ncsRespThAmp = filterLpHp(ncs(:,1),fs(1),opts2);
ncsRespThPh = filterLpHp(ncs(:,2),fs(1),opts2);
ncsRespAbdAmp = filterLpHp(ncs(:,3),fs(1),opts2);
ncsRespAbdPh = filterLpHp(ncs(:,4),fs(1),opts2);

% Filtering heartbeat
opts0.filtType = 'LpHp';
opts0.orderHP = 5; opts0.Ast = 20;
opts0.f3db = 0.7; opts0.fpLP = 1.5; opts0.fstLP = 1.8;
ncsHeartThAmp = filterLpHp(ncs(:,1),fs(1),opts0);
ncsHeartThPh = filterLpHp(ncs(:,2),fs(1),opts0);
ncsHeartAbdAmp = filterLpHp(ncs(:,3),fs(1),opts0);
ncsHeartAbdPh = filterLpHp(ncs(:,3),fs(1),opts0);


bioHeartBP = filterLpHp(bio(:,1),fs(2),opts0);
opts0.filtType = 'Hp';
opts0.f3db = 0.1; % Only removing baseline
bioHeart = filterLpHp(bio(:,1),fs(2),opts0);

figure
plot(tNcs,ncsHeartThAmp.*100,tNcs,bio(:,1));

close all;

%% ------------------------------------------------------------------------
% Again determining which waveform to select for respiration in thorax and
% abdomen, mostly should remain the same as earlier (before
% synchronization), but repeating it.
% Next step is to determine which NCS thorax and abdomen waveforms to be
% used for breath volume estimation: Based on better correlation in
% specified time duration (should change to overall?)

% Re calculating - you can change the time
tCorr = [2,120]; % Correlate NCS respiration with Biopac
fprintf('\nRepeat selection respiration NCS Amp/Ph after sync: \ncorrelating window [%d,%d]s\n',tCorr(1),tCorr(2));

[r1,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),2),ncsRespThAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r1 = r1(1,2);
[r2,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),2),ncsRespThPh(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r2 = r2(1,2); 
[r3,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),3),ncsRespAbdAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r3 = r3(1,2); 
[r4,~] = corrcoef(bio(tCorr(1)*fs(2):tCorr(2)*fs(2),3),ncsRespAbdPh(tCorr(1)*fs(1):tCorr(2)*fs(1)));
r4 = r4(1,2); 

[~,idxRespCorr] = max(abs([r1,r2,rC1,rC2]));
thMinCorr = 0.6;
if ((idxRespCorr == 1) && (abs(rC1) >= thMinCorr)) || ((idxRespCorr == 3) && (abs(r1) >= thMinCorr))
% Strict conditions to use amplitude as respiration, default is phase
    rRespThBio = r1;    rRespCalibThBio = rC1;
    ncsRespTh = ncsRespThAmp;
    ncsCalibRespTh = ncsCalibRespThAmp;
    fprintf('Thorax Resp: Amp, r: %3.2f.\n',rRespThBio);

else
    rRespThBio = r2; rRespCalibThBio = rC2; 
    ncsRespTh = ncsRespThPh; 
    ncsCalibRespTh = ncsCalibRespThPh;
    fprintf('Thorax Resp: Ph, r: %3.2f.\n',rRespThBio);
end

[~,idxRespCorr] = max(abs([r3,r4,rC3,rC4]));
thMinCorr = 0.6;
if ((idxRespCorr == 1) && (abs(rC3) >= thMinCorr)) || ((idxRespCorr == 3) && (abs(r3) >= thMinCorr))
% Strict conditions to use amplitude as respiration, default is phase
    rRespAbdBio = r3; rRespCalibAbdBio = rC3; 
    ncsRespAbd = ncsRespAbdAmp;
    ncsCalibRespAbd = ncsCalibRespAbdAmp;
    fprintf('Abdomen Resp: Amp, r: %3.2f.\n',rRespAbdBio);
else
    rRespAbdBio = r4; rRespCalibAbdBio = rC4; 
    ncsRespAbd = ncsRespAbdPh; 
    ncsCalibRespAbd = ncsCalibRespAbdPh;
    fprintf('Abdomen Resp: Ph, r: %3.2f.\n',rRespAbdBio);
end

%% ------------------------------------------------------------------------
% Determining which NCS waveform to be used for heartbeat: Amp/Ph
% Thorax has usually more heartbeat than abdomen, but we are making sure
% Correlate NCS heart with bandpassed ECG (Biopac) to determine

tCorr = [5,120]; 
fprintf('\nSelecting heart NCS Amp/Ph: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));

[r1,~] = corrcoef(bioHeartBP(tCorr(1)*fs(2):tCorr(2)*fs(2),1),ncsHeartThAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r1 = r1(1,2);
[r2,~] = corrcoef(bioHeartBP(tCorr(1)*fs(2):tCorr(2)*fs(2),1),ncsHeartThPh(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r2 = r2(1,2); 
[r3,~] = corrcoef(bioHeartBP(tCorr(1)*fs(2):tCorr(2)*fs(2),1),ncsHeartAbdAmp(tCorr(1)*fs(1):tCorr(2)*fs(1))); 
r3 = r3(1,2); 
[r4,~] = corrcoef(bioHeartBP(tCorr(1)*fs(2):tCorr(2)*fs(2),1),ncsHeartAbdPh(tCorr(1)*fs(1):tCorr(2)*fs(1)));
r4 = r4(1,2); 

rVecHeart = ([r1,r2,r3,r4]);
[~,idxHeartCorr] = max(abs(rVecHeart));
thMinCorr = 0.6;
% Default selection is NCS Th Amplitude and only change from here if other
% waveforms have not only higher, but more correlation than minimum
ncsHeart = ncsHeartThAmp; 
idxHeart = 1;
if ((idxHeartCorr == 2) && (abs(rVecHeart(2)) >= thMinCorr))
    idxHeart = idxHeartCorr;
    ncsHeart = ncsHeartThPh; corrFacHeart = rVecHeart(2);
    fprintf('Heart: NCS Th Ph, r: %3.2f.\n',rVecHeart(2));
elseif ((idxHeartCorr == 3) && (abs(rVecHeart(3)) >= thMinCorr))
    idxHeart = idxHeartCorr;
    ncsHeart = ncsHeartAbdAmp; corrFacHeart = rVecHeart(3);
    fprintf('Heart: NCS Abd Amp, r: %3.2f.\n',rVecHeart(3));
elseif ((idxHeartCorr == 4) && (abs(rVecHeart(4)) >= thMinCorr))
    idxHeart = idxHeartCorr;
    ncsHeart = ncsHeartAbdPh; corrFacHeart = rVecHeart(4);
    fprintf('Heart: NCS Abd Ph, r: %3.2f.\n',rVecHeart(4));
else
    idxHeart = 1;
    corrFacHeart = rVecHeart(1);
    fprintf('Heart: NCS Th Amp, r: %3.2f.\n',rVecHeart(1));
end
signNcsHeart = 1;
if (sign(corrFacHeart)) < 0
    signNcsHeart = -1;
    ncsHeart = ncsHeart .* signNcsHeart; % Changing the sign
end
    
ncsHeartFiltered = ncsHeart;

%% Now get heartbeat using wavelet

% Save filtered waveform differently, so as to not over-write or change
% remaining code.
    
% Get the wavelet component out of heart and see if it correlates well:
% Changing sampling rate, so that it is easier to separate heart from
% respiration for range 50-100 BPM
% With sample rate 348 Hz, 7th and 8th components will be having
% freq: [1.36 - 2.72] Hz & [0.68 - 1.36] Hz = 40-163 bpm
% Sample rate: 400 Hz, 46.8-187 BPM 
% 420 Hz: [0.82,1.64][1.64,3.28] Hz, [49.2,196.8] BPM

fsWav = 400; 
wName = 'db10';
ncsDS = resample(ncs(:,idxHeart),tNcs,fsWav);
tNcsDS = 0:1/fsWav:(length(ncsDS)-1)/fsWav;
figure
plot(tNcs,ncsHeart);hold on;plot(tNcsDS,ncsDS);

[ncsRecL7,~,~] = wavedecrec1(ncsDS,fsWav,wName,7,1);
[ncsRecL8,~,~] = wavedecrec1(ncsDS,fsWav,wName,8,1);

ncsHeartWavDS = ncsRecL7 + ncsRecL8;
ncsHeartWavDS = signNcsHeart.*ncsHeartWavDS;

ncsHeartWav = resample(ncsHeartWavDS, fs(1),fsWav); % Upsample back
ncsHeartWav = ncsHeartWav(1:length(ncsHeartFiltered));

opts2.filtType = 'Lp';
opts2.fpLP = 1.8; opts2.fstLP = 2.2;
ncsHeart = filterLpHp(ncsHeartWav,fs(1),opts2); % th amp

figure
plot(tNcs,ncsHeartFiltered)
hold on
plot(tNcs,ncsHeartWav)
plot(tNcs,ncsHeart)
legend('Heart Filtered','Heart Wavelet','Heart Wavelet Filtered')

% ncsHeart = ncsHeartWav;

%% Assigning data quality (respiration), once the waveforms are selected
% First divide the time into epochs of delT = 5s
qualEpochResp = 5; %Time in seconds
nSampEpoch = qualEpochResp*fs(1);
numEpoch = ceil((tNcs(end) - tNcs(1))/qualEpochResp);
qualIdxResp = ones(numEpoch,1);
tQualIdxResp = zeros(numEpoch,1); % Specifies start time of each epoch
qualRespConf = ones(numEpoch,3);
for i = 0:numEpoch-1
    nStart = i*nSampEpoch + 1;
    nEnd = (i+1)*nSampEpoch;
    tQualIdxResp(i+1) = tNcs(nStart);
    if nEnd > length(bio(:,3))
        nEnd = length(bio(:,3));
    end
    ncsRespThEpoch = ncsRespTh(nStart:nEnd);
    ncsRespAbdEpoch = ncsRespAbd(nStart:nEnd);
    bioRespThEpoch = bio(nStart:nEnd,2);
    bioRespAbdEpoch = bio(nStart:nEnd,3);
    
    minCorr = 0.4;
    % ---------------------------------------------------------------------
    % If Biopac Thorax and abdomen are uncorrelated -> low Quality
    % Reasons include motion, thoracoabdominal asynchrony
    [r1,p] = corrcoef(bioRespThEpoch,bioRespAbdEpoch);
    if (r1(1,2)<minCorr) && (p(1,2)<0.05) 
        % If correlation is low with high probability, low quality
        qualIdxResp(i+1) = 0; 
    end
    
    % ---------------------------------------------------------------------
    % If NCS abdomen is not correlated well to Biopac -> low quality
    % Reasons include motion, poor deep breathing signal
    % ----------------------------------
    % If NCS thorax is not correlated well to Biopac when abdomen is poor
    % -> low quality
    % Reasons include motion, poor deep breathing signal
    [r2,p2] = corrcoef(bioRespAbdEpoch,ncsRespAbdEpoch);
    [r3,p3] = corrcoef(bioRespThEpoch,ncsRespThEpoch);
    
    if (r2(1,2)<minCorr) && (p2(1,2)<0.05)
        qualIdxResp(i+1) = 0;
    end

    if (r3(1,2)<minCorr) && (p3(1,2)<0.05)
        qualIdxResp(i+1) = 0;
    end
    qualRespConf(i+1,:) = [r1(1,2),r2(1,2),r3(1,2)];
end

figure
ax(1) = subplot(3,1,1);
plot(tNcs,ncsRespTh); ylabel('NCS Th')
ax(2) = subplot(3,1,2);
plot(tNcs,ncsRespAbd); ylabel('NCS Abd')
ax(3) = subplot(3,1,3);
stairs(tQualIdxResp,qualIdxResp);
hold on
plot(tQualIdxResp, qualRespConf(:,1),tQualIdxResp, qualRespConf(:,2),tQualIdxResp, qualRespConf(:,3));
legend('Quality Index','r1','r2','r3')
ylabel('Resp data quality')
ylim([-0.25,1.25]);
linkaxes(ax,'x');

fprintf('Estimated NCS respiration quality index.\n');

%% Assigning data quality (heartbeat)
% First divide the time into epochs of 4 s
qualEpochHeart = 4; %Time in seconds
nSampEpoch = qualEpochHeart*fs(1);
numEpoch = ceil((tNcs(end) - tNcs(1))/qualEpochHeart);
tQualIdxHeart = zeros(numEpoch,1); % Specifies start time of each epoch

opts2.filtType = 'LpHp'; 
opts2.fpLP = 10; opts2.fstLP = 11;
opts2.f3db = 0.7;
ncsHeartHf = filterLpHp(ncs(:,idxHeart),fs(1),opts2); % Heartbeat with more high frequency content

% Implementing a similar motion detection algorithm as earlier:
% - for each epoch, instead of beat-by-beat basis, as this is faster.
% - No training time - it is exactly one class classification by fitting
% entire data & detecting outliers
% - features: [perPsdPower, relRmsPow, meanEp, varEp]
numFeat = 3; % Number of features
tTrain = [20,60]; % Selecting prior training time for heartbeat quality detection
epTrain = int64(tTrain./qualEpochHeart);
epTrain(2) = epTrain(2) - 1;
wind = hann(nSampEpoch); % Important to give less weightage to boundary points
X = zeros(numEpoch,numFeat); % Feature vector

rmsPowLastEp = 1; % Initialization for first epoch

for i = 0:numEpoch-1
    nStart = i*nSampEpoch + 1;
    nEnd = (i+1)*nSampEpoch;
    tQualIdxHeart(i+1) = tNcs(nStart);
    if nEnd > length(bio(:,1))
        nEnd = length(bio(:,1));
        wind = hann(nEnd-nStart+1);
    end
    ncsHeartHfEpoch = ncsHeartHf(nStart:nEnd);
    ncsHeartEpoch = ncsHeart(nStart:nEnd);
    
    % ---------------------------------------------------------------------
    % Feature estimation
    % normPsdPow
    [pxx,freqPxx] = periodogram(ncsHeartHfEpoch,wind,[],fs(1));
    pBand = bandpower(pxx,freqPxx,[0.7,5],'psd');
    pTot = bandpower(pxx,freqPxx,'psd');
    perPsdPow = pBand/pTot;
    
    % relRMSPow
    rmsPowEp = rms(ncsHeartHfEpoch);
    relRmsPow = rmsPowEp/rmsPowLastEp;
    rmsPowLastEp = rmsPowEp; % Only update after it is used once
    
    % meanEp 
    meanEp = mean(ncsHeartEpoch); % Not very indicative
    
    % relVarEp 
    varEp = std(ncsHeartEpoch).^2; % Slightly indicative but fails when waveform power becomes less - then variance is less
    
    % kurtosis: adding this, as should better indicate change in
    % probability distribution tails as we have some outlier values in case
    % of motion
    kurtEp = kurtosis(ncsHeartHfEpoch);
    
    % Update the feature vector:
    X(i+1,:) = [perPsdPow, relRmsPow, kurtEp]; % 
end
% Rel RMS power not defined for the 1st epoch, so substituting by 2nd epoch
X(1,2) = X(2,2); 
% Normalizing feature vector with only the training data mean and standard
% deviation information. If no training, then the entire.
Xnorm = (X - mean(X))./std(X); 

%% -------------------------------------------------------------------------
% fitcsvm for one-class SVM classification
% I'm not training at all
rng(1);
y = ones(numEpoch,1); % 1: No motion, 0: Motion
modelHeartMotion = fitcsvm(Xnorm,...
    y,'KernelFunction','rbf','KernelScale','auto','Nu',0.5,'OutlierFraction',0.1);
% Testing
crossValModel = crossval(modelHeartMotion,'kfold',5);
kFoldLossPred = kfoldLoss(crossValModel);

[~,scorePredOrig] = predict(modelHeartMotion,Xnorm);
scorePred = tanh(scorePredOrig); % Normalizing the scores

% Automated threshold variation based on scores in good period: less than
% 50% of the mean score is classified as motion - only less than is used,
% assuming score is high in the training period (not affected by motion)
scorePredCalib = mean(scorePred(epTrain));
threshPred = 0.2*scorePredCalib; 
qualIdxHeart = ones(numEpoch,1); % 1: Good quality, no motion
qualIdxHeart(scorePred<threshPred) = 0; % 0: If motion is detected

figure('position',[100,100,750,600])
ax = [];
nFig = 3;
ax(1) = subplot(nFig,1,1);
plot(tNcs,ncsHeart);ylabel('NCS Heartbeat (a.u.)')
ax(2) = subplot(nFig,1,2);
stairs(tQualIdxHeart,Xnorm(:,1));
hold on
stairs(tQualIdxHeart,Xnorm(:,2));
stairs(tQualIdxHeart,Xnorm(:,3));
% stairs(tQualIdxHeart,Xnorm(:,4));
legend({'Psd','relRMS','kurtEp'});
ylabel('Features');
ax(3) = subplot(nFig,1,3);
stairs(tQualIdxHeart,scorePred);hold on
stairs(tQualIdxHeart,qualIdxHeart);
stairs(tQualIdxHeart,threshPred.*ones(length(qualIdxHeart),1))
ylabel('Quality Index (heart)'); xlabel('Time (s)')
text (tQualIdxHeart(end-10),-0.11,'0: Poor Quality')
% ylim([-0.25,1.25])
linkaxes(ax,'x');

fprintf('Estimated NCS heartbeat quality index.\n');

%% ------------------------------------------------------------------------
% Calibration phase If calibration is selected, we calculate lung volume or 
% tidal volume from BIOPAC airflow data.

if ifCalib == 1
    fprintf('Performing volume calibration... \n');
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
%     volAirflow = volAirflow - .1;

    % Start and stop calibration times: Performing calibration during
    % normal breathing only: ONLY FOR calibType='vol'
    opts1.tCalib = [9,17]; 
    [beltCalibCoeff,vBeltCalib,~] = bioBeltVolCalib([bioCalib(:,2),bioCalib(:,3)],fs(2),tvAirflow,volAirflow,opts1); 
        
    opts2.tCalib = opts1.tCalib; 
    opts2.calibType = 'vol';     opts2.fitEqn = 'BiasedLinear'; 
    opts2.tWin = 4;              opts2.minInterceptDist = 0.2;
    [ncsCalibCoeff,vNcsCalib] = ncsVolCalib([ncsCalibRespTh, ncsCalibRespAbd],fs(1),tvAirflow,volAirflow,opts2);
       
    % Calibration phase: Plot Biopac and NCS volumes compared to airflow volume
    clear ax1
    fig(4) = figure('Position',[400 200 800 600]);
    nFig = 3;
    tBioCalib = tNcsCalib;
    ax1(1) = subplot(nFig,1,1);
    plot(tBioCalib,bioCalib(:,2).*beltCalibCoeff(1))
    hold on; plot(tBioCalib,bioCalib(:,3).*beltCalibCoeff(2));         
    plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
    
    ax1(2) = subplot(nFig,1,2);
    plot(tNcsCalib,ncsCalibRespTh.*ncsCalibCoeff(1));
    hold on; plot(tNcsCalib,ncsCalibRespAbd.*ncsCalibCoeff(2));
    plotCute1('Time (s)','NCS',ax1(2),[],{'Th','Abd'},1);
    
    ax1(3) = subplot(nFig,1,3);
    plot(tBioCalib,volAirflow,'color',[79, 209, 98]/256,'LineWidth',2); hold on;
    plot(tBioCalib,vBeltCalib,'color',[15, 36, 193]/256); 
    plot(tNcsCalib,vNcsCalib,'color',[229, 42, 25]/256)
    leg = {'Airflow','Belt','NCS'};
    plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1,'Horizontal');
    linkaxes(ax1,'x')
    
    fprintf('Vol calibration Done.\n');
end

%%
pause(nSec);

close all

%% ------------------------------------------------------------------------
% Use generated calibration coefficients to estimate Biopac and NCS Vol/TV

fprintf('Performing Volume estimation... ');

opts1.tWin = 4; opts1.minInterceptDist = 0.2;
vBelt = bioBeltVol(bio(:,2:3),beltCalibCoeff,fs(2),opts1);

opts2.tWin = 4;
vNcs = ncsVol([ncsRespTh,ncsRespAbd],ncsCalibCoeff,fs(1),opts2);

fprintf('Done.\n');

% Plot calibrated Biopac and NCS volumes (both using calib coefficients)
clear ax1
fig(5) = figure('Position',[400 200 600 600]);
nFig = 4;
ax1(1) = subplot(nFig,1,1);
plot(tBio,bio(:,2).*beltCalibCoeff(1))
hold on; plot(tBio,bio(:,3).*beltCalibCoeff(2));
plotCute1('Time (s)','Bio Belt (V)',ax1(1),[],{'Bio Th Belt','Bio Abd Belt'},1);
ax1(2) = subplot(nFig,1,2);
plot(tNcs,ncsRespTh.*ncsCalibCoeff(1));
hold on; plot(tNcs,ncsRespAbd.*ncsCalibCoeff(2));
plotCute1('Time (s)','NCS',ax1(2),[],{'Th','Abd'},1);
ax1(3) = subplot(nFig,1,3);
plot(tNcs,vBelt); hold on;
plot(tNcs,vNcs)
respConf = qualRespConf(:,1).*qualRespConf(:,2).*qualRespConf(:,3);
plot(tQualIdxResp,respConf,'color',[0.5,.9,0.5],'LineWidth',1)
leg = {'V Belt','V NCS','Confidence'};
plotCute1('Time (s)','Volume (L)',ax1(3),[],leg,1);

ax1(4) = subplot(nFig,1,4);
qualRespFs = repmat((qualIdxResp(:))',qualEpochResp*fs(1),1); % Making sure qualIdxResp is initially a column, converted to a row vector
qualRespFs = qualRespFs(:);
qualRespFs = qualRespFs(1:length(vBelt)); % As the last quality index may be for less than 5 whole seconds
vBeltCorr = vBelt.*qualRespFs;
vNcsCorr = vNcs.*qualRespFs;

plot(tNcs,vBeltCorr);hold on;
plot(tNcs,vNcsCorr)
respConf = qualRespConf(:,1).*qualRespConf(:,2).*qualRespConf(:,3);
plot(tQualIdxResp,respConf,'color',[0.5,.9,0.5],'LineWidth',1)
leg = {'V Belt Corr','V NCS Corr','Confidence'};
plotCute1('Time (s)','Volume Corrected (L)',ax1(4),[],leg,1);

linkaxes(ax1,'x')

%%
pause(nSec);
close all;
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


vBelt = vBelt(:);vNcs = vNcs(:);
    
rmseVol = sqrt(mean((vBelt - vNcs).^2));
rCorr = corrcoef(vBelt,vNcs);
fprintf('------------------- Post Calibration ------------------- \n');
fprintf('rmse Vol(ncs,belt) = %3.2f, corr Vol(ncs,belt) = %3.2f. \n',...
    rmseVol,rCorr(1,2));

idxVolNonZero = vBeltCorr ~=0;
vBeltCorr = vBeltCorr(idxVolNonZero);
vNcsCorr = vNcsCorr(idxVolNonZero);
rmseVolCorrect = sqrt(mean((vBeltCorr - vNcsCorr).^2));
rCorrCorrect = corrcoef(vBeltCorr,vNcsCorr);
fprintf('rmse Vol(ncsCorrect,beltCorrect) = %3.2f, corr Vol(ncsCorrect,beltCorrect) = %3.2f. \n',...
    rmseVolCorrect,rCorrCorrect(1,2));

%% Respiration BR estimation

fprintf('Performing BR estimation... \n');
ncsResp = ncsRespAbd;

opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.3;

[brNcs,ncsRespPk] = brEst(ncsResp,fs(1),opts3);
pkMaxNcs = ncsRespPk.idx(ncsRespPk.ind == 1);
pkMinNcs = ncsRespPk.idx(ncsRespPk.ind == 0);
pkMaxNcs = pkMaxNcs(ncsRespPk.idxValidPk);
pkMinNcs = pkMinNcs(ncsRespPk.idxValidPk);

bioResp = bio(:,2);
[brBio,bioRespPk] = brEst(bioResp,fs(1),opts3);
pkMaxBio = bioRespPk.idx(bioRespPk.ind == 1);
pkMinBio = bioRespPk.idx(bioRespPk.ind == 0);
pkMaxBio = pkMaxBio(bioRespPk.idxValidPk);
pkMinBio = pkMinBio(bioRespPk.idxValidPk);

rmseBR = sqrt(mean((brBio-brNcs).^2));
fprintf('Done. RMSE BR: %3.2f\n',rmseBR);

%Peforming correction
% --------------------------------------------------------------------------------
% For breath rate estimation, each time instant depends on past 15
% sec of data. So poor data quality should result in wrong
% estimation for a larger window
qualRespBR = qualIdxResp;
counter = 1;
while counter <= length(qualRespBR)
    if qualRespBR(counter) == 0
        % So our quality epochs are 5 s, and BR estimation depends
        % on past 15 s, so a poor quality for 5 s epoch affects
        % BR estimation at next two epochs as well, so those two
        % epochs are also rejected for BR consideration.
        % NOT Ideal?
        if counter + 2 <= length(qualRespBR)
            qualRespBR(counter:counter+2) = 0; 
        else
            qualRespBR(counter:end) = 0; % Boundary: There is only one epoch after second last
        end
        counter = counter + 3;
    else
        counter = counter + 1;
    end
end
        
qualRespBRfs = repmat(qualRespBR(:)',[qualEpochResp*fs,1]); % Making sure qualIdxResp is initially a column, converted to a row vector
qualRespBRfs = qualRespBRfs(:);
qualRespBRfs = qualRespBRfs(1:size(ncsRespAbd,1));
rmseBRcorr = sqrt(mean(((brBio-brNcs).*qualRespBRfs).^2));
fprintf('RMSE BR after correction: %3.2f\n',rmseBRcorr);

figure('Position',[100,100,600,600]);
nFig = 3;
ax1 = [];
ax1(1) = subplot(nFig,1,1);
plot(tNcs,ncsResp);hold on
plot(tNcs(pkMaxNcs),ncsResp(pkMaxNcs),'o');
plot(tNcs(pkMinNcs),ncsResp(pkMinNcs),'*');

ax1(2) = subplot(nFig,1,2);
plot(tNcs,bioResp); hold on
plot(tNcs(pkMaxBio), bioResp(pkMaxBio),'o');
plot(tNcs(pkMinBio),bioResp(pkMinBio),'*');

ax1(3) = subplot(nFig,1,3);
plot(tNcs,brNcs);
hold on
plot(tNcs,brBio);

linkaxes(ax1,'x')

%% 

pause(nSec);
close all;

%% Estimating Heart Rate
% opts2.filtType = 'LpHp'; 
% opts2.orderHP = 5; opts2.Ast = 20;
% opts2.f3db = 0.7; opts2.fpLP = 1.5; opts2.fstLP = 1.8;

fprintf('Performing HR estimation... \n');

% Filtering ECG
opts2.filtType = 'LpHp';
opts2.f3db = 4; opts2.fpLP = 20; opts2.fstLP = 25;
ecgFilt = filterLpHp(bio(:,1),fs(2),opts2);

opts3.tWinHR = 4; % Same for ECG and NCS
opts3.tWin = 0.5;
opts3.minInterceptDist = 0.05;
% [hrNcs, ncsHeartPk] = hrEst(ncsHeart,fs(1),opts3);
% pkMaxNcs = ncsHeartPk.idx(ncsHeartPk.ind == 1);
% opts3.minPkHt = 0.25;
% [hrBio, pkMaxBio] = ecgHR(ecgFilt,fs(2),opts3);
% 
% rmseHR = sqrt(mean((hrBio - hrNcs).^2));
% fprintf('\n RMSE HR: %3.2f\n',rmseHR);
% 
% %Peforming correction
% qualHeartFs = repmat(qualIdxHeart(:)',[qualEpochHeart*fs,1]); % Making sure qualIdxResp is initially a column, converted to a row vector
% qualHeartFs = qualHeartFs(:);
% qualHeartFs = qualHeartFs(1:size(ncsHeart,1));
% tOffHr = [5,5]; %  sec from start and end
% qualHeartFs(1:tOffHr(1)*fs) = 0;
% qualHeartFs(end-tOffHr(2)*fs+1:end) = 0;
% 
% hrBioNcsCorr = [hrBio, hrNcs].*qualHeartFs;
% hrBioNcsCorr = hrBioNcsCorr(qualHeartFs~=0,:);
% 
% errHrCorr = hrBioNcsCorr(:,1)-hrBioNcsCorr(:,2);
% rmseHRcorr = sqrt(mean(errHrCorr.^2));
% fprintf('RMSE HR after correction: %3.2f\n',rmseHRcorr);
% 
% figure
% clear ax1
% nFig = 3;
% ax1(1) = subplot(nFig,1,1);
% plot(tNcs,ncsHeart);hold on
% plot(tNcs(pkMaxNcs),ncsHeart(pkMaxNcs),'^');
% 
% ax1(2) = subplot(nFig,1,2);
% plot(tNcs,ecgFilt); hold on
% plot(tNcs(pkMaxBio),ecgFilt(pkMaxBio),'^');
% 
% ax1(3) = subplot(nFig,1,3);
% plot(tNcs,hrNcs); hold on
% plot(tNcs,hrBio);
% legend({'NCS','Ecg'});
% linkaxes(ax1,'x');

% --------------------------------------------------------------------
% Using RR interval
opts.tWinHR = 3;
opts.tWin = 0.5;
opts.minInterceptDist = 0.1;
opts.minEcgPkHt = 0.35;
opts.pkAmpRelRejThresh = 0.2; 
opts.tRRthresh = [0.4,1.2]; 

% Getting harmonic of heartbeat
opts2.filtType = 'LpHp';
opts2.f3db = 1.8; opts2.fpLP = 3.0; opts2.fstLP = 3.3;
ncsHeart2Harm = filterLpHp(ncs(:,idxHeart),fs(1),opts2);

opts.harmNum = 1;
[hrvNcs,hrvBio,~,~] = hrvFeatureEst(ncsHeart,ecgFilt,fs,opts);

opts.harmNum = 2;

[hrvNcs2ndHarm,~,~,~] = hrvFeatureEst(ncsHeart2Harm,[],fs,opts);

hrBio = resample(hrvBio.hr,hrvBio.tRR,fs(1));
hrNcs = resample(hrvNcs2ndHarm.hr,hrvNcs2ndHarm.tRR,fs(1));
% As last tRR above may be different, even after resampling length issue
% may exist
if length(hrBio)>length(hrNcs)
    hrBio = hrBio(1:length(hrNcs));
    tRR_resample = 0:1/fs(1):(length(hrNcs)-1)/fs(1);
else
    hrNcs = hrNcs(1:length(hrBio));
    tRR_resample = 0:1/fs(1):(length(hrBio)-1)/fs(1);
end
rmseHR = sqrt(mean((hrNcs-hrBio).^2));
fprintf('\n RMSE HR: %3.2f\n',rmseHR);

%Peforming correction
qualHeartFs = repmat(qualIdxHeart(:)',[qualEpochHeart*fs,1]); % Making sure qualIdxResp is initially a column, converted to a row vector
qualHeartFs = qualHeartFs(:);
qualHeartFs = qualHeartFs(1:size(hrBio,1));
tOffHr = [5,5]; %  sec from start and end
qualHeartFs(1:tOffHr(1)*fs) = 0;
qualHeartFs(end-tOffHr(2)*fs+1:end) = 0;

hrBioNcsCorr = [hrBio, hrNcs].*qualHeartFs;
hrBioNcsCorr = hrBioNcsCorr(qualHeartFs~=0,:);
errHrCorr = hrBioNcsCorr(:,1)-hrBioNcsCorr(:,2);
tRR_resample = tRR_resample(qualHeartFs~=0);
rmseHRcorr = sqrt(mean(errHrCorr.^2));
fprintf('RMSE HR after correction: %3.2f\n',rmseHRcorr);

figure
ax1(1) = subplot(3,1,1);
plot(hrvNcs.tRR,hrvNcs.hr)
hold on
plot(hrvNcs2ndHarm.tRR,hrvNcs2ndHarm.hr)
plot(hrvBio.tRR,hrvBio.hr)
legend({'NCS','NCS 2^{nd} Harm','ECG'})
ylabel('Heart Rate (BPM)')
xlabel('Time (s)')
ax1(2) = subplot(3,1,2);
plot(tQualIdxHeart,qualIdxHeart)

ax1(3) = subplot(3,1,3);
plot(tRR_resample,hrBioNcsCorr(:,1));hold on
plot(tRR_resample,hrBioNcsCorr(:,2));
linkaxes(ax1,'x');
% -------------------------------------------------------------------------
fprintf('Done.\n');

%%
pause(nSec);
close all

%% ------------------------------------------------------------------------
% savefig(fig, [dataPath,'Analysis\',ncsFile,'_linCalibVol.fig']);
% save([dataPath,'Analysis\',ncsFile,'_linCalibVol.mat'], 'ncs','bio','fs',...
%     'ncsCalib','bioCalib','meanDev','ncsRespAbd','ncsRespTh','tvAirflow',...
%     'volAirflow','beltCalibCoeff','vBeltCalib','ncsCalibCoeff','vNcsCalib',...
%     'vBelt','vNcs','tBio','tNcs','tBioCalib','tNcsCalib','opts1','opts2',...
%     'tStartEndOff','tManualOff','ncsCalibFile','tCalibStartEndOff',...
%     'signNcs','signNcsCalib','tManualOffCalib');

save(['C:\Research\NCS\HumanStudyData\','Analysis2_Sept2019\case',num2str(caseNum),'_',ncsFile(12:end),'_Vol_BR_HR.mat'], 'ncs','bio','fs',...
    'ncsCalib','bioCalib','ncsRespAbd','ncsRespTh','ncsHeart','tvAirflow',...
    'volAirflow','beltCalibCoeff','ncsCalibCoeff',...
    'tNcs','tNcsCalib','opts1','opts2',...
    'tOffNcsStEnd','tOffNcsStEndCalib','ncsCalibFile',...
    'signNcs','signNcsCalib',...
    'vBeltCalib','vNcsCalib','vBelt','vNcs',...
    'brNcs','ncsRespPk','brBio','bioRespPk','tQualIdxResp','qualIdxResp','qualRespConf',...
    'hrNcs','hrBio','tQualIdxHeart','qualIdxHeart','threshPred',...
    'X','modelHeartMotion','scorePred','kFoldLossPred','scorePredOrig','scorePredCalib',...
    'tDevCalib','tDevCalibThAbd','tDev','tDevThAbd',...
    'rRespThBio','rRespAbdBio','rRespCalibThBio','rRespCalibAbdBio','rVecHeart','corrFacHeart');

fprintf(' **************************** \nCase %d ONE Routine completed & saved..\n **************************** \n',caseNum);

end
