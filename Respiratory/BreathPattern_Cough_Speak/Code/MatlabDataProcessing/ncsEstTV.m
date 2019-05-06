% This function calibrates tidal volume from NCS against RIP.
% April 12, 2018
% Pragya Sharma, ps847@cornell.edu

function [tvCoeffAmpPhSum,tvCoeffAmp,tvCoeffPh,ncsTV] = ...
    ncsEstTV(ncsResp,inExAmp,inExPh,hxTV,ncsSampRate,hxSampRateTV,fitFunc,tOffsetTV,calibTime)
%% Input and Output:
% ncsResp: [ncs amp, ncs ph]
% hxTv: Hexoskin Tidal volume estimation
% ncsDownSampRate: sample frequency of input amp and ph data
% hxSampRateTV: Tidal volume sample frequency for input hxTv 
% calibTime: Time in seconds, used for calibration: [tStart, tEnd]
% tvCoeffAmpPhSum: Scaling coeff for getting tv = alpha (AMP) + beta (PH),
%                  arranged as [alpha, beta]
% tvCoeffAmp (or Ph): Scaling coefficient for tidal volume estimation from 
%          amplitude (or phase).
% ncsTV: TV calculated, without calibration, from amp and ph arranged as 
%        [ampTV,  phTV]. 
% tTV: ncsTV calculation points. Sampling time is defined by hxSampRateTV

%% ------------------------------------------------------------------------
% Checking correct NCS data input format
[dataRow, dataCol] = size(ncsResp);
if (dataRow < 1) || (dataCol ~=2)
    fprintf('In ncsEstTV(), input data is expected in [amp,ph] format.');
    return;
end 

%% ------------------------------------------------------------------------
% Should ensure that there is an exhalation for every inhalation.
tData = (0:(1/ncsSampRate):((length(ncsResp(:,1))-1)/ncsSampRate))';
tTV = (0:(1/hxSampRateTV):((length(hxTV)-1)/hxSampRateTV))'+tOffsetTV;

if inExAmp(1,2) == inExAmp(end,2)
    % If data starts and stops with the same event, truncating it such that
    % first event is inhalation.
    if inExAmp(1,2) == 0
        inExAmp = inExAmp(2:end,:);
    else
        inExAmp = inExAmp(1:end-1,:);
    end
    ampInhaleIdx = inExAmp(1:2:length(inExAmp(:,1))-1,1);
    ampExhaleIdx = inExAmp(2:2:length(inExAmp(:,1)),1);
else
    ampInhaleIdx = inExAmp((inExAmp(:,2) == 1),1);
    ampExhaleIdx = inExAmp((inExAmp(:,2) == 0),1);
end
tInhaleAmp = tData(ampInhaleIdx);
tExhaleAmp = tData(ampExhaleIdx);

if inExPh(1,2) == inExPh(end,2)
    % If data starts and stops with the same event, truncating it such that
    % first event is inhalation.
    if inExPh(1,2) == 0
        inExPh = inExPh(2:end,:);
    else
        inExPh = inExPh(1:end-1,:);
    end
    phInhaleIdx = inExPh(1:2:length(inExPh(:,1))-1,1);
    phExhaleIdx = inExPh(2:2:length(inExPh(:,1)),1);
else
    phInhaleIdx = inExPh((inExPh(:,2) == 1),1);
    phExhaleIdx = inExPh((inExPh(:,2) == 0),1);    
end
tInhalePh = tData(phInhaleIdx);
tExhalePh = tData(phExhaleIdx);

ampTV = ncsResp(ampInhaleIdx,1)-ncsResp(ampExhaleIdx,1); %#ok<*FNDSB>
% tAmp = (tInhaleAmp + tExhaleAmp)./2;
tAmp = tInhaleAmp;

phTV = ncsResp(phInhaleIdx,2)-ncsResp(phExhaleIdx,2);
% tPh = (tInhalePh + tExhalePh)./2;
tPh = tInhalePh;

%% ------------------------------------------------------------------------
% Correcting wrong TV if wrong inspiration and respiration are detected
% This is the case when small wrong peaks are detected as inspiration and
% expiration. Or there is a motion artefact and huge jump.

%% ------------------------------------------------------------------------
% Extrapolation: Need to find tidal volume (TV) at same sample frequency as 
% Hexoskin TV for comparison
ampTVextrap = interp1(tAmp,ampTV,tTV,'previous','extrap');
phTVextrap = interp1(tPh,phTV,tTV,'previous','extrap');
ncsTV = [ampTVextrap(:), phTVextrap(:)]; % This is output

%% ------------------------------------------------------------------------
% Plot
figure
ax(1) = subplot(3,1,1);
plot(tTV,hxTV); 
ylabel('Hexoskin TV (mL)'); xlabel('Time (sec)'); grid on

ax(2) = subplot(3,1,2);
plot(tTV,ampTVextrap);
ylabel('NCS Amp TV (au)'); xlabel('Time (sec)'); grid on

ax(3) = subplot(3,1,3);
plot(tTV,phTVextrap)
ylabel('NCS Ph TV (au)'); xlabel('Time (sec)'); grid on

linkaxes(ax,'x')

%% ------------------------------------------------------------------------
% Finding calibration coefficient

% *********************************************************************** %
% Ignore first and last very few seconds of ampTV, phTV - although peak
% detection code is corrected to perform peak detection here, best to leave
% the first and last peak. This can be taken care of in calibration time.
% *********************************************************************** %

if (length(calibTime) ~= 2) || (calibTime(2) == 0)
    calibTime = [tTV(1), tTV(2)];
end
calibStartIdx = find((abs(tTV - calibTime(1)) < (1/hxSampRateTV)),1);
calibStopIdx = find((abs(tTV - calibTime(2)) < (1/hxSampRateTV)),1);

fprintf('Calibration time is %f - %f sec.\n',tTV(calibStartIdx),...
        tTV(calibStopIdx));

ampTVcalib = ampTVextrap(calibStartIdx:calibStopIdx);
phTVcalib = phTVextrap(calibStartIdx:calibStopIdx);
hxTVcalib = hxTV(calibStartIdx:calibStopIdx);

% Exclude hxTV == 0 and NaN resulting from extrapolation in amp and ph TV
% estimation.
idxFit = find((hxTVcalib ~= 0) & (~isnan(ampTVcalib)) & (~isnan(phTVcalib))); 

ampTVcalibTrunc = ampTVcalib(idxFit); % Truncated amplitude and phase data to fit
phTVcalibTrunc = phTVcalib(idxFit);
hxTVcalibTrunc = hxTVcalib(idxFit);

% Find Coefficient A, and B in: HxTV = A*ncsAmpTV + B*ncsPhTV
% Using MATLAB's curve/ surface fitting tool.
fitAmpPh = @(A,B,x,y) A.*x + B.*y;
startPt = [0.5e5, 40]; % Approximate A and B points to start 

modelCalibrationAmpPhSum = fit([ampTVcalibTrunc(:) phTVcalibTrunc(:)],...
                               hxTVcalibTrunc,fitAmpPh, ...
                               'StartPoint',startPt);
                  
tvCoeffAmpPhSum = coeffvalues(modelCalibrationAmpPhSum);

% Fitting Amp
switch(fitFunc)
    case 'linear'
        fitAmp = @(A,x) A.*x;
        modelCalibrationAmp = fit(ampTVcalibTrunc(:),hxTVcalibTrunc(:),fitAmp,...
                          'StartPoint',startPt(1));
    case 'quad'
        fitAmp = @(A,B,C,x) A.*(x.^2) + B.*x + C;
        modelCalibrationAmp = fit(ampTVcalibTrunc(:),hxTVcalibTrunc(:),fitAmp,...
                          'StartPoint',[0.5e2,0.5e5,0]);

    otherwise
        fprintf('Enter correct fitting function. \n'); 
end
tvCoeffAmp = coeffvalues(modelCalibrationAmp);

% Fitting Phase
fitPh = @(B,x) B.*x;
modelCalibrationPh = fit(phTVcalibTrunc(:),hxTVcalibTrunc(:),fitPh,...
                         'StartPoint',startPt(2));
tvCoeffPh = coeffvalues(modelCalibrationPh);

end