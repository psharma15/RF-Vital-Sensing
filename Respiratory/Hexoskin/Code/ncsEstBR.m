% This function estimates Breath Rate (BR) from NCS 
% April 18, 2018
% Pragya Sharma, ps847@cornell.edu

function [ncsBR,tBR] = ...
    ncsEstBR(ncsRespDS,inExAmp,inExPh,hxBR,ncsDownSampRate,hxSampRateBR)
%% Input and Output:
% ncsRespDS: [ncs amp, ncs ph]
% hxBR: Hexoskin Breath Rate estimation
% ncsDownSampRate: sample frequency of input amp and ph data
% hxSampRateBR: BR sample frequency for input Hx data
% ncsBR: BR calculated from NCS amp and ph [ampBR, phBR]
% tBR: ncsBR calculation points. Sampling time is defined by hxSampRateTV

%% ------------------------------------------------------------------------
% Checking correct NCS data input format
[dataRow, dataCol] = size(ncsRespDS);
if (dataRow < 1) || (dataCol ~=2)
    fprintf('In ncsEstBR(), input data is expected in [amp,ph] format.');
    return;
end

%% ------------------------------------------------------------------------
% Should ensure that there is an exhalation for every inhalation.
tData = (0:(1/ncsDownSampRate):((length(ncsRespDS(:,1))-1)/ncsDownSampRate))';
tBR = (0:(1/hxSampRateBR):((length(hxBR)-1)/hxSampRateBR))';

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

% Using inhalation point for BR calculation, as peak is more accurate.
tAmpBR = tData(ampInhaleIdx); 

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

tPhBR = tData(phInhaleIdx);

%% ------------------------------------------------------------------------
% Finding NCS BR
ncsAmpBR = zeros(length(tBR),1);
ncsPhBR = zeros(length(tBR),1);

tCycle = 0;
idxAmpBRold = 1; % Previous cycle inhalation/ exhalation.
% Iteration starts from 2: 2 inhalation points for 1 breath cycle.
for iter = 2:length(tAmpBR)
    if iter <= 7
        nCycle = iter - 1;
        tCycle = tCycle + tAmpBR(iter) - tAmpBR(iter-1);
    else
        nCycle = 7;
        tCycle = tAmpBR(iter) - tAmpBR(iter-7);
    end
    idxAmpBR = find(tBR >= tAmpBR(iter),1);
    ncsAmpBR(idxAmpBRold:idxAmpBR-1) = ncsAmpBR(idxAmpBRold);
    ncsAmpBR(idxAmpBR) = nCycle/tCycle;
    idxAmpBRold = idxAmpBR;
end

tCycle = 0;
idxPhBRold = 1;
for iter = 2:length(tPhBR)
    if iter <= 7
        nCycle = iter - 1;
        tCycle = tCycle + tPhBR(iter) - tPhBR(iter-1);
    else
        nCycle = 7;
        tCycle = tPhBR(iter) -tPhBR(iter-7);
    end
    idxPhBR = find(tBR >= tPhBR(iter),1);
    ncsPhBR(idxPhBRold:idxPhBR-1) = ncsPhBR(idxPhBRold);
    ncsPhBR(idxPhBR) = nCycle/tCycle;
    idxPhBRold = idxPhBR;
end

ncsBR = [60.*ncsAmpBR(:), 60.*ncsPhBR(:)];
end
