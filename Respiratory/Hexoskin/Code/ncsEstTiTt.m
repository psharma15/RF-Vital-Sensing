% This function estimates fractional inspiratory time (Ti/Tt, ratio of
% inspiratory to total breath time). 
% Inhalation: Defined as end of inhalation
% Exhalation: Defined as end of exhalation

function [ampTiTt,phTiTt] = ncsEstTiTt(ncsRespDS,inExAmp,inExPh,ncsDownSampRate)
% Input:
%     ncsRespDS: Downsampled NCS respiratory waveform, arranged as
%                [amp(:), ph(:)]. 
%     inExAmp: Amplitude inhalation and exhalation indices in first column,
%              corresponding indicator in second column.
%     inExPh: Phase inhalation and exhalation indices in first column,
%             corresponding indicator in second column.
%     ncsDownSampRate: Sampling rate of NCS data in ncsRespDS
% Output:
%     ampTiTt: Ti/Tt for each breath in second column for amplitude 
%              waveform, corresponding inhalation time (end of inspiration)
%              in first column. 
%     phTiTt: Ti/Tt for each breath in second column for phase waveform, 
%             corresponding inhalation time (end of inspiration) in first 
%             column.

%% ------------------------------------------------------------------------
% Finding breath start point (Exhalation point: Inhalation beginning), as
% finding inhalation point requires that information

% Defining tData as a column vector
tData = (0:(1/ncsDownSampRate):((length(ncsRespDS(:,1))-1)/ncsDownSampRate))';

if inExAmp(1,2) == 1
    % If data starts with inhalation (exhalation beginning), remove that 
    % and start measurement from next cycle. Truncate if last
    % point is inhalation (inhalation end). Keep last point as exhalation
    % to have complete number of cycles
    if inExAmp(end,2) == 1
        inExAmp = inExAmp(2:end-1,:);
    else
        inExAmp = inExAmp(2:end,:);
    end
else
    if inExAmp(end,2) == 1
        % If starting from exhalation (inhalation beginning), keep that,
        % but if ending at inhalatin (exhalation beginning), remove that.
        inExAmp = inExAmp(1:end-1,:);
    end
end
ampInhaleIdx = inExAmp((inExAmp(:,2) == 1),1);
ampExhaleIdx = inExAmp((inExAmp(:,2) == 0),1);

tInhaleAmp = tData(ampInhaleIdx);
tExhaleAmp = tData(ampExhaleIdx);

if inExPh(1,2) == 1
    if inExPh(end,2) == 1
        inExPh = inExPh(2:end-1,:);
    else
        inExPh = inExPh(2:end,:);
    end
else
    if inExPh(end,2) == 1
        inExPh = inExPh(1:end-1,:);
    end
end
phInhaleIdx = inExPh((inExPh(:,2) == 1),1);
phExhaleIdx = inExPh((inExPh(:,2) == 0),1);    

tInhalePh = tData(phInhaleIdx);
tExhalePh = tData(phExhaleIdx);

%% ------------------------------------------------------------------------
% Finding inspiration time (Ti) and total time (Tt) for each breath and
% Ti/Tt for both amplitude and phase.
ampTi = tInhaleAmp - tExhaleAmp; 
ampTt = diff(tExhaleAmp);
ampTiTt = [tInhaleAmp, ampTi./ampTt];

phTi = tInhalePh - tExhalePh;
phTt = diff(tExhalePh);
phTiTt = [tInhalePh, phTi./phTt];


