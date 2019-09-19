% This algorithm finds maxima and minima of a timeseries using AMPD
% algorithm. 
% "An efficient algorithm for automatic peak detection in noisy periodic
% and quasi-periodic signals", F Scholkmann.

function [minimaIdx, maximaIdx] = peakDetAMPD(signal, freqRange, fs)
% Input: 
%       signal: Input signal, 1-D array
%       freqRange: Expected frequency range of input signal [min, max] Hz
%       fs: Sampling frequency of signal in Hz
% Output:
%       peakLoc: Index of the peak
% ------------------------------------------------------------------------
% Dividing the signal in segments 
% Length: 2 * n + 1
% Overlap: n + 1 
% n is given by minimum expected frequency. Need few quasi-periodic cycles
% in a segment. Taking number of quasi periodic cycles at MINIMUM frequency
% to be 4.
% This is an important parameter that affects detection of multiple wrong
% peaks when breath rate is very low. 

signal = signal(:); % Make sure that signal is a column vector
nCyclesPerSegment = 1; 

% *********************************************************************** %
% Change n based on instantaneous frequency of the signal, as there is
% issue if there is a sudden frequency shift. Acc to AMPD paper, it can
% give sufficient performance, if fmax < 4*fmin, which is not the case from
% a slow to fast-shallow breath shift. The freq range is fmax < 10*fmin.
% *********************************************************************** %
n = (nCyclesPerSegment*fs/freqRange(1));
lenSegment = 2*n + 1;
lenOverlap = n + 1;

% ------------------------------------------------------------------------
% Padding the signal at the startign and end so that peak detection is
% performed for the entire signal range. Because in each segment, AMPD
% method only finds peaks for signal in the range i = k+2:lenSignalDtr-k+1 
% (see below).
padLength = lenOverlap-1;
signalPadded = [signal(1).*ones(padLength,1);...
                signal;...
                signal(end).*ones(padLength,1)];
lenSignalPadded = length(signalPadded);
% ------------------------------------------------------------------------
% Starting the loop
t = (0:1/fs:((lenSignalPadded- 1)/fs))';
iterStartStop = [1, lenSegment]; % Start and stop indices of a segment
maxSegment = 2*ceil(length(signalPadded)/lenSegment)-2; % Check this
maximaIdx = [];
minimaIdx = [];

for segmentNum = 1:maxSegment
    
    signalSegment = signalPadded(iterStartStop(1):iterStartStop(2));
    tSegment = t(iterStartStop(1):iterStartStop(2));
    
    % ---------------------------------------------------------------------
    % Starting with detrending, using linear fit method. Using weighted
    % detrending for the overlap region.
    [linearFitCoeff,linearFitErr] = polyfit(tSegment,signalSegment,1);
    [linearFit,~] = polyval(linearFitCoeff,tSegment,linearFitErr);
    if segmentNum > 1 
        for sampleNum = 1 : lenOverlap
            wts = [1-((sampleNum-1)/n); (sampleNum-1)/n];
            linearFit(sampleNum) = wts(1)*linearFitOld(sampleNum+n) + ...
                wts(2)*linearFit(sampleNum);
        end
    end
    signalDtr = signalSegment - linearFit;
    lenSignalDtr = length(signalDtr);
    
    % ---------------------------------------------------------------------
    % Finding Local Maxima and Minima Scalogram (LMS): modifying the paper
    
    % ******************************************************************* %
    % REMEMBER: No implementation yet for overlap region. They might have
    % multiple resulting peaks because of that (Slightly shifted in time).
    % Currently taking care of this in the post-processing.
    % ******************************************************************* %
    
    numScaleMin = 1;
%     numScaleMax = ceil(fs/freqRange(1));
    numScaleMax = ceil(lenSignalDtr/2) - 1; % Number of scales
    numScale = numScaleMax - numScaleMin + 1;
    rng('default');
    lMAXs = zeros(numScale,lenSignalDtr)+rand(numScale,lenSignalDtr);
    lMINs = zeros(numScale,lenSignalDtr)+rand(numScale,lenSignalDtr);
    
    % Iterations to find local maxima and minima at different scales
    % 1: Condition satisfied (minima/ maxima)
    % 0+rand: Condition not satified 
    for k = numScaleMin:numScaleMax
        for i = k+2:lenSignalDtr-k+1
            if (signalDtr(i-1) > signalDtr(i-k-1)) && ...
                    (signalDtr(i-1) > signalDtr(i+k-1))
                lMAXs(k-numScaleMin+1,i) = 1;
            end
            if (-signalDtr(i-1) > -signalDtr(i-k-1)) && ...
                    (-signalDtr(i-1) > -signalDtr(i+k-1))
                lMINs(k-numScaleMin+1,i) = 1;
            end
        end
    end                
    
    % ---------------------------------------------------------------------
    % Finding minima and maxima based on lMINs and lMAXs respectively
    % Gamma contains information about the sclae-dependent distribution of
    % ones (maximas/ minimas). This is used for scaling of the LMS
    % matrices. Lambda gives scale for which maximum number of peaks occur,
    % and any window size greater than that is removed.
    gammaMAX = sum(lMAXs,2);
    [~,lambdaMAX] = max(gammaMAX);
    lMAXs = lMAXs(1:lambdaMAX,:);
    gammaMIN = sum(lMINs,2);
    [~,lambdaMIN] = max(gammaMIN);
    lMINs = lMINs(1:lambdaMIN,:);
    
    % Finding standard deviation of each column, which will be 0, if that
    % point is always a maxima or a minima for all the scales.
    stdDevMAX = std(lMAXs);
    maximaPt = find(stdDevMAX == 0);
    stdDevMIN = std(lMINs);
    minimaPt = find(stdDevMIN == 0);
    
    % Correcting index to represent index in unpadded signal vector. 
    maximaPt = iterStartStop(1) + maximaPt - 1 - padLength; 
    minimaPt = iterStartStop(1) + minimaPt - 1 - padLength;
    maximaPt = maximaPt(:);
    minimaPt = minimaPt(:);
    
    maximaPtUpdated = [];
    minimaPtUpdated = [];
    for i = 1:length(maximaPt)
        if (maximaPt(i) > 0) && (maximaPt(i) < length(signal))
            maximaPtUpdated = [maximaPtUpdated; maximaPt(i)];
        end
    end
    for i = 1:length(minimaPt)
        if (minimaPt(i) > 0) && (minimaPt(i) < length(signal))
            minimaPtUpdated = [minimaPtUpdated; minimaPt(i)];
        end
    end       

    maximaIdx = [maximaIdx; maximaPtUpdated];
    minimaIdx = [minimaIdx; minimaPtUpdated];
    
    % ---------------------------------------------------------------------
    % For next iteration
    linearFitOld = linearFit;
    iterStartStop = iterStartStop + n;
    if iterStartStop(2) > length(signalPadded)
        iterStartStop(2) = length(signalPadded);
    end    

end
