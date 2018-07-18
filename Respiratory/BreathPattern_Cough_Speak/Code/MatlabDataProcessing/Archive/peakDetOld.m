% This algorithm finds maxima and minima of a timeseries using AMPD
% algorithm. 
% "An efficient algorithm for automatic peak detection in noisy periodic
% and quasi-periodic signals", F Scholkmann.
% *********************************************************************** %
% Leaves out peak detection for beginning and end data. Use updated peak
% detection algorithm.
% *********************************************************************** %

function [maximaIdx, minimaIdx] = peakDetOld(signal, freqRange, fs)
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
nCyclesPerSegment = 4; 
n = (nCyclesPerSegment*fs/freqRange(1));
lenSegment = 2*n + 1;
lenOverlap = n + 1;

% ------------------------------------------------------------------------
% Starting the loop
t = (0:1/fs:((length(signal)- 1)/fs))';
iterStartStop = [1, lenSegment]; % Start and stop indices of a segment
maxSegment = 2*ceil(length(signal)/lenSegment)-2; % Check this
maximaIdx = [];
minimaIdx = [];

for segmentNum = 1:maxSegment
    
    signalSegment = signal(iterStartStop(1):iterStartStop(2));
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
    
    maximaIdx = [maximaIdx, (iterStartStop(1) + maximaPt - 1)];
    minimaIdx = [minimaIdx, (iterStartStop(1) + minimaPt - 1)];
    
    % ---------------------------------------------------------------------
    % For next iteration
    linearFitOld = linearFit;
    iterStartStop = iterStartStop + n;
    if iterStartStop(2) > length(signal)
        iterStartStop(2) = length(signal);
    end    

end
