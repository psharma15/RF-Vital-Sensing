% Updated on 07 May 2019 
% -------------------------------------------------------------------------
% This function finds minima and maxima pairs, corresponding to end of
% inspiration and end of expiration in respiration waveforms.
% The peak detection is performed by AMPD (automatic multiscale based peak
% detection). The source is: "An efficient algorithm for automatic peak
% detection in noisy periodic and quasi-periodic signals", F Scholkmann.
% The data is then processed to find pairs of maxima and minima based on
% certain rules. Refer to comments for details.
% -------------------------------------------------------------------------
% Pragya Sharma, ps847@cornell.edu
% April 15th, 2018.

function pkIdxInd = findMaxMin(data,fs,opts)
% Input:
%   data: [Can be a vector or a matrix. For a matrix, each column will be
%   considered separately.
%   fs: Sampling frequency of data
%   freqSeg: A segment's frequency
% Output:
%   pkIdxInd: Structure with pkIdxInd(i).idx, pkIdxInd(i).ind giving
%   indices and max min indicator for i-th column of data.
% -------------------------------------------------------------------------

[nRow,nCol] = size(data);
t = (0:(nRow-1))./fs;

pkIdxInd = struct([]); %Initialize output struct
f1 = figure;    title('Uncorrected Peaks');
f2 = figure;    title('Corrected Peaks');

for colNum = 1:nCol

    % Finding peaks: both maximum and minimum
    % [locMax,locMin] = peakDetAMPD(data(:,colNum),fs,opts);
    [locMin,locMax] = peakDet3(data(:,colNum),fs,opts);
    
    % Somehow repeated inhale exhale points - most likely due to uncorrected
    % overlap segments. Can be corrected if nearby. But some issue if exact
    % point repeated. So checking that
    for j = 2:length(locMax)
        if locMax(j) == locMax(j-1)
            locMax(j) = 0;
        end
    end
    locMax = locMax(locMax ~= 0);

    for k = 2:length(locMin)
        if locMin(k) == locMin(k-1)
            locMin(k) = 0;
        end
    end
    locMin = locMin(locMin ~= 0);

    figure(f1);
    ax(colNum) = subplot(nCol,1,colNum); %#ok<*AGROW>
    plot(t,data(:,colNum)); 
    xlabel('Time (sec)');     ylabel(['column: ',colNum]);
    hold on
    plot(t(locMax),data(locMax,colNum),'^');
    plot(t(locMin),data(locMin,colNum),'v');
    legend('data','max','min');
    ax(colNum).XGrid = 'on';

    % Forming vector of peak indices: combining sorted maxima & minima pts.  
    % Corresponding indicator: 1: Max, 0: Min
    pkInd = [ones(length(locMax),1); zeros(length(locMin),1)];

    [pkIdx, sortIdx] = sort([locMax(:); locMin(:)]);
    pkInd = pkInd(sortIdx);

    % Keeping maxima and minima indices in pairs, dropping if consecutive
    % max or min points are found
    pkIdxInd(colNum).idx = zeros(length(pkIdx),1);
    pkIdxInd(colNum).ind = zeros(length(pkIdx),1); % pkIdx and pkInd of same length
    pkIdxInd(colNum).idx(1) = pkIdx(1);
    pkIdxInd(colNum).ind(1) = pkInd(1);
    
    counter = 2;
    tempMaxMin = [0,0,0]; % saves [inExAmpIdx, inExAmpIndicator, counter]
    dataMaxMin = 1e10;

    for iter = 2:length(pkIdx)
        if pkInd(iter) ~= pkInd(iter-1)
            if tempMaxMin(1) ~= 0
                % Updating max/ min if there were any multiple consecutive
                % max/ min, with a different peak than the one saved for that
                % counter.
                pkIdxInd(colNum).ind(tempMaxMin(3)) = tempMaxMin(2);
                pkIdxInd(colNum).idx(tempMaxMin(3)) = tempMaxMin(1);
                tempMaxMin = [0,0,0];
            end
            pkIdxInd(colNum).idx(counter) = pkIdx(iter);
            pkIdxInd(colNum).ind(counter) = pkInd(iter);
            counter = counter + 1;
            dataMaxMin = data(pkIdx(iter),colNum); 
        end
        if (pkInd(iter) == pkInd(iter-1))
            if dataMaxMin == 1e10
                % This is the max/ min data pt among consecutive max/ min
                dataMaxMin = data(pkIdx(iter-1),colNum); 
            end
            if ((pkInd(iter) == 0) && (data(pkIdx(iter),colNum) > dataMaxMin)) ...
                    || ((pkInd(iter) == 1) && (-data(pkIdx(iter),colNum) > -dataMaxMin))
                dataMaxMin = data(pkIdx(iter),colNum);
                tempMaxMin = [pkIdx(iter), pkInd(iter), ...
                    counter-1];
            end
        end
    end
    
    pkIdxInd(colNum).idx = pkIdxInd(colNum).idx(1:(counter-1));
    pkIdxInd(colNum).ind = pkIdxInd(colNum).ind(1:(counter-1));

    maxIdx = pkIdxInd(colNum).idx(pkIdxInd(colNum).ind == 1);
    minIdx = pkIdxInd(colNum).idx(pkIdxInd(colNum).ind == 0);

    % Plotting corrected locations
    figure(f2);
    ax1(colNum) = subplot(nCol,1,colNum);
    plot(t,data(:,colNum)); 
    xlabel('Time (sec)'); ylabel(['column: ',num2str(colNum)]);
    hold on
    plot(t(maxIdx),data(maxIdx,colNum),'^');
    plot(t(minIdx),data(minIdx,colNum),'v');
    ax1(colNum).XGrid = 'on';
    
end
linkaxes(ax,'x');
linkaxes(ax1,'x');
