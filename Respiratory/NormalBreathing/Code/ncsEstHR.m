%% Find Heart rate from NCS amp and Ph
% 26 April 2018
% Pragya Sharma, ps847@cornell.edu

function ncsHR = ...
    ncsEstHR(ncsHeart,heartAmpMinMax,heartPhMinMax,ncsHeartSampRate,tHR)
% Derived entirely from BR estimate

%% ------------------------------------------------------------------------
% Checking correct NCS data input format
[dataRow, dataCol] = size(ncsHeart);
if (dataRow < 1) || (dataCol ~=2)
    fprintf('In ncsEstHR(), input data is expected in [amp,ph] format.');
    return;
end
    
%% ------------------------------------------------------------------------
tData = (0:(1/ncsHeartSampRate):((length(ncsHeart(:,1))-1)/ncsHeartSampRate))';


ampMaxIdx = heartAmpMinMax((heartAmpMinMax(:,2) == 1),1);
ampMinIdx = heartAmpMinMax((heartAmpMinMax(:,2) == 0),1);
 
phMaxIdx = heartPhMinMax((heartPhMinMax(:,2) == 1),1);
phMinIdx = heartPhMinMax((heartPhMinMax(:,2) == 0),1);    

% Using minimas for a period, that is more sharp, hence accurate
tAmpHR = tData(ampMinIdx);
tPhHR = tData(phMinIdx);

%% ------------------------------------------------------------------------
% Finding NCS HR
ncsAmpHR = zeros(length(tHR),1);
ncsPhHR = zeros(length(tHR),1);

tCycle = 0;
idxAmpHRold = 1; % Previous cycle inhalation/ exhalation.
% Iteration starts from 2: 2 inhalation points for 1 breath cycle.
nCycleAvg = 16; % Number of cycles except at the beginning.
for iter = 2:length(tAmpHR)
    if iter <= nCycleAvg
        nCycle = iter - 1;
        tCycle = tCycle + tAmpHR(iter) - tAmpHR(iter-1);
    else
        nCycle = nCycleAvg;
        tCycle = tAmpHR(iter) - tAmpHR(iter-nCycleAvg);
    end
    idxAmpHR = find(tHR >= tAmpHR(iter),1);
    ncsAmpHR(idxAmpHRold:idxAmpHR-1) = ncsAmpHR(idxAmpHRold);
    ncsAmpHR(idxAmpHR) = nCycle/tCycle;
    idxAmpHRold = idxAmpHR;
end

tCycle = 0;
idxPhHRold = 1;
for iter = 2:length(tPhHR)
    if iter <= nCycleAvg
        nCycle = iter - 1;
        tCycle = tCycle + tPhHR(iter) - tPhHR(iter-1);
    else
        nCycle = nCycleAvg;
        tCycle = tPhHR(iter) -tPhHR(iter-nCycleAvg);
    end
    idxPhHR = find(tHR >= tPhHR(iter),1);
    ncsPhHR(idxPhHRold:idxPhHR-1) = ncsPhHR(idxPhHRold);
    ncsPhHR(idxPhHR) = nCycle/tCycle;
    idxPhHRold = idxPhHR;
end

ncsHR = [60.*ncsAmpHR(:), 60.*ncsPhHR(:)];
end
