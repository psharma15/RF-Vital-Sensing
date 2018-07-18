% This code reads NCS data, takes acquisition start time as input, find
% section of Hexoskin RAW data synched with that *(nearest 1-2 seconds)*.
% April 24, 2018
% Pragya Sharma, ps847@cornell.edu

function [ncsDataTrunc,ncsTimeTrunc,hxDataTrunc,hxDateTimeTrunc,hxSampRate] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNum,ncsFileName,ncsSampRate,...
    manualTimeOffset,dataDuration,ncsTstart)
    
%% Input to function
% dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
%     'Hexoskin\Data\2'];
% hxFolder = '\user_13412';
% hxDataNum = 2; % dataNum: 1 = resp_abd, 2 = resp_thrx, 3 = ecg
% ncsFileName = 'freq2G_1v2a'; % data at different time instants
% manualTimeOffset = 20.5; % sec: This is by observation.
% dataDuration = 0;
% ncsTstart = 0*60; % Time is relative to NCS in seconds

%% Reading Hexoskin data
[hxData, hxSampRate, hxDateTime] = ...
    readHxData([dataPath, hxFolder],hxDataNum);

%% Reading NCS data
[ncsData, ncsStartDateTime] = ...
    readNcsData(dataPath,ncsFileName);

%% Find Hx data start offset, assuming Hx data started before NCS data
% Compute the time elapsed to 0.01-second accuracy.
% IMP: Remember issue with the way Hx data is saved. It saves previous data
% records in the same file as new data sometimes.
hxSameDateDataStart = find(hxDateTime(:,3) == ncsStartDateTime(3),1);
hxData = hxData(hxSameDateDataStart:end);
hxDateTime = hxDateTime(hxSameDateDataStart:end,:);
hxOffsetSeconds = -1*etime(hxDateTime(1,:),ncsStartDateTime);

if hxOffsetSeconds < 0
    % If NCS started before, error and exit
    fprintf('NCS data collection started before Hx. ERROR.\n');
    return
end

hxOffsetSeconds = hxOffsetSeconds-manualTimeOffset; 
ncsTime =  0 : 1/ncsSampRate:(length(ncsData)-1)/ncsSampRate;

%%
if dataDuration == 0 
    dataDuration = ncsTime(end); % Data duration in seconds
end

if (dataDuration + ncsTstart) > ncsTime(end) 
    dataDuration = ncsTime(end) - ncsTstart;
    fprintf(['Data duration is modified to %f sec, as otherwise it ',...
        'exceeds total data time.\n'],dataDuration);
    if dataDuration <= 1/hxSampRate
        fprintf('Data duration is less than allowed. \n');
    end
end

% Add 1, as idx 1 == time 0 sec
ncsStartIdx = ncsTstart*ncsSampRate + 1; 
ncsEndIdx = (ncsTstart + dataDuration)*ncsSampRate + 1;

% Truncating time, but losing original time. 
ncsTimeTrunc = ncsTime(ncsStartIdx:ncsEndIdx) - ncsTime(ncsStartIdx);
ncsDataTrunc = ncsData(ncsStartIdx:ncsEndIdx,:);

hxTime = etime(hxDateTime,hxDateTime(1,:)); % Make starting time 0.
hxStartIdx = find(hxTime > (hxOffsetSeconds + ncsTstart),1,'first');
% hxStartIdx = uint64((hxOffsetSeconds + ncsTstart)*hxSampRate)+1;
hxEndIdx = find(hxTime > (hxOffsetSeconds + ncsTstart + dataDuration),1,'first');
% hxEndIdx = uint64((hxOffsetSeconds + ncsTstart + dataDuration)*hxSampRate)+1;

if (hxEndIdx > length(hxData))
    fprintf('Hexoskin data ended befor NCS, try with less data duration.\n');
    return
end

hxDataTrunc = hxData(hxStartIdx:hxEndIdx);

% This is the absolute time in date-time vector format, not starting from 0.
hxDateTimeTrunc = hxDateTime(hxStartIdx:hxEndIdx,:);


