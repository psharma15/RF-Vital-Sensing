% This code reads NCS data, takes acquisition start time as input, find
% section of Hexoskin RAW data synched with that *(nearest 1-2 seconds)*.
% April 05, 2018
% Pragya Sharma, ps847@cornell.edu

function [ncsDataTrunc,ncsTimeTrunc,ncsSampRate,hxDataTrunc,hxTimeTrunc,hxSampRate] = ...
    ncsHxRawSync(dataPath,hxFolder,hxDataNum,ncsDataNum,...
    manualTimeOffset,dataDuration,ncsTstart)
    
%% Input to function
% dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
%     'Hexoskin\Data\2'];
% hxFolder = '\user_13412';
% hxDataNum = 2; % dataNum: 1 = resp_abd, 2 = resp_thrx, 3 = ecg
% ncsDataNum = 2; % data at different time instants
% manualTimeOffset = 20.5; % sec: This is by observation.
% dataDuration = 0;
% ncsTstart = 0*60; % Time is relative to NCS in seconds

%% Reading Hexoskin data
[hxData, hxSampRate, hxDateTime] = ...
    readHxData([dataPath, hxFolder],hxDataNum);

%% Reading NCS data
[ncsData, ncsSampRate, ncsStartDateTime] = ...
    readNcsData(dataPath,ncsDataNum);

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

hxTime =   0 : 1/hxSampRate : (length(hxData)-1)/hxSampRate ;
ncsTime =  0 : 1/ncsSampRate:(length(ncsData)-1)/ncsSampRate;

figure
nFig = 2;
ax(1) = subplot(nFig,1,1);
plot(hxTime,hxData); 
ylabel('Hx Data'); xlabel('Entire Time (sec)');
title(['Hx data collection started on: ',num2str(hxDateTime(1,1)),'/',...
    num2str(hxDateTime(1,2)),'/',num2str(hxDateTime(1,3)),' ',...
    num2str(hxDateTime(1,4)),':',num2str(hxDateTime(1,5)),':',...
    num2str(hxDateTime(1,6))])

ax(2) = subplot(nFig,1,2);
plot(ncsTime+hxOffsetSeconds,ncsData(:,1)); 
ylabel('NCS Amp Data'); xlabel('Time (sec)');
title(['Ncs data collection started on: ',num2str(ncsStartDateTime(1)),'/',...
    num2str(ncsStartDateTime(2)),'/',num2str(ncsStartDateTime(3)),' ',...
    num2str(ncsStartDateTime(4)),':',num2str(ncsStartDateTime(5)),':',...
    num2str(ncsStartDateTime(6))])

%%
if dataDuration == 0 
    dataDuration = ncsTime(end); % Data duration in seconds
end

if (dataDuration + ncsTstart) > ncsTime(end) 
    dataDuration = ncsTime(end) - ncsTstart;
    fprintf(['Data duration is modified to %f sec, as otherwise it ',...
        'exceeds total data time.\n'],dataDuration);
    if dataDuration <= 1/hxSampRate
        fprintf(['NCS data starting time is more than total duration',...
            'code stopping. \n']);
    end
end

% Add 1, as idx 1 == time 0 sec
ncsStartIdx = ncsTstart*ncsSampRate + 1; 
ncsEndIdx = (ncsTstart + dataDuration)*ncsSampRate + 1;

% Truncating time, but losing original time. 
ncsTimeTrunc = ncsTime(ncsStartIdx:ncsEndIdx) - ncsTime(ncsStartIdx);
ncsDataTrunc = ncsData(ncsStartIdx:ncsEndIdx,:);

hxStartIdx = uint64((hxOffsetSeconds + ncsTstart)*hxSampRate)+1;
hxEndIdx = uint64((hxOffsetSeconds + ncsTstart + dataDuration)*hxSampRate)+1;

hxTimeTrunc = hxTime(hxStartIdx:hxEndIdx) - hxTime(hxStartIdx); % Updating for the short duration
hxDataTrunc = hxData(hxStartIdx:hxEndIdx);

figure
nFig = 2;
ax(1) = subplot(nFig,1,1);
plot(hxTimeTrunc,hxDataTrunc); 
ylabel('Hx Data'); xlabel('Time (sec)');

ax(2) = subplot(nFig,1,2);
yyaxis left
plot(ncsTimeTrunc,ncsDataTrunc(:,1)); 
ylabel('NCS Amp');

yyaxis right
plot(ncsTimeTrunc,ncsDataTrunc(:,2));
ylabel('NCS Ph');

xlabel('Time (sec)');
 
linkaxes(ax,'x')

