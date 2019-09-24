% This program calculates Heart Rate (HR) and Heart Rate Variation (HRV)
% using NCS and benchmarks against 2 lead ECG.
% 27 Feb 2018
% ps847@cornell.edu

%% Reading data
dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Data\Mar03';
fileName = 'freq2G1';
tSynchStabilize = 0;
FsOld = 512;
FsNew = 200;
tDelay = 90e-3; % ECG and NCS delay

[ncsAmpData,ncsPhData,ecgData,t] = readEcgNcs(dataPath,fileName,tSynchStabilize,FsOld,FsNew,tDelay);
nSample = length(ncsAmpData);
if length(ecgData) ~= nSample
    fprintf('Ecg and Ncs have different length. \nStopping...\n');
    return
end

%% Perform Post Processing of NCS to view heartbeat
[ncsAmpFilt,ncsPhFilt,~,~,~] = postProcess(0.8,5,6,[ncsAmpData,ncsPhData],FsNew,[1,1]);
[ecgFilt,~,~,~,~] = postProcess(0.9,80,81,[ecgData, zeros(nSample,1)],FsNew,[1,1]);

% Plotting scaled versions of NCS and ECG
ncsAmp = 0.025e4*ncsAmpFilt;
ncsPh = ncsPhFilt;
ecg = 0.25e-3*ecgFilt;

% Removing first 6 seconds of data as filter output takes time to stabilize
tFiltStabilize = 6;
ncsAmp = ncsAmp(tFiltStabilize*FsNew+1:end);
ncsPh = ncsPh(tFiltStabilize*FsNew+1:end);
ecg = ecg(tFiltStabilize*FsNew+1:end);
t = t(1:end-(tFiltStabilize*FsNew));
nSample = nSample - 6*FsNew;

figure
grid on
yyaxis left
plot(t,ncsAmp); 
ylabel('NCS')
yyaxis right
plot(t,ecg); 
ylabel('ECG');

%% Beat detection in NCS

% Taking NCS as NCS Amp or Ph: Better one would be a weighted combination  
% of AMP and PH that gives the best heartbeat waveform retrieval
ncs = ncsPh; 
minPeakInterval = round(FsNew*60/130); % Apriori Condition: max heart rate 200 beats/min

minPeakHeightNcs = 0.3;
% Wavelet decomposition and reconstruction: D7 gives the best correlation
% wName = 'db10';
% [recNcsd7,~,~] = wavedecrec1(ncs,t,wName,7,1); 
% Not using recNcsd7 for HRV
[~,locsNcs] = findpeaks(ncs,'MinPeakDistance',minPeakInterval, ...
    'MinPeakHeight',minPeakHeightNcs);
nBeatNcs = length(locsNcs)-1;
rrIntervalNcs = t(locsNcs(2:end))-t(locsNcs(1:end-1));

% figure
% plot(t,recNcsd7,t(locsNcs),recNcsd7(locsNcs),'*');
% xlabel('Time(sec)')
% ylabel('Amplitude')
% title('NCS - d7')
% grid on; 

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
nFig = 2;

ax1(1) = subplot(nFig,1,1);
plot(ax1(1),t,ncs,':');
hold on
plot(ax1(1),t(locsNcs),ncs(locsNcs),'*','color',[0.3,0.75,0.93]);
plotCute1('Time (sec)','Heartbeat NCS (a.u.)',ax1(1),[],[],0);

%% Find ECG peaks
minPeakHeightEcg = 0.3;
[~,locsEcg] = findpeaks(ecg,'MinPeakDistance',minPeakInterval,...
    'MinPeakHeight',minPeakHeightEcg);
nBeatsEcg = length(locsEcg)-1;
rrIntervalEcg = t(locsEcg(2:end))-t(locsEcg(1:end-1));

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),t,ecg);
hold on
plot(ax1(2),t(locsEcg),ecg(locsEcg),'*','color','k');
plotCute1('Time (sec)','Heartbeat ECG (a.u.)',ax1(2),[],[],0);
linkaxes(ax1(:),'x');

%% Heart Rate (HR) and Heart Rate Variability (HRV) calculation
% Time in seconds to calculate HR/ HRV
tInterval = 5;
% Heart Rate in BPM
hr = zeros(nSample,2); 
% Root Mean Square of Successive Difference between each heartbeat, given
% by: sqrt(mean((RR interval 1 - RR interval 2)^2 + ...)))
rmssd = zeros(nSample,2); 

for sample = tInterval*FsNew:length(ncs)
    % Start and End time for the window used to calculate heart beat
    tWindow = [t(sample-(tInterval*FsNew)+1),t(sample)];
    
    % Finding start & end indices of locs corresponding to first & last
    % peak in the window
    locsWindowNcsIdx = [find(t(locsNcs) >= tWindow(1),1), ...
        find(t(locsNcs) < tWindow(2),1,'last')];
    locsWindowEcgIdx = [find(t(locsEcg) >= tWindow(1),1), ...
        find(t(locsEcg) < tWindow(2),1,'last')];
       
    if (isempty(locsWindowNcsIdx) || isempty(locsWindowEcgIdx))
        fprintf('No peak detected in this period, stopping at sample: %d\n',sample);
        return
    end
    
    nBeatsNcsWindow = (locsWindowNcsIdx(2)-locsWindowNcsIdx(1))-1;
    hr(sample,1) = 60*nBeatsNcsWindow/tInterval; % First column NCS
    
    nBeatsEcgWindow = (locsWindowEcgIdx(2)-locsWindowEcgIdx(1))-1;
    hr(sample,2) = 60*nBeatsEcgWindow/tInterval; % Second column ECG    

    % successive differences between each rr interval
    sdNcs = rrIntervalNcs(locsWindowNcsIdx(1)+1:locsWindowNcsIdx(2)-1) ...
        - rrIntervalNcs(locsWindowNcsIdx(1):locsWindowNcsIdx(2)-2);
    sdEcg = rrIntervalEcg(locsWindowEcgIdx(1)+1:locsWindowEcgIdx(2)-1) ...
        - rrIntervalEcg(locsWindowEcgIdx(1):locsWindowEcgIdx(2)-2);
    rmssd(sample,:) = [sqrt(mean(sdNcs.^2)), sqrt(mean(sdEcg.^2))];

    % In many cases number of peaks are not enough to calculate RMSSD,
    % resulting in NaN. I can also detect NaN and replace by 0, this is
    % just to check 
    if ((locsWindowNcsIdx(2)-locsWindowNcsIdx(1))<3)
        fprintf(['NCS: Number of peaks in interval is less than 3. \n',...
            'Replacing rmssd with 0. \n']);
        rmssd(sample,1) = 0;
    else
        if (locsWindowEcgIdx(2)-locsWindowEcgIdx(1))<3
        fprintf(['ECG: Number of peaks in interval is less than 3. \n',...
            'Replacing rmssd with 0. \n']);
        rmssd(sample,2) = 0;
        end
    end
end

figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
nFig2 = 3;

ax2(1) = subplot(nFig2,1,1);
plot(ax2(1),t,hr(:,1),'color',[0.3,0.75,0.93]);
hold on
plot(ax2(1),t,hr(:,2),'color','k');
plotCute1('Time (sec)','Heartbeat (BPM)',ax2(1),['Data: ',fileName],{'NCS','ECG'},1);

ax2(2) = subplot(nFig2,1,2);
plot(ax2(2),t,rmssd(:,1),'color',[0.3,0.75,0.93]);
hold on
plot(ax2(2),t,rmssd(:,2),'color','k');
plotCute1('Time (sec)','RMSSD (sec)',ax2(2),[],{'NCS','ECG'},1);

ax2(3) = subplot(nFig2,1,3);
plot(ax2(3),t(locsNcs(2:end)),rrIntervalNcs,'color',[0.3,0.75,0.93]);
hold on
plot(ax2(3),t(locsEcg(2:end)),rrIntervalEcg,'color','k');
plotCute1('Time (sec)','RR Interval (sec)',ax2(3),[],{'NCS','ECG'},1);

linkaxes(ax2(:),'x');


