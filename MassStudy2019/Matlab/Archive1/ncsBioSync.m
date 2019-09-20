% This code synchronizes NCS with BIOPAC. NCS file is a subsection of
% BIOPAC data. 
% May 02 2019
% Pragya Sharma, ps847@cornell.edu

function [ncs,bio,mic,fs,tNcs,tBio,fig] = ...
    ncsBioSync(dataPath,ncsFile,bioFile,fsNcs,fsBio,tStartOff,tEndOff,...
    tManualOff,ifDownSamp,fsDS,signNcs,unwrapPh)
% -------------------------------------------------------------------------
% Inputs:
% dataPath: Path for NCS and BIOPAC data
% ncsFile: Name of NCS file
% bioFile: Name of BIOPAC file
% fsNcs: Sampling frequency of NCS in Hz. Usual is 50 kHz
% fsBio: Sampling frequency of BIOPAC in Hz. Usual is 2 kHz
% tStartOff: Start time offset to NCS waveform (to truncate NCS waveform)
% tEndOff: End time offset to NCS waveform (to truncate NCS waveform)
% tManualOFf: If NCS start time or BIOPAC start times are not accurate
% and needs any manual adjustment. Or difference only due to missing 
% fractional seconds information in BIOPAC time.
% -------------------------------------------------------------------------
% Outputs:
% ncs: ncs data in format [amp_ch1 ph_ch1 amp_ch2 ph_ch2]
% bio: biopac data in format [ch1, ch2, ch3, ch4] = [ecg, th_belt,
% abd_belt, pneumotach]
% mic:  MIC data.
% fsNcs, fsBio: Sampling frequency of NCS and BIOPAC in Hz. Same as input
% if no downsamping is performed.
% fs: [ncs,bio] Sampling frequencies in Hz, different from input if
% downsampled
% -------------------------------------------------------------------------

% Reading NCS and BIOPAC data
[ncsUnsync,mic, tNcsStart, ~] = readNcsData(dataPath,ncsFile);
[bioUnsync,tBioStart] = readBioData(dataPath,bioFile);

% If unwrap phase, do it now
if unwrapPh(1) == 1
    ncsUnsync(:,2) = rad2deg(unwrap(deg2rad(ncsUnsync(:,2))));
end
if unwrapPh(2) == 1
    ncsUnsync(:,4) = rad2deg(unwrap(deg2rad(ncsUnsync(:,4))));
end
% Optional downsample of NCS and/or BIOPAC to input sampling rate
if ifDownSamp(1) == 1 
    ncsUnsync = resample(ncsUnsync,fsDS(1),fsNcs);
    fs(1) = fsDS(1);
else 
    fs(1) = fsNcs;
end
if ifDownSamp(2) == 1
    bioUnsync = resample(bioUnsync,fsDS(2),fsBio);
    fs(2) = fsDS(2);
else
    fs(2) = fsBio;
end

% Finding starting index for NCS and BIOPAC
tRelNcsStart = tNcsStart + [0 0 0 0 0 tStartOff];
tRelBioStart = tRelNcsStart + [0 0 0 0 0 tManualOff]; 

tDiffNcs = etime(tRelNcsStart,tNcsStart);
tDiffBio = etime(tRelBioStart,tBioStart);

ncsStartSamp = int64(tDiffNcs*fs(1)+1);
bioStartSamp = int64(tDiffBio*fs(2)+1);

% Finding stopping index for NCS and BIOPAC
tNcsStop = tNcsStart + [0 0 0 0 0 1].*(length(ncsUnsync)/fs(1));
tRelNcsStop = tNcsStop - [ 0 0 0 0 0 tEndOff];
tBioStop = tRelNcsStop + [0 0 0 0 0 tManualOff]; 

tDiffNcs = etime(tRelNcsStop,tNcsStart);
tDiffBio = etime(tBioStop,tBioStart);

ncsStopSamp = int64(tDiffNcs*fs(1));
bioStopSamp = int64(tDiffBio*fs(2));

ncs = ncsUnsync(ncsStartSamp:ncsStopSamp,:);
if size(mic,1) ~=0
    mic = mic(ncsStartSamp:ncsStopSamp);
end
bio = bioUnsync(bioStartSamp:bioStopSamp,:);

% Sign of NCS
ncs = [signNcs(1).*ncs(:,1), signNcs(2).*ncs(:,2), signNcs(3).*ncs(:,3),...
       signNcs(4).*ncs(:,4)];
    
% Plotting figure
tNcs = ((0:(length(ncs)-1))/fs(1))';
tBio = ((0:(length(bio)-1))/fs(2))';

fig = figure('position',[100,100,900,800]);
nFig = 4;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(tNcs,ncs(:,1));xlabel('Time (s)'); ylabel('NCS Th Amp');grid on;
yyaxis right
plot(tNcs,ncs(:,2));  xlabel('Time (s)'); ylabel('NCS Th Ph');
% ylim([-200,200]);

ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(tNcs,ncs(:,3));xlabel('Time (s)'); ylabel('NCS Abd Amp');grid on;
yyaxis right
plot(tNcs,ncs(:,4)); xlabel('Time (s)'); ylabel('NCS Abd Ph');
% ylim([-200,200]);

ax1(3) = subplot(nFig,1,3);
plot(tBio,bio(:,2)); hold on;
plot(tBio,bio(:,3)); grid on;
xlabel('Time (s)'); ylabel('BIOPAC resp (mV)');
legend('Thorax','Abdomen'); hold off;
ax1(4) = subplot(nFig,1,4);
plot(tBio,bio(:,4));grid on;
xlabel('Time (s)'); ylabel('BIOPAC airflow (L/s)');
linkaxes(ax1,'x');

end