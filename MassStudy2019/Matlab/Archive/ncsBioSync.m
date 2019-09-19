% This code synchronizes NCS with BIOPAC. NCS file is a subsection of
% BIOPAC data. 
% May 02 2019
% Pragya Sharma, ps847@cornell.edu

function [ncs,bio,mic] = ncsBioSync(dataPath,ncsFile,bioFile,fsNcs,fsBio,tStartOff,tEndOff,tManualOff)
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

[ncsUnsync,mic, tNcsStart, ~] = readNcsData(dataPath,ncsFile);
[bioUnsync,tBioStart] = readBioData(dataPath,bioFile);

% Finding starting index for NCS and BIOPAC
tRelNcsStart = tNcsStart + [0 0 0 0 0 tStartOff];
tRelBioStart = tRelNcsStart + [0 0 0 0 0 tManualOff]; 

tDiffNcs = etime(tRelNcsStart,tNcsStart);
tDiffBio = etime(tRelBioStart,tBioStart);

ncsStartSamp = int64(tDiffNcs*fsNcs+1);
bioStartSamp = int64(tDiffBio*fsBio+1);

% Finding stopping index for NCS and BIOPAC
tNcsStop = tNcsStart + [0 0 0 0 0 1].*(length(ncsUnsync)/fsNcs);
tRelNcsStop = tNcsStop - [ 0 0 0 0 0 tEndOff];
tBioStop = tRelNcsStop + [0 0 0 0 0 tManualOff]; 

tDiffNcs = etime(tRelNcsStop,tNcsStart);
tDiffBio = etime(tBioStop,tBioStart);

ncsStopSamp = int64(tDiffNcs*fsNcs);
bioStopSamp = int64(tDiffBio*fsBio);

ncs = ncsUnsync(ncsStartSamp:ncsStopSamp,:);
mic = mic(ncsStartSamp:ncsStopSamp);
bio = bioUnsync(bioStartSamp:bioStopSamp,:);

end