function [ncsFilt,ncsUnfiltAmpTrunc,ncsUnfiltPhTrunc,tTrunc] = postProcess(f3db,fp,fst,data,fs,signData,varargin)
%POSTPROCESS function is for post processing NCS signal. It applies
%high pass, low pass filter to get dc-removed, respiratory, or heartbeat 
%signal.
% Example: [ncsAmpFilt,~,~,~,t] = postProcess(0.9,10,11,ncs,fs,[-1,-1],0,30);
% where, f3db = 3dB cutoff frequency of butterworth iir hp filter, it can
% be input as 0, if no HP filter needed.
% fp, fst = passband and stopband frequency of kaiserwin fir lp filter
% data = Amplitude and phase column vectors arranged as [ncsAmp, ncsPh]
% fs = sampling frequency of data
% signData = Takes care of polarity reversal for amplitude and phase with
%            sign arranged as [ampSign, phSign]
% varargin{1} = sampStart: Start time(in seconds) of truncated data
% varargin{2} = sampEnd: End time (in seconds) of truncated data

%%
[dataRow, dataCol] = size(data);
if (dataRow < 1000)
    fprintf('In post-processing, each column should have minimum data points ');
    return;
end
    
nSample = length(data(:,1));
idx = 1:nSample;
t = ((idx-1)/fs)';

%%
sampStart = 0*fs + 1;
sampEnd = nSample;

if nargin > 6
    sampStart = varargin{1}*fs+1;
    if nargin == 8
        sampEnd = varargin{2}*fs+1;
    end
end

%%
ncsUnfiltAmpTrunc = signData(1)*data(sampStart:sampEnd,1);
ncsUnfiltPhTrunc = signData(2)*data(sampStart:sampEnd,2);

% ncsUnfiltMagTrunc = ncsUnfiltAmpTrunc.*(cos(ncsUnfiltPhTrunc) + 1i*sin(ncsUnfiltPhTrunc));
tTrunc = t(sampStart:sampEnd);

% High Pass filter - for dc or dc and breath removal
if f3db > 0
    % Option to perform HP filter.
    orderHP = 3; % 20 original
    filtHP = fdesign.highpass('N,F3db',orderHP,f3db,fs); % Nonlinear phase filter but better filtering characteristics
    HdHP = design(filtHP,'butter'); % 'butter' with 'N,F3db' specifications
    ncsAmpHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltAmpTrunc);
    ncsPhHP = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,ncsUnfiltPhTrunc);
else
    ncsAmpHP = ncsUnfiltAmpTrunc;
    ncsPhHP = ncsUnfiltPhTrunc;
end

% Low pass filter 
filtLP = fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fst,0.1,10,fs); % Change Ast from 80 to 10 if filtfilt creates error
HdLP = design(filtLP,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
ncsAmpLP = filtfilt(HdLP.Numerator,1,ncsAmpHP); % Taking care of group delay
ncsPhLP = filtfilt(HdLP.Numerator,1,ncsPhHP); % Taking care of group delay

ncsFilt = [ncsAmpLP, ncsPhLP];

figure
fontSize = 14;
yyaxis left
plot(tTrunc,ncsFilt(:,1),'markers',24)
xlabel('Time (sec)','FontSize',fontSize)
ylabel('Filtered Amplitude','FontSize',fontSize)
yyaxis right
plot(tTrunc,unwrap(ncsFilt(:,2)),'markers',24)
ylabel('Filtered Phase','FontSize',fontSize)
grid on;
end
