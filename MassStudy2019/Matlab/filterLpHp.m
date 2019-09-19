function dataFilt = filterLpHp(data,fs,opts)
% filterLpHp function is for filtering any signal. It applies
% high pass and/or low pass filter to get dc-removed signal. If data has
% multiple channels as columns, it filters all of them separately.
% -------------------------------------------------------------------------
% data: Input data, each data is arranged in a column
% fs: sampling frequency of data
% Specify options: opts._:
% opts.filtType: Low pass 'Lp' or High pass 'Hp', or both 'LpHp'.
% opts.f3db: 3dB cutoff frequency of butterworth iir hp filter.
% fp, fst: passband and stopband frequency of kaiserwin fir lp filter
% signData: Takes care of polarity reversal for amplitude and phase with
%            sign arranged as [ampSign, phSign]
% Example: 
% fs = 500;
% opts.filtType = 'LpHp'; 
% opts.f3dbHP = 0.5;
% opts.fpLP = 5;
% opts.fstLP = 7;
% [dataFilt,~] = filterLpHp(data,opts);

%% Checking input size
[dataRow, dataCol] = size(data);
if (dataRow < dataCol)
    fprintf(['filterLpHp: Each channel should be a column. \nDoesn''t ',...
        'satisfy shape constraints. Stopping ...\n ']);
    return;
end    
% Specifying default opts
if ~isfield(opts,'filtType')
    opts.filtType = 'Lp'; % Only low pass filtering by default
    fprintf('No filtType specified, using low pass filter.\n');
end
if strcmp(opts.filtType,'Lp') || strcmp(opts.filtType,'LpHp')
    if ~isfield(opts,'fpLP') || ~isfield(opts,'fstLP')
        fprintf('No fpLP, fstLP specified, stopping ...\n'); return;
    end
end
if strcmp(opts.filtType,'Hp') || strcmp(opts.filtType,'LpHp')
    if ~isfield(opts,'f3db')
        fprintf('No f3dbHp specified, stopping ...\n');return;
    end
end
if ~isfield(opts,'Ast')
    % Only needed when filttype is 'Lp' or 'LpHp'
    opts.Ast = 10; % 10 dB attenuation by default 
    % fprintf('No Ast specified, using default %3.2f.\n',opts.Ast);
end
if ~isfield(opts,'orderHP')
    % Only needed when filttype is 'Hp' or 'LpHp'
    opts.orderHP = 3; % Order of HP filter by default
    % fprintf('No HP filter order specified, using %d.\n',opts.orderHP);
end

%% 
t = (0:(length(data)-1))/fs;

dataFilt = data;

if strcmp(opts.filtType,'Hp') || strcmp(opts.filtType,'LpHp')
    % High Pass filter
    filtHP = fdesign.highpass('N,F3db',opts.orderHP,opts.f3db,fs); % Nonlinear phase filter but better filtering characteristics
    HdHP = design(filtHP,'butter'); % 'butter' with 'N,F3db' specifications
    dataFilt = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,dataFilt);
end
if strcmp(opts.filtType,'Lp') || strcmp(opts.filtType,'LpHp')
    % Low pass filter 
    filtLP = fdesign.lowpass('Fp,Fst,Ap,Ast',opts.fpLP,opts.fstLP,0.1,opts.Ast,fs); % Change Ast from 80 to 10 if filtfilt creates error
    HdLP = design(filtLP,'kaiserwin'); %kaiserwin is much faster than equiripple without much degradation (?)
    dataFilt = filtfilt(HdLP.Numerator,1,dataFilt); % Taking care of group delay
end

figure
nFig = dataCol;
for i = 1:nFig
    ax(i) = subplot(nFig,1,i); %#ok<AGROW>
    plot(t,dataFilt(:,i));
end
linkaxes(ax,'x');
end
