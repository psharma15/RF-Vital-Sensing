% This function estimates heart rate from ECG (BIOPAC)
% Pragya Sharma, ps847@cornell.edu

function [hr,locPk] = ecgHR(data,fs,opts)

% ECG should be filtered ECG (f3db hp: 4 Hz, lp: [20,25] Hz) to get good
% peak detection using this code

if ~isfield(opts,'minPkHt')
    opts.minPkHt = 0.5;
    fprintf('Default ECG min peak height %3.2f V\n',opts.minPkHt);
end
if ~isfield(opts,'maxHR')
    opts.maxHR = 180; % Max HR in BPM
    fprintf('Default ECG max heart rate %d BPM\n',opts.maxHR);
end    
if ~isfield(opts,'tWinHR')
    opts.tWinHR = 5;
    fprintf('Default HR estimation window: %3.2f\n',opts.tWinHR);
end

t = (0:(length(data)-1))/fs;

minPkPkSamp = fs*(1/(opts.maxHR/60)); % Minimum number of samples between two peaks

% Simple findpeaks for peak detection
[~,locPk] = findpeaks(data,'MinPeakHeight',opts.minPkHt,'MinPeakDistance',minPkPkSamp);

figure
hold on;
plot(t,data);
plot(t(locPk),data(locPk),'^');

%% Find hreath rate in last tWinHR sec (or slightly less)
hr = zeros(length(t),1);
t = t(:);
tHRmax = t(locPk);

for iter = 1:length(t)

    idxHR = find((tHRmax > (t(iter)-opts.tWinHR))&(tHRmax <= t(iter)));

    if (size(idxHR,1) <= 1)
        hr(iter,1) = 0;
    else
        hr(iter,1) = (idxHR(2)- idxHR(1))/(tHRmax(idxHR(2))-tHRmax(idxHR(1)));
%         ncsAmpHR(iter) = (idxAmpHR(2)- idxAmpHR(1))/tWinHR;
    end

end

hr = 60.*hr;


end