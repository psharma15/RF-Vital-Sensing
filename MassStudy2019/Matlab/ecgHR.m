% This function estimates heart rate from ECG (BIOPAC)
% Pragya Sharma, ps847@cornell.edu

function [hr,locPk] = ecgHR(data,fs,opts)

if ~isfield(opts,'minPkHt')
    opts.minPkHt = 0.5;
    fprintf('Default ECG min peak height %3.2f\n',opts.minPkHt);
end
if ~isfield(opts,'tWinHR')
    opts.tWinHR = 5;
    fprintf('Default HR estimation window: %3.2f\n',opts.tWinHR);
end


t = (0:(length(data)-1))/fs;

% Simple findpeaks for peak detection
[~,locPk] = findpeaks(data,'MinPeakHeight',opts.minPkHt);

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