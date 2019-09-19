% This code converts airflow, F (L/s) from BIOPAC pneumotach transducer to
% volume, V (L). V is integration of F over time.
% May 03, 2019
% Pragya Sharma, ps847@cornell.edu
% -------------------------------------------------------------------------
function [vol,zcPosNeg,airflowFilt] = airflowToVol(airflowBio,fs,meanDev)
t = (0:(length(airflowBio)-1))/fs;

% Airflow baseline correction
airflowBio = airflowBio - meanDev;

% Filtering for true zero-crossing detection
opts.filtType = 'Lp';
opts.fpLP = 2; opts.fstLP = 2.5;
airflowFilt = filterLpHp(airflowBio,fs,opts);

% Minimum time between two zero crossing points
opts.minTime = 0.3; 

% Gives indices and pos/neg slope indication
zcPosNeg = zeroCrossDet(airflowFilt,fs,opts);

% Starting integration from inspire beginning only. and making sure ends
% with inspire beginning as well.
if zcPosNeg(1,2) == -1
    zcPosNeg = zcPosNeg(2:end,:);
end
if zcPosNeg(end,2) == -1
    zcPosNeg = zcPosNeg(1:end-1,:);
end

% Now find volume by integrating for each cycle
vol = zeros(length(airflowBio),1);

for i = 1:2:(length(zcPosNeg))-2
    idxStart = zcPosNeg(i,1);
    idxEnd = zcPosNeg(i+2,1);
    vol(idxStart:idxEnd-1) = cumtrapz(t(idxStart:idxEnd-1),airflowBio(idxStart:idxEnd-1));
end


end