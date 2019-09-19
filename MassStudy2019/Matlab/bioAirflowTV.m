% This function estimates Tidal volume for each breath (at the end of
% expiration) from the BIOPAC airflow data.
% May 04, 2019
% Pragya Sharma, ps847@cornell.edu

function [tv,vol, fig] = bioAirflowTV(airflowBio,fs,meanDev)
% Input: 
% airflowBio: BIOPAC airflow data in L/s
% fs: Sample rate of BIOPAC in Hz
% meanDev: Baseline correction of airflow
% -------------------------------------------------------------------------
t = (0:(length(airflowBio)-1))/fs;

% Estimating volume from airflow by integration over each cycle. 
[vol,zcPosNeg,airflowFilt] = airflowToVol(airflowBio,fs,meanDev);

% Estimating TV by calculating (vol_max - vol_min) for each breath
tv = zeros(length(t),1);
% tTV = zeros(length(t),1);
count = 1;
for i = 1:2:length(zcPosNeg)-2
    volMax = max(vol(zcPosNeg(i,1):zcPosNeg(i+2,1)));
    volMin = min(vol(zcPosNeg(i,1):zcPosNeg(i+2,1)));
    if (i+4)<=length(zcPosNeg)
        tv(zcPosNeg(i+2,1):(zcPosNeg(i+4,1)-1)) = volMax-volMin;
    else
        tv(zcPosNeg(i+2,1):end) = volMax-volMin;
    end
    % tv(count) = volMax-volMin;
    % tTV = t(zcPosNeg(i+2,1));
    count = count + 1;
end

fig = figure;
nFig = 2;
ax(1) = subplot(nFig,1,1);
plot(t,airflowBio);hold on;
plot(t,airflowFilt);
idxPos = zcPosNeg(:,2) > 0;
idxNeg = zcPosNeg(:,2) <= 0;
plot(t(zcPosNeg(idxPos,1)),airflowFilt(zcPosNeg(idxPos,1)),'^','MarkerSize',6,'color',[0 0 0],'MarkerFaceColor',[0 0 0]);
plot(t(zcPosNeg(idxNeg,1)),airflowFilt(zcPosNeg(idxNeg,1)),'v','MarkerSize',6,'color',[0.3 0.7 0.3],'MarkerFaceColor',[0.3 0.7 0.3]);
leg = {'airflow','airflow filtered','Inspire Start','Expire Start'};
plotCute1('Time (s)','Airflow (L/s)',ax(1),[],leg,1,'Horizontal');

ax(2) = subplot(nFig,1,2);
plot(t,vol);
hold on
% stem(tTV,tv,'filled');
plot(t,tv);
leg = {'Volume','TV'};
plotCute1('Time (s)','Volume (L)',ax(2),[],leg,1);

linkaxes(ax,'x');

end