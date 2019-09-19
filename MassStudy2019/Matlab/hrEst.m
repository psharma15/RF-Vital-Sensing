% This function estimates NCS heart rate by peak finding algorithm
% Pragya sharma, ps847@cornell.edu
% 17 June 2019

function [hr,pk] = hrEst(data,fs,opts)
% Input: 
% ncsData: NCS data arranged in column format (1 or 2 columns corresponding
% to different channels). 
% fs: Sampling frequency
% opts: Options related to peak detection.


if ~isfield(opts,'tWinHR')
    opts.tWinHR = 5;
    fprintf('Default HR estimation window: %3.2f\n',opts.tWinHR);
end

t = (0:(length(data)-1))/fs;

% -------------------------------------------------------------------------
% Minima and maxima detection on NCS thorax and abdomen data
% -------------------------------------------------------------------------
pk = findMaxMin(data,fs,opts);

if pk(1).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 1
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end


pkMax2 = pk(1).idx(pk(1).ind == 1);
pkMin2 = pk(1).idx(pk(1).ind == 0);

% This is the change in Pk-Pk signal, could be used for rejecting a peak
del2 = zeros(length(data),1);

% So for a cycle, considering there exists 2 minima and 2 maxima point:
% Calculation is peformed using the difference between maxima and first
% minima. 
for i = 1:length(pkMax2)
    if i<length(pkMax2)
        del2(pkMin2(i+1):pkMin2(i+2)) = abs(data(pkMax2(i),1)-data(pkMin2(i),1)); % Making it positive always
    else
        del2(pkMin2(i+1):end) = abs(data(pkMax2(i),1)-data(pkMin2(i),1));
    end
end

fig1 = figure;
nFig = size(data,2)+1;
if size(data,2) == 2
    if pk(2).ind(1) == 1
        pk(2).ind = pk(2).ind(2:end);
        pk(2).idx = pk(2).idx(2:end);
    end

    if pk(2).ind(end) == 1
        pk(2).ind = pk(2).ind(1:end-1);
        pk(2).idx = pk(2).idx(1:end-1);
    end

    pkMax1 = pk(2).idx(pk(2).ind == 1);
    pkMin1 = pk(2).idx(pk(2).ind == 0);

    del1 = zeros(length(data),1);

    for i = 1:length(pkMax1)
        if i<length(pkMax1)
            del1(pkMin1(i+1):pkMin1(i+2)) = abs(data(pkMax1(i),2)-data(pkMin1(i),2)); % Making it positive always
        else
            del1(pkMin1(i+1):end) = abs(data(pkMax1(i),2)-data(pkMin1(i),2));
        end
    end

    %% Find heart rate in last tWinHR sec (or slightly less)
    hr = zeros(length(t),2);
    tHRmax1 = t(pkMax2);
    tHRmax2 = t(pkMax1);
    
    for iter = 1:length(t)

        idxHR = find((tHRmax1 > (t(iter)-opts.tWinHR))&(tHRmax1 <= t(iter)));

        if (size(idxHR,1) <= 1)
            hr(iter,1) = 0;
        else
            hr(iter,1) = (idxHR(2)- idxHR(1))/(tHRmax1(idxHR(2))-tHRmax1(idxHR(1)));
    %         ncsAmpHR(iter) = (idxAmpHR(2)- idxAmpHR(1))/tWinHR;
        end

        idxHR2 = find((tHRmax2 > (t(iter)-opts.tWinHR))&(tHRmax2 <= t(iter)));

        if (size(idxHR2,1) <= 1) 
            hr(iter,2) = 0;
        else
            hr(iter,2) = (idxHR2(2)- idxHR2(1))/(tHRmax2(idxHR2(2))-tHRmax2(idxHR2(1)));
    %         ncsPhHR(iter) = (idxPhHR(2)- idxPhHR(1))/tWinHR;

        end

    end

    hr = 60.*hr;

    figure(fig1);
    ax(1) = subplot(nFig,1,1);
    plot(t,data(:,1)); hold on;
    plot(t(pkMax2),data(pkMax2,1),'^',...
        t(pkMin2),data(pkMin2,1),'v');
    leg = {'Th NCS','Max','Min'};
    plotCute1('Time (s)','NCS (mV)',ax(1),[],leg,1,'Horizontal');
    ax(2) = subplot(nFig,1,2);
    plot(t,data(:,2)); hold on;
    plot(t(pkMax1),data(pkMax1,2),'^',...
        t(pkMin1),data(pkMin1,2),'v');
    leg = {'Abd NCS','Max','Min'};
    plotCute1('Time (s)','NCS (mV)',ax(2),[],leg,1,'Horizontal');
    ax(3) = subplot(nFig,1,3);
    plot(t,hr(:,1)); hold on
    plot(t,hr(:,2));
    leg = {'Th HR','Abd HR'};
    plotCute1('Time (s)','HR (BPM)',ax(3),[],leg,1,'Horizontal');
    linkaxes(ax,'x');


else
    
    %% Find hreath rate in last tWinHR sec (or slightly less)
    hr = zeros(length(t),1);
    t = t(:);
    tHRmax1 = t(pkMax2);
    
    for iter = 1:length(t)

        idxHR = find((tHRmax1 > (t(iter)-opts.tWinHR))&(tHRmax1 <= t(iter)));

        if (size(idxHR,1) <= 1)
            hr(iter,1) = 0;
        else
            hr(iter,1) = (idxHR(2)- idxHR(1))/(tHRmax1(idxHR(2))-tHRmax1(idxHR(1)));
    %         ncsAmpHR(iter) = (idxAmpHR(2)- idxAmpHR(1))/tWinHR;
        end

    end

    hr = 60.*hr;

    
    figure(fig1);
    ax(1) = subplot(nFig,1,1);
    plot(t,data); hold on;
    plot(t(pkMax2),data(pkMax2,1),'^',...
        t(pkMin2),data(pkMin2,1),'v');
    leg = {'NCS','Max','Min'};
    plotCute1('Time (s)','NCS (mV)',ax(1),[],leg,1,'Horizontal');
    ax(2) = subplot(nFig,1,2);
    plot(t,hr(:,1)); 
    leg = {'HR'};
    plotCute1('Time (s)','HR (BPM)',ax(2),[],leg,1,'Horizontal');
    linkaxes(ax,'x');
    
end

