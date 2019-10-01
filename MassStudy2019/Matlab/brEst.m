function [br,pk] = brEst(ncsData,fs,opts)

if ~isfield(opts,'tWinBR')
    opts.tWinBR = 8;
    fprintf('Default BR estimation window: %3.2f\n',opts.tWinBR);
end
if ~isfield(opts,'calibPk')
    opts.calibPk = 0; % Do not calibrate by default
else
    if ~isfield(opts,'calibT')
        opts.calibT = [20,40]; 
        fprintf('Default BR estimation peak calibration window [%d, %d].\n',opts.calibT(1),opts.calibT(2));
    end
    if ~isfield(opts,'calibMinPkRatio')
        opts.calibMinPkRatio = 0.4;
        fprintf('Default min peak height can be %d%% of avg peak height in calibration window.\n',100*opts.calibMinPkRatio);
    end
end

t = (0:(length(ncsData)-1))/fs;

% -------------------------------------------------------------------------
% Minima and maxima detection on NCS thorax and abdomen data
% -------------------------------------------------------------------------
pk = findMaxMin(ncsData,fs,opts);

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

if opts.calibPk == 1
    fprintf('\nBR: Performing calibration of pk-pk height.\n')
    % This is the change in Pk-Pk thorax and abdomen signal
    del2 = zeros(length(pkMax2),1);

    % So for a cycle, considering there exists 2 minima and 2 maxima point:
    % Calculation is peformed using the difference between maxima and first
    % minima. The update is performed at the end of cycle to be consistent with
    % TV from airflow calculation. And this value is held until next update or
    % the end of the waveform.
    for i = 1:length(pkMax2)
            del2(i) = abs(ncsData(pkMax2(i),1)-ncsData(pkMin2(i),1)); % Making it positive always
    end

    tBRmax1 = t(pkMax2);
    idx2Calib = ((tBRmax1 >= opts.calibT(1))&(tBRmax1 <= opts.calibT(2)));
    delNcs2Calib = mean(del2(idx2Calib)); 
%     fprintf('delNCS2Calib = %f\n',delNcs2Calib);

    idxPk = (del2 >= opts.calibMinPkRatio*delNcs2Calib);

    pkMax2 = pkMax2(idxPk); % Only keep indices that satisfy, otherwise ignore
    pkMin2 = pkMin2(idxPk);
    
    pk(1).idxValidPk = idxPk;
end

fig1 = figure;
nFig = size(ncsData,2)+1;

%%
if size(ncsData,2) == 2
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

    if opts.calibPk == 1
        fprintf('BR: Performing calibration of pk-pk height.\n');
        del1 = zeros(length(pkMax1),1);
        for i = 1:length(pkMax1)
                del1(i) = abs(ncsData(pkMax1(i),1)-ncsData(pkMin1(i),1)); % Making it positive always
        end

        tBRmax2 = t(pkMax1);
        idx1Calib = ((tBRmax2 >= opts.calibT(1))&(tBRmax2 <= opts.calibT(2)));
        delNcs1Calib = mean(del1(idx1Calib)); 
%         fprintf('delNCS1Calib = %f\n',delNcs1Calib);

        idxPk = (del1 >= opts.calibMinPkRatio*delNcs1Calib);

        pkMax1 = pkMax1(idxPk); % Only keep indices that satisfy, otherwise ignore
        pkMin1 = pkMin1(idxPk);
        pk(2).idxValidPk = idxPk;
    end
    %% Find breath rate in last tWinBR sec (or slightly less)
    br = zeros(length(t),2);
    
    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinBR))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            br(iter,1) = 0;
        else
            br(iter,1) = (idxBR(2)- idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
    %         ncsAmpBR(iter) = (idxAmpBR(2)- idxAmpBR(1))/tWinBR;
        end

        idxBR2 = find((tBRmax2 > (t(iter)-opts.tWinBR))&(tBRmax2 <= t(iter)));

        if (size(idxBR2,1) <= 1) 
            br(iter,2) = 0;
        else
            br(iter,2) = (idxBR2(2)- idxBR2(1))/(tBRmax2(idxBR2(2))-tBRmax2(idxBR2(1)));
    %         ncsPhBR(iter) = (idxPhBR(2)- idxPhBR(1))/tWinBR;

        end

    end

    br = 60.*br;

    figure(fig1);
    ax(1) = subplot(nFig,1,1);
    plot(t,ncsData(:,1)); hold on;
    plot(t(pkMax2),ncsData(pkMax2,1),'^',...
        t(pkMin2),ncsData(pkMin2,1),'v');
    leg = {'Th','Max','Min'};
    plotCute1('Time (s)','Resp (mV)',ax(1),[],leg,1,'Horizontal');
    ax(2) = subplot(nFig,1,2);
    plot(t,ncsData(:,2)); hold on;
    plot(t(pkMax1),ncsData(pkMax1,2),'^',...
        t(pkMin1),ncsData(pkMin1,2),'v');
    leg = {'Abd','Max','Min'};
    plotCute1('Time (s)','Resp (mV)',ax(2),[],leg,1,'Horizontal');
    ax(3) = subplot(nFig,1,3);
    plot(t,br(:,1)); hold on
    plot(t,br(:,2));
    leg = {'Th BR','Abd BR'};
    plotCute1('Time (s)','BR (BPM)',ax(3),[],leg,1,'Horizontal');
    linkaxes(ax,'x');


else
%% Only one NCS column    
    % Find breath rate in last tWinBR sec (or slightly less)
    br = zeros(length(t),1);
    t = t(:);
    tBRmax1 = t(pkMax2);
    
    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinBR))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            br(iter,1) = 0;
        else
            br(iter,1) = (idxBR(2)- idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
    %         ncsAmpBR(iter) = (idxAmpBR(2)- idxAmpBR(1))/tWinBR; % This is wrong
  
        end

    end

    br = 60.*br;

    
    figure(fig1);
    ax(1) = subplot(nFig,1,1);
    plot(t,ncsData); hold on;
    plot(t(pkMax2),ncsData(pkMax2,1),'^',...
        t(pkMin2),ncsData(pkMin2,1),'v');
    leg = {'Resp','Max','Min'};
    plotCute1('Time (s)','Resp (mV)',ax(1),[],leg,1,'Horizontal');
    ax(2) = subplot(nFig,1,2);
    plot(t,br(:,1)); 
    leg = {'BR'};
    plotCute1('Time (s)','BR (BPM)',ax(2),[],leg,1,'Horizontal');
    linkaxes(ax,'x');
    
end

