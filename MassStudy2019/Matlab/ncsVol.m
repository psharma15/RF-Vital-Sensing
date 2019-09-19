% This function uses calibration coefficient from belt/airflow TV to 
% estimate NCS lung volume.
% May 08, 2019
% Pragya Sharma, ps847@cornell.edu

function vNcs = ncsVol(ncsData,calibCoeff,fs,opts)
% Input:
% ncsData: Ncs data in column format [Th Amp, Abd Amp] @fs
% calibCoeff: Calibration coefficients, fit equation in opts
% fs: Sampling frequency of the above data in Hz
% opts: Additional options related to calibration, peak detection etc.
% Output:
% tvNcs: Calibrated NCS tidal volume
%
% -------------------------------------------------------------------------

t = (0:(length(ncsData)-1))/fs;

switch (opts.calibType)
    
    case 'tv'

        % -------------------------------------------------------------------------
        % Minima and maxima detection on NCS thorax and abdomen data
        % -------------------------------------------------------------------------
        ncsPk = findMaxMin(ncsData,fs,opts);

        if ncsPk(1).ind(1) == 1
            % If first peak is inhalation peak, skip it
            ncsPk(1).ind = ncsPk(1).ind(2:end); 
            ncsPk(1).idx = ncsPk(1).idx(2:end);
        end

        if ncsPk(1).ind(end) == 1
            ncsPk(1).ind = ncsPk(1).ind(1:end-1);
            ncsPk(1).idx = ncsPk(1).idx(1:end-1);
        end


        pkMax2 = ncsPk(1).idx(ncsPk(1).ind == 1);
        pkMin2 = ncsPk(1).idx(ncsPk(1).ind == 0);

        % This is the change in Pk-Pk thorax and abdomen signal
        del2 = zeros(length(ncsData),1);

        % So for a cycle, considering there exists 2 minima and 2 maxima point:
        % Calculation is peformed using the difference between maxima and first
        % minima. The update is performed at the end of cycle to be consistent with
        % TV from airflow calculation. And this value is held until next update or
        % the end of the waveform.
        for i = 1:length(pkMax2)
            if i<length(pkMax2)
                del2(pkMin2(i+1):pkMin2(i+2)) = abs(ncsData(pkMax2(i),1)-ncsData(pkMin2(i),1)); % Making it positive always
            else
                del2(pkMin2(i+1):end) = abs(ncsData(pkMax2(i),1)-ncsData(pkMin2(i),1));
            end
        end

        fig1 = figure;
        nFig = size(ncsData,2);
        if size(ncsData,2) == 2
            if ncsPk(2).ind(1) == 1
                ncsPk(2).ind = ncsPk(2).ind(2:end);
                ncsPk(2).idx = ncsPk(2).idx(2:end);
            end

            if ncsPk(2).ind(end) == 1
                ncsPk(2).ind = ncsPk(2).ind(1:end-1);
                ncsPk(2).idx = ncsPk(2).idx(1:end-1);
            end

            pkMax1 = ncsPk(2).idx(ncsPk(2).ind == 1);
            pkMin1 = ncsPk(2).idx(ncsPk(2).ind == 0);

            del1 = zeros(length(ncsData),1);

            for i = 1:length(pkMax1)
                if i<length(pkMax1)
                    del1(pkMin1(i+1):pkMin1(i+2)) = abs(ncsData(pkMax1(i),2)-ncsData(pkMin1(i),2)); % Making it positive always
                else
                    del1(pkMin1(i+1):end) = abs(ncsData(pkMax1(i),2)-ncsData(pkMin1(i),2));
                end
            end

            figure(fig1);
            ax(1) = subplot(nFig,1,1);
            plot(t,ncsData(:,1)); hold on;
            plot(t(pkMax2),ncsData(pkMax2,1),'^',...
                t(pkMin2),ncsData(pkMin2,1),'v');
            leg = {'Th NCS','Max','Min'};
            plotCute1('Time (s)','NCS (mV)',ax(1),[],leg,1,'Horizontal');
            ax(2) = subplot(nFig,1,2);
            plot(t,ncsData(:,2)); hold on;
            plot(t(pkMax1),ncsData(pkMax1,2),'^',...
                t(pkMin1),ncsData(pkMin1,2),'v');
            leg = {'Abd NCS','Max','Min'};
            plotCute1('Time (s)','NCS (mV)',ax(2),[],leg,1,'Horizontal');
            linkaxes(ax,'x');

        else
            figure(fig1);
            ax(1) = subplot(nFig,1,1);
            plot(t,ncsData); hold on;
            plot(t(pkMax2),ncsData(pkMax2,1),'^',...
                t(pkMin2),ncsData(pkMin2,1),'v');
            leg = {'NCS','Max','Min'};
            plotCute1('Time (s)','NCS (mV)',ax(1),[],leg,1,'Horizontal');
        end

        % -------------------------------------------------------------------------
        % Perform calibration now using available coefficients
        % -------------------------------------------------------------------------
        % Select calibration type
        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del1,del2);
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    del1,del2);
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),del1,del2);
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del1);
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del2);
            case 'AbdQuad'
                fitFunc = @(A,B,C,x) A.*(x.^2)+ B.*x + C;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),del2);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end
        
    case 'vol'
        % Calibrate directly the respiration waveforms to get instantaneous
        % lung volume due to inhalation and exhalation. Does not include
        % the default lung volume information.
        ncs1 = ncsData(:,1); 
        if size(ncsData,2) == 2
            ncs2 = ncsData(:,2);
        end

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1,ncs2);
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    ncs1,ncs2);
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),ncs1,ncs2);
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1);
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1);
            case 'AbdQuad'
                fitFunc = @(A,B,C,x) A.*(x.^2)+ B.*x + C;
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),ncs1);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        fig1 = figure;

        if size(ncsData,2) == 2
            figure(fig1);
            ax1(1) = subplot(2,1,1);
            plot(t,ncsData(:,1),t,ncsData(:,2));
            leg = {'NCS Th','NCS Abd'};
            plotCute1('Time (s)','mV',ax1(1),[],leg,1);

        else
            figure(fig1);
            ax1(1) = subplot(2,1,1);
            plot(t,ncsData(:,1));
            leg = {'NCS'};
            plotCute1('Time (s)','mV',ax1(1),[],leg,1);
        end
        ax1(2) = subplot(2,1,2);
        plot(t,vNcs)
        leg = {'Calibrated NCS'};
        plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        linkaxes(ax1,'x');

end