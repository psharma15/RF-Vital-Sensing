% This function calibrates NCS using biopac belt/airflow meter (pneumotach)
% May 08, 2019
% Pragya Sharma, ps847@cornell.edu

function [calibCoeff,vNcs] = ncsVolCalib(ncsData,fs,tvAirflow,volAirflow,opts)
% Input:
% data: Ncs data in column format [Th Amp, Abd Amp] @fs
% fs: Sampling frequency of the above data in Hz
% opts: Additional options related to calibration, peak detection etc.
% Output:
% calibCoeff: Calibration coefficient of NCS
% tvNcs: Calculated TV from NCS after calibration
% ncsPk: Peaks of NCS
% -------------------------------------------------------------------------
if ~isfield(opts,'fitEqn')
    opts.fitEqn = 'BiasedLinear';    
    fprintf('Default calibration fit equation used: %s\n',opts.fitEqn);
end
if ~isfield(opts,'calibType')
    opts.calibType = 'vol';    
    fprintf('Default calibration type used: %s\n',opts.calibType);
end
if ~isfield(opts,'tCalib')
    opts.tCalib = [0, (length(bioData)-1)/fs];
end

t = (0:(length(ncsData)-1))/fs;

switch(opts.calibType)
    
    case 'tv'
        % -----------------------------------------------------------------
        % Minima and maxima detection on NCS thorax and/or abdomen data
        % -----------------------------------------------------------------
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
                del2(pkMin2(i+1):pkMin2(i+2)) = abs(ncsData(pkMax2(i),1)-ncsData(pkMin2(i),1)); 
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

        % -----------------------------------------------------------------
        % Perform calibration now.
        % -----------------------------------------------------------------
        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                calibFit = fit([del1(:),del2(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.4,0.6],'Lower',[0,0]);
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                calibFit = fit([del1(:),del2(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.4,0.6,0.01],'Lower',[0,0,-Inf]);
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                calibFit = fit([del1(:),del2(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Lower',[-Inf,-Inf,0,0,-Inf]);
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                calibFit = fit(del1(:),tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Lower',[0,-Inf]);
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                calibFit = fit(del2(:),tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Lower',[0,-Inf]);
            case 'AbdQuad'
                fitFunc = @(A,B,C,x) A.* (x.^2) + B.*x + C;
                calibFit = fit(del2(:),tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.1],'Lower',[0,0,-Inf]);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        calibCoeff = coeffvalues(calibFit);

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del1,del2);

            case 'BiasedLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    del1,del2);
            case 'Quad'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),del1,del2);
            case 'ThLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del1);
            case 'AbdLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),del2);
            case 'AbdQuad'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),del2);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        fig2 = figure;

        if size(ncsData,2) == 2
            figure(fig2);
            ax1(1) = subplot(2,1,1);
            plot(t,volAirflow);
            hold on
            plot(t,ncsData(:,1),t,ncsData(:,2));
            leg = {'Volume','NCS Th','NCS Abd'};
            plotCute1('Time (s)','Volume (L)',ax1(1),[],leg,1);

            ax1(2) = subplot(2,1,2);
            plot(t,tvAirflow); hold on;
            plot(t,vNcs)
            leg = {'TV Airflow','TV NCS'};
            plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        else
            figure(fig2);
            ax1(1) = subplot(2,1,1);
            plot(t,volAirflow);
            hold on
            plot(t,ncsData(:,1));
            leg = {'Airflow Volume','NCS'};
            plotCute1('Time (s)','Volume (L)',ax1(1),[],leg,1);

            ax1(2) = subplot(2,1,2);
            plot(t,tvAirflow); hold on;
            plot(t,vNcs)
            leg = {'TV Airflow','TV NCS'};
            plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        end
        
    case 'vol'
        % -----------------------------------------------------------------
        % Perform calibration now.
        % -----------------------------------------------------------------
        nStart = opts.tCalib(1)*fs+1;
        nEnd = opts.tCalib(2)*fs+1;
        ncs1 = ncsData(:,1); 
        if size(ncsData,2) == 2
            ncs2 = ncsData(:,2);
        end
        
        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
%                 calibFit = fit([ncs1(nStart:nEnd),ncs2(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.4,0.6],'Lower',[0,0],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
                calibFit = fit([ncs1(nStart:nEnd),ncs2(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.4,0.6],'Algorithm','Levenberg-Marquardt' );
            
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                calibFit = fit([ncs1(nStart:nEnd),ncs2(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.4,0.6,0.01],'Lower',[0,0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
%                 calibFit = fit([ncs1(nStart:nEnd),ncs2(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.4,0.6,0.01],'Algorithm','Levenberg-Marquardt');
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
%                 calibFit = fit([ncs1(nStart:nEnd),ncs1(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Lower',[-Inf,-Inf,0,0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
                calibFit = fit([ncs1(nStart:nEnd),ncs1(nStart:nEnd)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Algorithm','Levenberg-Marquardt');
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
%                 calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.9,0.01],'Lower',[0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
                 calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Algorithm','Levenberg-Marquardt');
           case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Lower',[0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
%                 calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.9,0.01],'Algorithm','Levenberg-Marquardt' );
            case 'AbdQuad'
                fitFunc = @(A,B,C,x) A.* (x.^2) + B.*x + C;
%                 calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.2,0.4,0.1],'Lower',[0,0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',1000 );
                calibFit = fit(ncs1(nStart:nEnd),volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.1],'Lower',[0,0,-Inf],'Algorithm','Levenberg-Marquardt');
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        calibCoeff = coeffvalues(calibFit);

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1,ncs2);

            case 'BiasedLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    ncs1,ncs2);
            case 'Quad'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),ncs1,ncs2);
            case 'ThLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1);
            case 'AbdLinear'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),ncs1);
            case 'AbdQuad'
                vNcs = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),ncs1);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        fig2 = figure;

        if size(ncsData,2) == 2
            figure(fig2);
            ax1(1) = subplot(2,1,1);
            plot(t,ncsData(:,1),t,ncsData(:,2));
            leg = {'NCS Th','NCS Abd'};
            plotCute1('Time (s)','mV',ax1(1),[],leg,1);

        else
            figure(fig2);
            ax1(1) = subplot(2,1,1);
            plot(t,ncsData(:,1));
            leg = {'NCS'};
            plotCute1('Time (s)','mV',ax1(1),[],leg,1);
        end
        ax1(2) = subplot(2,1,2);
        plot(t,volAirflow); hold on;
        plot(t,vNcs)
        leg = {'Airflow','NCS'};
        plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        linkaxes(ax1,'x');
        

end