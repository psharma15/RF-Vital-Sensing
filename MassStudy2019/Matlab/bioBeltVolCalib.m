% This function uses airflow TV to calibrate chest belts
% May 07, 2019
% Pragya Sharma, ps847@cornell.edu

function [calibCoeff,vBelt,gof] = bioBeltVolCalib(bioData,fs,tvAirflow,volAirflow,opts)
% Input:
% bioData: Biopac data in column format [Th belt, Abd belt] @fs
% tvAirflow: Tidal/ Lung volume extracted from airflow @fs
% fs: Sampling frequency of the above data in Hz
% calibType: Select what is the calibration equation (linear/ quadratic
% etc.)
% opts: Additional options related to calibration, peak detection etc.
% Output:
% calibCoeff: Calibration coefficient
% vBelt: Calibrated belt tidal volume/volume
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

t = (0:(length(bioData)-1))/fs;

switch (opts.calibType)
    
    case 'tv'
        % Extracting tidal volume for each breath cycle and using that to
        % calibrate. So there are multiple steps
        % -----------------------------------------------------------------
        % STEP 1: Finding peaks in data
        % -----------------------------------------------------------------
        % Minima and maxima detection on BIOPAC thorax and abdomen data
        beltPk = findMaxMin(bioData(:,1:2),fs,opts);
        % -----------------------------------------------------------------
        % STEP 2: Get Uncalibrated TV estimate after peak detection
        % -----------------------------------------------------------------
        % Starting and ending data with expiration end. Not assuming all
        % peaks are available/ detected for thorax and abdomen waveforms.
        if beltPk(1).ind(1) == 1
            % If first peak is inhalation peak, skip it
            beltPk(1).ind = beltPk(1).ind(2:end); 
            beltPk(1).idx = beltPk(1).idx(2:end);
        end
        if beltPk(2).ind(1) == 1
            beltPk(2).ind = beltPk(2).ind(2:end);
            beltPk(2).idx = beltPk(2).idx(2:end);
        end

        if beltPk(1).ind(end) == 1
            beltPk(1).ind = beltPk(1).ind(1:end-1);
            beltPk(1).idx = beltPk(1).idx(1:end-1);
        end
        if beltPk(2).ind(end) == 1
            beltPk(2).ind = beltPk(2).ind(1:end-1);
            beltPk(2).idx = beltPk(2).idx(1:end-1);
        end


        thPkMax = beltPk(1).idx(beltPk(1).ind == 1);
        thPkMin = beltPk(1).idx(beltPk(1).ind == 0);
        abdPkMax = beltPk(2).idx(beltPk(2).ind == 1);
        abdPkMin = beltPk(2).idx(beltPk(2).ind == 0);

        % length(thPkMax)
        % length(thPkMin)
        % length(abdPkMax)
        % length(abdPkMin) 

        figure
        nFig = 2;
        ax(1) = subplot(nFig,1,1);
        plot(t,bioData(:,1)); hold on;
        plot(t(thPkMax),bioData(thPkMax,1),'^',...
            t(thPkMin),bioData(thPkMin,1),'v');
        leg = {'Thorax belt','Max','Min'};
        plotCute1('Time (s)','Belt (mV)',ax(1),'Final Peaks',leg,1,'Horizontal');
        ax(2) = subplot(nFig,1,2);
        plot(t,bioData(:,2)); hold on;
        plot(t(abdPkMax),bioData(abdPkMax,2),'^',...
            t(abdPkMin),bioData(abdPkMin,2),'v');
        leg = {'Abdomen belt','Max','Min'};
        plotCute1('Time (s)','Belt (mV)',ax(2),[],leg,1,'Horizontal');

        % This is the change in Pk-Pk thorax and abdomen signal
        delTh = zeros(length(bioData),1);
        delAbd = zeros(length(bioData),1);

        % So for a cycle, considering there exists 2 minima and 2 maxima point:
        % Calculation is peformed using the difference between maxima and first
        % minima. The update is performed at the end of cycle to be consistent with
        % TV from airflow calculation. And this value is held until next update or
        % the end of the waveform.
        for i = 1:length(thPkMax)
            if i<length(thPkMax)
                delTh(thPkMin(i+1):thPkMin(i+2)) = abs(bioData(thPkMax(i),1)-bioData(thPkMin(i),1)); % Making it positive always
            else
                delTh(thPkMin(i+1):end) = abs(bioData(thPkMax(i),1)-bioData(thPkMin(i),1));
            end
        end
        for i = 1:length(abdPkMax)
            if i<length(abdPkMax)
                delAbd(abdPkMin(i+1):abdPkMin(i+2)) = abs(bioData(abdPkMax(i),2)-bioData(abdPkMin(i),2)); % Making it positive always
            else
                delAbd(abdPkMin(i+1):end) = abs(bioData(abdPkMax(i),2)-bioData(abdPkMin(i),2));
            end
        end

        % ----------------------------------------------------------------
        % STEP 3: Performing calibration to get calib coefficients and TV
        % ----------------------------------------------------------------

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                calibFit = fit([delTh(:),delAbd(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.4,0.6],'Lower',[0,0],...
                                  'Robust', 'LAR' );
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                calibFit = fit([delTh(:),delAbd(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.1,0.1,0.01],'Lower',[0,0,-Inf],...
                                  'Robust', 'LAR' );
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                calibFit = fit([delTh(:),delAbd(:)],tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Lower',[-Inf,-Inf,0,0,-Inf],...
                                  'Robust', 'LAR' );
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                calibFit = fit(delTh(:),tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Lower',[0,-Inf],...
                                  'Robust', 'LAR' );
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                calibFit = fit(delAbd(:),tvAirflow(:),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Lower',[0,-Inf],...
                                  'Robust', 'LAR' );
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        calibCoeff = coeffvalues(calibFit);

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delTh,delAbd);

            case 'BiasedLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    delTh,delAbd);
            case 'Quad'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),delTh,delAbd);
            case 'ThLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delTh);
            case 'AbdLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delAbd);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        figure
        ax1(1) = subplot(2,1,1);
        plot(t,volAirflow);
        hold on
        plot(t,bioData(:,1),t,bioData(:,2));
        leg = {'Volume','Bio Th','Bio Abd'};
        plotCute1('Time (s)','Volume (L)',ax1(1),[],leg,1);

        ax1(2) = subplot(2,1,2);
        plot(t,tvAirflow); hold on;
        plot(t,vBelt)
        leg = {'TV Airflow','TV Belt'};
        plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        
    case 'vol'
        % Calibrate directly the respiration waveforms to get instantaneous
        % lung volume due to inhalation and exhalation. Does not include
        % the default lung volume information.
 
        % ----------------------------------------------------------------
        % Perform calibration now to get calibration coefficients and Vol
        % ----------------------------------------------------------------
        nStart = opts.tCalib(1)*fs+1;
        nEnd = opts.tCalib(2)*fs+1;
        
        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
%                 [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.5,0.5],'Lower',[0,0],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',2000 );
                 [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.5,0.5],'Algorithm','Levenberg-Marquardt' );
           case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.1,0.1,0.01],'Lower',[0,0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',2000 );
%                 [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.1,0.1,0.01],'Algorithm','Levenberg-Marquardt');

            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
%                 [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Lower',[-Inf,-Inf,0,0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',2000 );
                [calibFit, gof] = fit([bioData(nStart:nEnd,1),bioData(nStart:nEnd,2)],volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.2,0.4,0.2,0.4,0.01],'Algorithm','Levenberg-Marquardt');
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
%                 [calibFit, gof] = fit(bioData(nStart:nEnd,1),volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.9,0.01],'Lower',[0,-Inf],'Robust', 'LAR','TolFun',10e-10,'TolX',10e-10,'MaxIter',2000 );
                [calibFit, gof] = fit(bioData(nStart:nEnd,1),volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.9,0.01],'Algorithm','Levenberg-Marquardt');
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
%                 [calibFit, gof] = fit(bioData(nStart:nEnd,2),volAirflow(nStart:nEnd),fitFunc,...
%                                   'StartPoint',[0.001,0.01],'Lower',[0,-Inf],'Robust', 'LAR','TolFun',10e-8,'TolX',10e-8,'MaxIter',2000);
                [calibFit, gof] = fit(bioData(nStart:nEnd,2),volAirflow(nStart:nEnd),fitFunc,...
                                  'StartPoint',[0.001,0.01],'Algorithm','Levenberg-Marquardt');
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        calibCoeff = coeffvalues(calibFit);

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,1),...
                    bioData(:,2));
            case 'BiasedLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    bioData(:,1),bioData(:,2));
            case 'Quad'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),bioData(:,1),bioData(:,2));
            case 'ThLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,1));
            case 'AbdLinear'
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,2));
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

%         fprintf(gof);
        
        figure
        ax1(1) = subplot(2,1,1);
        plot(t,bioData(:,1),t,bioData(:,2));
        leg = {'Bio Th','Bio Abd'};
        plotCute1('Time (s)','Volume (L)',ax1(1),[],leg,1);

        ax1(2) = subplot(2,1,2);
        plot(t,volAirflow); hold on;
        plot(t,vBelt)
        leg = {'Airflow','Belt'};
        plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);

        
    otherwise 
        fprintf('bioBeltVolCalib: Specify valid calibration type.\n');
end

end