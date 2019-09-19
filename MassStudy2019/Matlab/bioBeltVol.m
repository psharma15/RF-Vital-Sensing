% This function estimates Tidal volume from Biopac belts using calibration
% coefficients.
% May 08, 2019
% Pragya Sharma, ps847@cornell.edu

function vBelt = bioBeltVol(bioData,calibCoeff,fs,opts)
% Input:
% bioData: Biopac data in column format [Th belt, Abd belt] @fs
% calibCoeff: Calibration coefficient extracted earlier
% fs: Sampling frequency of the above data in Hz
% calibType: Select what is the calibration equation (linear/ quadratic
% etc.)
% opts: Additional options related to peak detection etc.
% Output:
% tvBelt: Calibrated belt tidal volume
% -------------------------------------------------------------------------

t = (0:(length(bioData)-1))/fs;

switch (opts.calibType)
    
    case 'tv'

        % -------------------------------------------------------------------------
        % Minima and maxima detection on BIOPAC thorax and abdomen data
        % -------------------------------------------------------------------------
        beltPk = findMaxMin(bioData(:,1:2),fs,opts);

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
        plotCute1('Time (s)','Belt (mV)',ax(1),[],leg,1,'Horizontal');
        ax(2) = subplot(nFig,1,2);
        plot(t,bioData(:,2)); hold on;
        plot(t(abdPkMax),bioData(abdPkMax,2),'^',...
            t(abdPkMin),bioData(abdPkMin,2),'v');
        leg = {'Abdomen belt','Max','Min'};
        plotCute1('Time (s)','Belt (mV)',ax(2),[],leg,1,'Horizontal');
        linkaxes(ax,'x')
        % -------------------------------------------------------------------------
        % Get Uncalibrated TV estimate after peak detection
        % -------------------------------------------------------------------------

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

        % -------------------------------------------------------------------------
        % Perform calibration now.
        % -------------------------------------------------------------------------

        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delTh,delAbd);
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    delTh,delAbd);
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),delTh,delAbd);
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delTh);
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),delAbd);
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end
        
    case 'vol'
        % Calibrate directly the respiration waveforms to get instantaneous
        % lung volume due to inhalation and exhalation. Does not include
        % the default lung volume information.
        switch(opts.fitEqn)
            case 'UnbiasedLinear'
                fitFunc = @(A,B,x,y) A.*x + B.*y;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,1),bioData(:,2));
            case 'BiasedLinear'
                fitFunc = @(A,B,C,x,y) A.*x + B.*y + C;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    bioData(:,1),bioData(:,2));
            case 'Quad'
                fitFunc = @(A,B,C,D,E,x,y) A.*(x.^2) + B.*(y.^2) + C.*x + D.*y + E;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),calibCoeff(3),...
                    calibCoeff(4),calibCoeff(5),bioData(:,1),bioData(:,2));
            case 'ThLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,1));
            case 'AbdLinear'
                fitFunc = @(A,B,x)A.*x + B;
                vBelt = fitFunc(calibCoeff(1),calibCoeff(2),bioData(:,2));
            otherwise
                fprintf('Enter correct fitting function. \n'); 
        end

        figure
        ax1(1) = subplot(2,1,1);
        plot(t,bioData(:,1),t,bioData(:,2));
        leg = {'Bio Th','Bio Abd'};
        plotCute1('Time (s)','Volume (L)',ax1(1),[],leg,1);

        ax1(2) = subplot(2,1,2);
        plot(t,vBelt)
        leg = {'Calibrated Belt'};
        plotCute1('Time (s)','Volume (L)',ax1(2),[],leg,1);
        
    otherwise
        fprintf('bioBeltVol: Specify valid calibration type.\n');
end

end