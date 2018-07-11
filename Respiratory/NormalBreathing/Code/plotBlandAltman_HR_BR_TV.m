% This code has different sections - saving the HR, TV, BR data for
% different experiments (only needed once), loading the data, dividing the
% segmenst in normal, deep and shallow breaths, combining the information,
% and plotting Bland Altman plots for these.

%% ------------------------------------------------------------------------
% ONLY RUN THIS ONCE FOR EACH DATA
% Run each experiment and save the data.
% HR, BR, TV are all calculated at same sampling rate, with synchronized
% time.
% saving [time, Hx_HR, NCS_HR_Amp, NCS_HR_Ph, Hx_BR, NCS_BR_Amp, NCS_BR_Ph, Hx_TV, NCS_Tv]
% NCS HR and BR saves both, maybe use just amp or phase. For TV, just save
% ncsTVAmpPhSum. FInal dimenstion will be (nRow x 9). nRow = 419
% data15HrBrTv = [tHR(3:end), hxHR(3:end), ncsHR(3:end,:) , hxBR(3:end), ...
%     ncsBR(3:end,:), hxTV(3:end), ncsTVAmpPhSum(3:end)];
% save('data15HrBrTv.mat','data15HrBrTv');

%% ------------------------------------------------------------------------
% Loading the data
dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\',...
    'Hexoskin\Data\3'];
addpath(dataPath);
% Comments indicate what should be used Amp/ Ph. '*' just shows preferred,
% when both are similar.
load('data08HrBrTv.mat'); % HR: Amp*, BR: Ph*
load('data09HrBrTv.mat'); % HR: Amp, BR: Ph*
load('data11HrBrTv.mat'); % HR: Amp, BR: Ph*
load('data12HrBrTv.mat'); % HR: Ph, BR: Ph
load('data13HrBrTv.mat'); % HR: Ph, BR, Ph*
load('data15HrBrTv.mat'); % HR: Ph, BR: Amp
rmpath(dataPath);

% Truncate first 5 seconds and last 2 seconds of data
nStart = 30;
nEnd = 30;
data08HrBrTv = data08HrBrTv(nStart:(end-nEnd),:);
data09HrBrTv = data09HrBrTv(nStart:(end-nEnd),:);
data11HrBrTv = data11HrBrTv(nStart:(end-nEnd),:);
data12HrBrTv = data12HrBrTv(nStart:(end-nEnd),:);
data13HrBrTv = data13HrBrTv(nStart:(end-nEnd),:);
data15HrBrTv = data15HrBrTv(nStart:(end-nEnd),:);

% dataHxHR = [data08HrBrTv(:,2); data09HrBrTv(:,2); data11HrBrTv(:,2); ...
%             data12HrBrTv(:,2); data13HrBrTv(:,2); data15HrBrTv(:,2)];
% dataNcsHR = [data08HrBrTv(:,3); data09HrBrTv(:,3); data11HrBrTv(:,3); ...
%             data12HrBrTv(:,4); data13HrBrTv(:,4); data15HrBrTv(:,4)];  
% dataHxBR = [data08HrBrTv(:,5); data09HrBrTv(:,5); data11HrBrTv(:,5); ...
%             data12HrBrTv(:,5); data13HrBrTv(:,5); data15HrBrTv(:,5)];
% dataNcsBR = [data08HrBrTv(:,7); data09HrBrTv(:,7); data11HrBrTv(:,7); ...
%              data12HrBrTv(:,7); data13HrBrTv(:,7); data15HrBrTv(:,6)];
% dataHxTV = [data08HrBrTv(:,8); data09HrBrTv(:,8); data11HrBrTv(:,8); ...
%             data12HrBrTv(:,8); data13HrBrTv(:,8); data15HrBrTv(:,8)];
% dataNcsTV = [data08HrBrTv(:,9); data09HrBrTv(:,9); data11HrBrTv(:,9); ...
%              data12HrBrTv(:,9); data13HrBrTv(:,9); data15HrBrTv(:,9)];
         
%% ------------------------------------------------------------------------
% Plotting them to see correctness of data. Not needed after marking above
% either amp or phase should be used.
% dataToPlot = data15HrBrTv;
% t = dataToPlot(:,1);
% figure('Units', 'pixels', 'Position', [100 100 900 600]);
% subplot(3,1,1)
% plot(t,dataToPlot(:,2),t,dataToPlot(:,3),t,dataToPlot(:,4)); 
% legend('Hx HR','NCS Amp HR','NCS Ph HR')
% 
% subplot(3,1,2)
% plot(t,dataToPlot(:,5),t,dataToPlot(:,6),t,dataToPlot(:,7));
% legend('Hx BR','NCS Amp BR','NCS Ph BR')
% 
% subplot(3,1,3)
% plot(t,dataToPlot(:,8),t,dataToPlot(:,9));
% legend('Hx TV', 'NCS Amp Ph Sum TV')

%%
dataHxHR = [data15HrBrTv(:,2); data12HrBrTv(:,2);data13HrBrTv(:,2)]; %data15HrBrTv(:,2); 
dataNcsHR = [data15HrBrTv(:,4); data12HrBrTv(:,4);data13HrBrTv(:,4)]; %data15HrBrTv(:,4);
dataHxBR = [data15HrBrTv(:,5); data12HrBrTv(:,5);data13HrBrTv(:,5)]; %data15HrBrTv(:,5); 
dataNcsBR = [data15HrBrTv(:,6); data12HrBrTv(:,7);data13HrBrTv(:,7)]; %data15HrBrTv(:,6);
dataHxTV = [data15HrBrTv(:,8); data12HrBrTv(:,8);data13HrBrTv(:,8)]; %data15HrBrTv(:,8); 
dataNcsTV = [data15HrBrTv(:,9); data12HrBrTv(:,9);data13HrBrTv(:,9)]; %data15HrBrTv(:,9); 

%% ------------------------------------------------------------------------
% Bland-Altman Plot for HR
baPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\Respiratory\Hexoskin\Code\BlandAltman\BlandAltman';
addpath(baPath);

colors = [0.3,0.75,0.93; ...
          0, 0, 1; ...
          1, 1,0];

labelHR = {'Hexoskin HR','NCS HR','BPM'};
[hrBA, hrFig, hrStats] = BlandAltman(dataHxHR, dataNcsHR, labelHR, [],[],'colors',colors,'on');

labelBR = {'Hexoskin BR','NCS BR','BPM'};
[brBA, brFig, brStats] = BlandAltman(dataHxBR, dataNcsBR, labelBR, [],[],'colors',colors,'on');

labelTV = {'Hexoskin TV','NCS TV','ml'};
[brTV, tvFig, tvStats] = BlandAltman(dataHxTV, dataNcsTV, labelTV, [],[],'colors',colors,'on');

rmpath(baPath);

% Standard deviations are stats.differenceSTD

%% 
% The RMSE is calculated after linear fitting, which is not desired. Hence,
% calculating separately.
rmseTV = sqrt(sum((dataHxTV-dataNcsTV).^2)/length(dataNcsTV));
rmseHR = sqrt(sum((dataHxHR-dataNcsHR).^2)/length(dataNcsHR));
rmseBR = sqrt(sum((dataHxBR-dataNcsBR).^2)/length(dataNcsBR));
