%% MOTIONDETECT1.m: NCS motion detection code by ABNORMAL signal detection.
% This code only works on amplitude or phase. No windowing.
% Parts of code:
%              - Beat Extraction
%              - Feature Extraction
%              - SVM Training - Semi-supervised learning
%              - SVM abnormal signal classification

clearvars; 
rng(0); % SVM classifier 'KernelScale' = 'auto' selects scale heuristically, needs repeatability.
fs = 500;
threshPrediction = 0; % 0 is the standard.
f3db = 0.6;
fp = 10; fst = 11;

%% PREPROCESSING
dataPath = ['D:\Research\SummerFall17Spring18\CnC\NCS\MotionDetectionSleep\',...
    'MyExperiment\DataCollected\Dec02'];
codePostProcessPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\MotionDetectionSleep\MyExperiment';
addpath(codePostProcessPath);

% This section combines the entire data set. Specify start (and stop) time
% if required. Skip the initial 2 deep breaths from the data set
% (approximately 30-60 seconds). Check that the data is arranged in
% appropriate format.
[ncs25,~,ncsUnfiltAmp25,ncsUnfiltPh25,~] = postProcess2(f3db,fp,fst,dataPath,'b8_5.mat',-1);
[ncs24,~,ncsUnfiltAmp24,ncsUnfiltPh24,~] = postProcess2(f3db,fp,fst,dataPath,'b8_4.mat',-1,420,535);
[ncs23,~,ncsUnfiltAmp23,ncsUnfiltPh23,~] = postProcess2(f3db,fp,fst,dataPath,'b8_4.mat',-1,140,290);
[ncs22,~,ncsUnfiltAmp22,ncsUnfiltPh22,~] = postProcess2(f3db,fp,fst,dataPath,'b8_3.mat',-1,0,265);
[ncs21,~,ncsUnfiltAmp21,ncsUnfiltPh21,~] = postProcess2(f3db,fp,fst,dataPath,'b8_1.mat',-1,0);
[ncs20,~,ncsUnfiltAmp20,ncsUnfiltPh20,~] = postProcess2(f3db,fp,fst,dataPath,'b6_6.mat',-1,0);
[ncs19,~,ncsUnfiltAmp19,ncsUnfiltPh19,~] = postProcess2(f3db,fp,fst,dataPath,'b6_5.mat',-1,0); % More like slight torso motion
[ncs18,~,ncsUnfiltAmp18,ncsUnfiltPh18,~] = postProcess2(f3db,fp,fst,dataPath,'b6_4.mat',-1,0);
[ncs17,~,ncsUnfiltAmp17,ncsUnfiltPh17,~] = postProcess2(f3db,fp,fst,dataPath,'b6_2.mat',-1,0,320);
[ncs16,~,ncsUnfiltAmp16,ncsUnfiltPh16,~] = postProcess2(f3db,fp,fst,dataPath,'b6_1.mat',-1,0);
[ncs15,~,ncsUnfiltAmp15,ncsUnfiltPh15,~] = postProcess2(f3db,fp,fst,dataPath,'b5_5.mat',1,0);
[ncs14,~,ncsUnfiltAmp14,ncsUnfiltPh14,~] = postProcess2(f3db,fp,fst,dataPath,'b5_3.mat',-1,355);
[ncs13,~,ncsUnfiltAmp13,ncsUnfiltPh13,~] = postProcess2(f3db,fp,fst,dataPath,'b5_2.mat',-1,0);
[ncs12,~,ncsUnfiltAmp12,ncsUnfiltPh12,~] = postProcess2(f3db,fp,fst,dataPath,'b5_1.mat',-1,0);
[ncs11,~,ncsUnfiltAmp11,ncsUnfiltPh11,~] = postProcess2(f3db,fp,fst,dataPath,'b1_11.mat',-1,0); % Long
[ncs10,~,ncsUnfiltAmp10,ncsUnfiltPh10,~] = postProcess2(f3db,fp,fst,dataPath,'b1_10.mat',-1,0);
[ncs9,~,ncsUnfiltAmp9,ncsUnfiltPh9,~] = postProcess2(f3db,fp,fst,dataPath,'b1_9.mat',-1,0);
[ncs8,~,ncsUnfiltAmp8,ncsUnfiltPh8,~] = postProcess2(f3db,fp,fst,dataPath,'b1_8.mat',-1,0);
[ncs7,~,ncsUnfiltAmp7,ncsUnfiltPh7,~] = postProcess2(f3db,fp,fst,dataPath,'b1_7.mat',-1,0);
[ncs6,~,ncsUnfiltAmp6,ncsUnfiltPh6,~] = postProcess2(f3db,fp,fst,dataPath,'b1_6.mat',-1,46);
[ncs5,~,ncsUnfiltAmp5,ncsUnfiltPh5,~] = postProcess2(f3db,fp,fst,dataPath,'b1_5.mat',-1,50);
[ncs2,~,ncsUnfiltAmp2,ncsUnfiltPh2,~] = postProcess2(f3db,fp,fst,dataPath,'b1_2.mat',-1,125);

close all;

% Make sure that the data is arranged in such a way that first two minutes
% has no motion at all. It should be normal spontaneous breathing.
ncs = [ncs2;ncs5;ncs6;ncs7;ncs8;ncs9;ncs10;ncs11;ncs12;ncs13;ncs14;ncs15;ncs16;...
    ncs17;ncs18;ncs19;ncs20;ncs21;ncs22;ncs23;ncs24;ncs25];

ncsAmpUnfilt = [ncsUnfiltAmp2;ncsUnfiltAmp5;ncsUnfiltAmp6;ncsUnfiltAmp7;...
    ncsUnfiltAmp8;ncsUnfiltAmp9;ncsUnfiltAmp10;ncsUnfiltAmp11;ncsUnfiltAmp12;...
    ncsUnfiltAmp13;ncsUnfiltAmp14;ncsUnfiltAmp15;ncsUnfiltAmp16;ncsUnfiltAmp17;...
    ncsUnfiltAmp18;ncsUnfiltAmp19;ncsUnfiltAmp20;ncsUnfiltAmp21;ncsUnfiltAmp22;...
    ncsUnfiltAmp23;ncsUnfiltAmp24;ncsUnfiltAmp25];

ncsPhUnfilt = [ncsUnfiltPh2;ncsUnfiltPh5;ncsUnfiltAmp6;ncsUnfiltPh7;...
    ncsUnfiltPh8;ncsUnfiltPh9;ncsUnfiltPh10;ncsUnfiltPh11;ncsUnfiltPh12;...
    ncsUnfiltPh13;ncsUnfiltPh14;ncsUnfiltPh15;ncsUnfiltPh16;ncsUnfiltPh17;...
    ncsUnfiltPh18;ncsUnfiltPh19;ncsUnfiltPh20;ncsUnfiltPh21;ncsUnfiltPh22;...
    ncsUnfiltPh23;ncsUnfiltPh24;ncsUnfiltPh25];

clearvars -except fs ncs ncsAmpUnfilt ncsPhUnfilt codePostProcessPath 

nSample = length(ncs);
idx = 1:nSample;
t = ((idx-1)/fs)';

figure
yyaxis left
plot(t,ncsAmpUnfilt)
yyaxis right
plot(t,ncsPhUnfilt)

rmpath(codePostProcessPath);

%%
codeMDPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\MotionDetectionSleep\MDCodes';
addpath(codeMDPath);

ncsSample = ncs;
tSample = t;

% Wavelet decomposition and D-8 reconstruction
wName = 'db10';
[recNcsd8,~,~] = wavedecrec1(ncsSample,tSample,wName,8,1);
[recNcsd4,~,~] = wavedecrec1(ncsSample,tSample,wName,4,0);

figure
ax1 = gca;
plot(ax1,tSample,recNcsd8);
xlabel(ax1,'Time(sec)')
ylabel(ax1,'Amplitude')
title(ax1,'NCS - d8')
grid on; 

%% 
nStart = 1;
windowLen = 1000*60*fs; % Entire data set. No windowing in this code.
minPeakInterval = round(fs*60/200); % Apriori Condition: max heart rate 180 or 200 beats/min

nFeature = 7;
% nBeatTrain = 70; % Number of beats for training from beginning
tTrain = 120; % Time for training from beginning
beatWindowDetect = 2; % Symmetric 2 beats before and after current beat.

X = [];

figure('Units', 'pixels', ...
    'Position', [100 100 500 300]);
ax2 = gca;
h1 = plot(ax2,t,ncs,':');
hold on

h2 = [];
beatTimeLast = 60/90; % Starting from assuming 90 beats/min;

while (nStart <= length(tSample)) && (tSample(nStart) <= tSample(end))
    nEnd = nStart + windowLen - 1;
    if nEnd > length(tSample)
        nEnd = length(tSample);
    end
    if nEnd-nStart <= 1*fs
        fprintf('Window length too small for analysis, nStart: %d.\n',nStart);
        break
    end
    
    ncsWindow = ncs(nStart:nEnd);
    recNcsd8Window = recNcsd8(nStart:nEnd);
    recNcsd4Window = recNcsd4(nStart:nEnd);
    tWindow = tSample(nStart:nEnd);
    
    %% Peak Detection
    % Simple First-order Derivative, with just one apriori condition
    [~,locs] = findpeaks(recNcsd8Window,'MinPeakDistance',minPeakInterval);
    
    nStart = nEnd + 1;
    
    nBeatTrain = find(t(locs)>tTrain,1) - 1;

    if length(locs) < 2
        fprintf('No beat detected in this time frame. \n');
        if nEnd ~= length(tSample)
            continue
        else 
            break
        end
    end
    
    nBeatsWindow = length(locs) - 1;
    X = zeros(nBeatsWindow,nFeature);
    relIdx = zeros(nBeatsWindow,2); % Start and End Idx relative to ncs original
    
    %% BEAT LEVEL PROCESSING
    for i = 1:length(locs)-1     
        if i>3 && i<=length(locs)-3
            ncsWindowPSD = ncsWindow(round((locs(i-3)+locs(i-2))/2):round((locs(i+2)+locs(i+3))/2)); % PSD over a window
        elseif i>length(locs)-3
            ncsWindowPSD = ncsWindow(locs(i-1):locs(i+1));
        elseif i<=3
            ncsWindowPSD = ncsWindow(locs(i):locs(i+2));
        end
%         ncsWindowStats = ncsWindowPSD;
        if i>1 && i<=length(locs)-2
            ncsWindowStats = ncsWindow(round((locs(i-1)+locs(i))/2):round((locs(i+1)+locs(i+2))/2)); % HOS over a window
        elseif i> length(locs) - 2
            ncsWindowStats = ncsWindow(round((locs(i-1)+locs(i))/2):locs(i+1));
        elseif i<=1
            ncsWindowStats = ncsWindow(round((locs(i+1)+locs(i+2))/2):locs(i+2));
        end
        
        recNcsd4Beat = recNcsd4Window(locs(i):locs(i+1));
        tBeat = tWindow(locs(i):locs(i+1));
        
        %% Feature detection per beat

        % Feature a2: Freq Domain: Should be on NCS, not wavelet reconstructed
        psdest = psd(spectrum.periodogram,ncsWindowPSD,'Fs',fs,'NFFT',length(ncsWindowPSD));
        psdPower = avgpower(psdest,[0.8,10]);
        normPsdPower = psdPower/avgpower(psdest);

        % Feature Time Domain:
        rmsBeat = rms(recNcsd4Beat);   
        
        if i == length(locs)-1
            relRmsBeat = rms(ncsWindow(locs(i):locs(i+1)))/rms(ncsWindow(locs(i-1):locs(i)));
        else
            relRmsBeat = rms(ncsWindow(locs(i):locs(i+1)))/rms(ncsWindow(locs(i+1):locs(i+2)));
        end
        
        beatTime = tWindow(locs(i+1))-tWindow(locs(i));
        relIBI = beatTime/beatTimeLast;
        
        meanBeat = mean(ncsWindowStats);
        varBeat = std(ncsWindowStats);
        skewnessBeat = skewness(ncsWindowStats); % Skewness
        kurtosisBeat = kurtosis(ncsWindowStats); % Kurtosis
        entropyBeat = wentropy(ncsWindowStats,'shannon'); % Shannon's entropy

        beatTimeLast = beatTime;
        
        %% IV.a: Feature collection
        % Possible features:
        % rmsBeat,beatTime,psdPower,normPsdPower,skewnessBeat,kurtosisBeat,
        % entropyBeat, meanBeat, relIBI, varBeat, relRmsPower
        if (nFeature == 3)
            X(i,:) = [normPsdPower,beatTime,entropyBeat];
        elseif (nFeature == 4)
            X(i,:) = [normPsdPower,relRmsBeat,meanBeat,kurtosisBeat];
        elseif nFeature == 5
            X(i,:) = [varBeat,meanBeat,normPsdPower,relRmsBeat,kurtosisBeat];
        elseif (nFeature == 6)
            X(i,:) = [varBeat,meanBeat,normPsdPower,relRmsBeat,kurtosisBeat,...
                skewnessBeat];
        elseif nFeature == 7
            X(i,:) = [relIBI,normPsdPower,relRmsBeat,meanBeat,varBeat,skewnessBeat,...
                kurtosisBeat];
        elseif nFeature == 8
            X(i,:) = [relIBI,normPsdPower,relRmsBeat,meanBeat,varBeat,skewnessBeat,...
                kurtosisBeat,entropyBeat];
        else
            fprintf('Wrong number of features selected.\n')
            return
        end
        
        relIdx(i,:) = [find(tSample == tBeat(1)), find(tSample == tBeat(end))];
        
    end

    %%  Abnormality detection
    X(1,1) = X(2,1); % If first feature is IBI - it is duplicated for first instance;
    y = ones(nBeatsWindow,1);
    rng(1);
    XNorm = (X - mean(X(1:nBeatTrain,:)))./std(X(1:nBeatTrain,:)); %Normalizing only with training data

    [coef,pcaScore,~,~,explained] = pca(XNorm);

    fprintf('%% of variance explained by each principal component:\n');
    fprintf('%f \n',explained');

    figure
    if (nFeature == 3)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3'});
    elseif (nFeature == 4)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4'});
    elseif (nFeature == 5)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4','v_5'});
    elseif (nFeature == 6)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4','v_5','v_6'});
    elseif (nFeature == 7)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4','v_5','v_6','v_7'});
    elseif (nFeature == 8)
        biplot(coef(:,1:2),'scores',pcaScore(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4','v_5','v_6','v_7','v_8'});
    else
        fprintf('Figure: Wrong number of features selected.\n')
    end

    figure
    scatter3(pcaScore(:,1),pcaScore(:,2),pcaScore(:,3))
    axis equal
    xlabel('1st PC')
    ylabel('2nd PC')
    zlabel('3rd OC')
        
%     % Using PCA
%     model = fitcsvm(pcaScore(1:nBeatTrain,1:3),y(1:nBeatTrain),'KernelScale','auto');
%     [~,scorePred] = predict(model,pcaScore(nBeatTrain+1:end,1:3));
        
    % First nBeatTrain number of beats for training
    model = fitcsvm(XNorm(1:nBeatTrain,:),y(1:nBeatTrain),'KernelScale','auto',...
        'Nu',0.5,'Prior',0.98);

%     model = fitcsvm(XNorm(1:nBeatTrain,:),y(1:nBeatTrain),'KernelScale','auto',...
%         'scoreTransform','symmetriclogit');

    [~,scorePred] = predict(model,XNorm(nBeatTrain+1:end,:));

    %% 
    motionStartTime = [tSample(relIdx(nBeatTrain+1,1));1913.5;1920;1974;1980    ;...
        2016;2022;2069.5;2071;2075;2083;2136;2142;2196;2201;2257;2262;2279      ;...
        2285;2318;2323;2378;2387;2394;2394.5;2438;2448;2480;2481;2498;2505      ;...
        2557;2563;2618;2623;2678;2683;2737;2744;2798;2804;2858;2863;2865        ;...
        2866;2918;2923;2978;2984;3038;3044;3098;3104;3157;3160;3217;3222        ;...
        3277;3282;3337;3342;3398;3401;3421;3422;3458;3461;3518;3523;3580        ;...
        3583;3640;3644;3700;3702;3759;3764;3819;3823;3839;3841;3899;3902        ;...
        3959;3965;4019.5;4023;4080;4084;4139;4143;4259;4264;4319;4323           ;...
        4379;4383;4438;4442;4500;4505;4559;4564;4620;4623;4637;4641;4697        ;...
        4702;4729;4730;4758;4762;4819;4823;4879;4881;4939;4940;4998;5003        ;...
        5058;5064;5119;5120;5179;5183;5299;5303;5359;5362;5419;5423;5479        ;...
        5485;5539;5542;5598;5604;5718;5723;5778;5782;5899;5903;5958.5;5965      ;...
        6018;6023;6079;6091;6144;6149;6203;6209;6253;6260;6314;6319;6489        ;...
        6494;6550;6552;6609;6611;6670;6673;6729;6733;6789;6792;6848;6853        ;...
        6909;6912;tSample(relIdx(end,2))];
    classLabelBeat = zeros(length(scorePred)+nBeatTrain,1);
    motionStatus = zeros(length(motionStartTime)-1,1);
    beatNum = zeros(length(motionStartTime),1);
    idx = zeros(length(motionStartTime),1);
    motionStatus(1) = 1;
    for i = 1:length(motionStartTime)
        idx(i) = find(tSample>= motionStartTime(i),1);
        beatNum(i) = find(relIdx(:,1) <= idx(i),1,'last');
        if i>1
            classLabelBeat(beatNum(i-1):beatNum(i)-1) = motionStatus(i-1);
            motionStatus(i) = motionStatus(i-1) * -1;
        end
    end
     
%% PREDICTION
    normScorePred = tanh((scorePred));
    % Depends on if you have used 'scoreTransform' in SVM or not. If used,
    % no need to do anything. Otherwise use tanh normalization. Both
    % functions are slightly different. g(x) = 1/(1+e^-x). tanh(x) =
    % 2g(2x)-1. 'symmetriclogit' = 2g(x) - 1. There is slight difference in
    % answer. 
%     normScorePred = scorePred; 
    w1 = 0.65; w2 = 0.4;
    classPredBeat = zeros(length(scorePred),1);
    classPredTimeSeries = zeros(length(tSample),1); % Zero indicates training

    % Prediction for one beat depends on nearby beats
    for i = 1:length(scorePred)
        if i < beatWindowDetect+1
            classPredBeat(i) = [1,w1,w2]*normScorePred(i:i+beatWindowDetect);
        elseif i>length(scorePred)-beatWindowDetect
            classPredBeat(i) = [w2,w1,1]*normScorePred(i-beatWindowDetect:i);
        else
            classPredBeat(i) = [w2,w1,1,w1,w2]*normScorePred(i-beatWindowDetect:i+beatWindowDetect);
        end
        relIdxCount = i+nBeatTrain; 
        if classPredBeat(i) < threshPrediction
            % MOTION
            classPredTimeSeries(relIdx(relIdxCount,1):relIdx(relIdxCount,2)) = -1; % Motion
            h2 = plot(ax2,tSample(relIdx(relIdxCount,1):relIdx(relIdxCount,2)),...
                ncsSample(relIdx(relIdxCount,1):relIdx(relIdxCount,2)),'color',[0.3,0.75,0.93]); %[0.3,0.75,0.93]
        else
            classPredTimeSeries(relIdx(relIdxCount,1):relIdx(relIdxCount,2)) = 1; % Motion
            h3 = plot(ax2,tSample(relIdx(relIdxCount,1):relIdx(relIdxCount,2)),...
                ncsSample(relIdx(relIdxCount,1):relIdx(relIdxCount,2)),'k');
        end        
    end
   classPredBeat = [zeros(nBeatTrain,1);sign(classPredBeat)];
    
  %% Moving average heart rate
   windowSize = 30; % 30 beats to calculate hr
   nTotBeat = size(relIdx,1);
   hr = zeros(nTotBeat,1);
   hrMotionCorrect = zeros(nTotBeat,1);
   tWindow = zeros(nTotBeat,1);
   tWindowMotionCorrect = zeros(nTotBeat,1);
   
   tBuffer = 0;
   for i = nBeatTrain+1+windowSize:nTotBeat
       tWindow(i) = tSample(relIdx(i,2))-tSample(relIdx(i-windowSize+1,1));
       hr(i) = 60*windowSize/tWindow(i);
       
       windowFill = 0; % #Beats in buffer
       counter = 1;
       while (windowFill < windowSize) && (i-windowSize+counter)<=nTotBeat
           if classPredBeat(i-windowSize+counter) == 1
               windowFill = windowFill + 1;
               tWindowMotionCorrect(i) = tWindowMotionCorrect(i) + tSample(relIdx(i-windowSize+counter,2))-tSample(relIdx(i-windowSize+counter,1));
           end
           counter = counter + 1;
       end
       hrMotionCorrect(i) = 60*windowSize/tWindowMotionCorrect(i);
   end

end

hXLabel = xlabel(ax2,'Time (sec)','FontSize',10);
hYLabel = ylabel(ax2,'Heartbeat (a.u.)','FontSize',10); 
classPredTimeSeriesMotion = (classPredTimeSeries == -1);
k = find(classPredTimeSeriesMotion);
% h4 = plot(ax2,tSample(k),0.0035*classPredTimeSeriesMotion(k),'o','MarkerSize',2.5,'MarkerEdgeColor' , 'none'      , ...
%   'MarkerFaceColor' , [0.3,0.75,0.93]);
h4 = plot(ax2,tSample,2.25e-3+0.0015*classPredTimeSeriesMotion,'color' , [0.3,0.75,0.93]);

% tFalseRest = 799.5:1/fs:803;
% h5 = plot(ax2,tFalseRest,-0.0035*ones(length(tFalseRest),1),'o','MarkerSize',2.5,'MarkerEdgeColor' , 'none'      , ...
%   'MarkerFaceColor' , [.75 .75 1]);
hLabel = legend(h4,'Motion');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hLabel, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );

set(ax2, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , -0.008:0.002:0.008, ...
  'LineWidth'   , 1         );
grid(ax2,'on')

figure;
ax3 = gca;
plot(t,ncsAmpUnfilt); grid on;
linkaxes([ax2,ax3],'x')

grid(ax3,'on')

figure
tBeat = tSample(relIdx(:,1));
grid on
plot(tBeat(nBeatTrain+1+windowSize:end),hr(nBeatTrain+1+windowSize:end),'k','LineWidth',2)
hold on;
plot(tBeat(nBeatTrain+1+windowSize:end),hrMotionCorrect(nBeatTrain+1+windowSize:end),'Color',[0.85,0.33,0.1]); 
ylabel('HR(bpm)')
xlabel('Time(sec)')
title('Heart Rate Motion Correction','FontSize',12)
legend('HR: No motion Correction','HR: Motion Correction')
hold off; 

figure
ax4(1) = subplot(2,1,1);
nTotBeat = nBeatTrain+length(scorePred);
plot(1:nTotBeat,classLabelBeat); grid on; 
ylim([-1.25,1.25])
ax4(2) = subplot(2,1,2);
plot(1:nTotBeat,classPredBeat); grid on
ylim([-1.25,1.25])
linkaxes(ax4,'x')


TP = sum((classLabelBeat == -1) & (classPredBeat == classLabelBeat)); % True motion detected
FP = sum((classLabelBeat == 1) & (classPredBeat ~= classLabelBeat)); % False motion detected
TN = sum((classLabelBeat == 1) & (classPredBeat == classLabelBeat)); % True rest state detected
FN = sum((classLabelBeat == -1) & (classPredBeat ~= classLabelBeat)); % Flase rest state detected

acc = 100*(TP+TN)/(TP+FP+FN+TN);
sensitivity = 100*TP/(TP+FN);
specificity = 100*TN/(TN+FP);
f1Score = 2*TP/(2*TP+FP+FN);

fprintf('TP,FP,TN,FN = %d, %d, %d, %d\n',TP,FP,TN,FN);
fprintf('Accuracy = %2.2f%%, \nSensitivity = %2.2f%%, \nSpecificity = %2.2f%% \n',acc,sensitivity,specificity);

rmpath(codeMDPath);
