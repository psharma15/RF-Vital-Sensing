% Result plot: This code plots the final figures for different papers. 
% Results are obtained from ncsHxCompare.m

%%
tStart = 205; tEnd = 285;

% Knowing that sample rate is 1Hz, indicator is defined for normal (1),
% deep (2), fast shallow (3) and rapid shallow (4).
% indicator = [2.* ones(26,1); 1.* ones(129,1); 2.*ones(58,1); 1.*ones(120,1);...
%     3.*ones(30,1); 1.*ones(151,1); 4.* ones(33,1); 1.*ones(19,1); 4.* ones(20,1);1.*ones(15,1)]; %0215_1423
% indicator = [2.* ones(28,1); 1.* ones(132,1); 2.*ones(53,1); 1.*ones(59,1);...
%     3.*ones(30,1); 1.*ones(89,1); 4.* ones(29,1); 1.*ones(181,1)]; %0215_1410

% indicator = 5.*ones(tEnd-tStart+1,1);

% ncs_BR_TV_hx_BR_TV = [indicator(tStart:tEnd),...
%                       ncsBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),...
%                       ncsTVAmp(tStart*hxSampRateTV:tEnd*hxSampRateTV),...
%                       hxBR2(tStart*hxSampRateBR:tEnd*hxSampRateBR,2),...
%                       hxTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)];
% save(['E:\NCS\Respiratory\BreathPattern_Cough_Speak\Data\Pragya\JanFeb2019\Feb15',...
%     '\brtv_ataxic_0215_1552_quadTV.mat'],'ncs_BR_TV_hx_BR_TV');
% save(['E:\NCS\Respiratory\NormalBreathing\Data\v2\Feb15_2019','\brtv_0215_1410_quadCalib.mat'],'ncs_BR_TV_hx_BR_TV');
figure('Units', 'pixels','Position', [100 100 400 320]);
nFig = 2;
ax1(1) = subplot(nFig,1,1);
yyaxis left

ncs1Min = min(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Max = max(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Norm = -1 + (2.*((ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));

plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    ncs1Norm,...
    'color','k','LineWidth',1);
plotCute1('Time (s)','NCS (a.u.)',ax1(1),[],[],0);
ylim([-1.1,1.8])
yyaxis right
plot(ax1(1),tBR(tStart*hxSampRateBR:tEnd*hxSampRateBR)-tStart,...
    ncsBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),':','color',...
    [0.3,0.75,0.93],'LineWidth',2.5);
plotCute1('Time (s)','BR (BPM)',ax1(1),[],{'NCS','BR'},1,'Horizontal');
ylim([0,120])
yticks(0:40:120)

xlim([tStart-tStart,tEnd-tStart])

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),tTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)-tStart,...
    hxTV(tStart*hxSampRateTV:tEnd*hxSampRateTV,1),':','color',...
    [0.3,0.75,0.93],'LineWidth',2.5);
hold on
plot(ax1(2),tTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)-tStart,...
    ncsTVAmp(tStart*hxSampRateTV:tEnd*hxSampRateTV,1),'color',...
    'k','LineWidth',1.5);
plotCute1('Time (s)','LV (mL)',ax1(2),[],{'Hx','NCS'},1,'Horizontal');
ylim([0,800])
yticks(0:200:800)
xlim([tStart-tStart,tEnd-tStart])
linkaxes(ax1(:),'x');



%%
tStart = 600; tEnd = 700;

ncs_BR_TV_hx_BR_TV = [ncsBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),...
                      ncsTVAmp(tStart*hxSampRateTV:tEnd*hxSampRateTV),...
                      hxBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),...
                      hxTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)];
% save(['E:\NCS\Respiratory\NormalBreathing\Data\v2\Feb15_2019','\brtv_0215_1423_2.mat'],'ncs_BR_TV_hx_BR_TV');
% save('deepNormal_1813_A.mat','ncs_BR_TV_hx_BR_TV'); % So that we don't overwrite by mistake

figure('Units', 'pixels','Position', [100 100 400 320]);
nFig = 2;
ax1(1) = subplot(nFig,1,1);

ncs1Min = min(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Max = max(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Norm = -1 + (2.*((ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));

plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    ncs1Norm,...
    'color','k','LineWidth',1);
plotCute1('Time (s)','NCS (a.u.)',ax1(1),[],[],0);
ylim([-1.3,1.3])

xlim([tStart-tStart,tEnd-tStart])

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),tTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)-tStart,...
    hxTV(tStart*hxSampRateTV:tEnd*hxSampRateTV,1),':','color',...
    [0.3,0.75,0.93],'LineWidth',2.5);
hold on
plot(ax1(2),tTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)-tStart,...
    ncsTVAmp(tStart*hxSampRateTV:tEnd*hxSampRateTV,1),'color',...
    'k','LineWidth',1.5);
plotCute1('Time (s)','LV (mL)',ax1(2),[],{'Hx','NCS'},1,'Horizontal');
ylim([0,1600])
yticks(0:400:1600)
xlim([tStart-tStart,tEnd-tStart])
linkaxes(ax1(:),'x');

%%
tStart = 460; tEnd = 520;

% ncs_BR_TV_hx_BR_TV = [ncsBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),...
%                       ncsTVAmp(tStart*hxSampRateTV:tEnd*hxSampRateTV),...
%                       hxBR(tStart*hxSampRateBR:tEnd*hxSampRateBR),...
%                       hxTV(tStart*hxSampRateTV:tEnd*hxSampRateTV)];
% save([dataPath,'\normal_1855_A.mat'],'ncs_BR_TV_hx_BR_TV');

figure('Units', 'pixels','Position', [100 100 400 320]);
nFig = 2;
ax1(1) = subplot(nFig,1,1);

ncs1Min = min(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Max = max(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
hxThMin = min(hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxThMax = max(hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxAbdMin = min(hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxAbdMax = max(hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate));

ncs1Norm = -1 + (2.*((ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));
hxThNorm = -1 + (2.*((hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate) - hxThMin)/(hxThMax-hxThMin)));
hxAbdNorm = -1 + (2.*((hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate) - hxAbdMin)/(hxAbdMax-hxAbdMin)));

plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    hxThNorm,':',...
    'LineWidth',2);
hold on
plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    hxAbdNorm,'--',...
    'LineWidth',1.5);

plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    ncs1Norm,'color','k',...
    'LineWidth',1);
plotCute1('Time (s)','Respiration (a.u.)',ax1(1),[],{'Hx Th','Hx Abd','Ncs'},1,'Horizontal');
xlim([tStart-tStart,tEnd-tStart])
ylim([-1.8,1])

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),tBR(tStart*hxSampRateBR:tEnd*hxSampRateBR)-tStart,...
    hxBR(tStart*hxSampRateBR:tEnd*hxSampRateBR),':','color',...
    [0.3,0.75,0.93],'LineWidth',2.5);
hold on
plot(ax1(2),tBR(tStart*hxSampRateBR:tEnd*hxSampRateBR)-tStart,...
    ncsBR(tStart*hxSampRateBR:tEnd*hxSampRateBR,1),'color',...
    'k','LineWidth',1.5);
plotCute1('Time (s)','BR (BPM)',ax1(2),[],{'Hx','NCS'},1,'Horizontal');
ylim([10,20])

xlim([tStart-tStart,tEnd-tStart])
linkaxes(ax1(:),'x');

%% Cough plot 1
tStart = 600; tEnd = 700;

figure('Units', 'pixels','Position', [100 100 400 320]);
nFig = 2;

ncs1Min = min(ncsRespFiltered(tStart*ncsSampRate+1:tEnd*ncsSampRate,1));
ncs1Max = max(ncsRespFiltered(tStart*ncsSampRate+1:tEnd*ncsSampRate,1));

ncsLPNorm = -1 + (2.*((ncsRespFiltered(tStart*ncsSampRate+1:tEnd*ncsSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));
ncsHPNorm = -1 + (2.*((ncsRespHP(tStart*ncsSampRate+1:tEnd*ncsSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));
ncsHPNorm = ncsHPNorm - mean(ncsHPNorm);

ncsHFenergy = ncsHPNorm.^2;

t = 0:1/ncsSampRate:(length(ncsHFenergy)-1)/ncsSampRate;
tCalib = 0:10;
calibEnergy = mean(ncsHFenergy(tCalib(1)*ncsSampRate+1:tCalib(2)*ncsSampRate));
ind = ncsHFenergy>100*calibEnergy;

ax1(1) = subplot(nFig,1,1);
plot(ax1(1),TFilt(tStart*ncsSampRate+1:tEnd*ncsSampRate)-tStart,...
    ncsLPNorm,...
    'LineWidth',1);
plotCute1('Time (s)','NCS_L_F(a.u.)',ax1(1),[],[],0);
ylim([-1,1])

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),TFilt(tStart*ncsSampRate+1:tEnd*ncsSampRate)-tStart,...
    ncsHPNorm,...
    'LineWidth',1);
hold on
% plot(ax1(2),t,ncsHFenergy);
plot(ax1(2),t,0.02.*ind+0.05);
ylim([-0.06,0.08])

plotCute1('Time (s)','NCS_H_F(a.u.)',ax1(2),[],[],0);


%% This is 2 sensor plot: Thorax and abdomen
tStart = 5; tEnd = 65;

ncs1Min = min(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs1Max = max(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1));
ncs2Min = min(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,2));
ncs2Max = max(ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,2));

hxThMin = min(hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxThMax = max(hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxAbdMin = min(hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate));
hxAbdMax = max(hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate));

ncs1Norm = -1 + (2.*((ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,1) - ncs1Min)/(ncs1Max-ncs1Min)));
ncs2Norm = -1 + (2.*((ncsResp(tStart*hxRespSampRate:tEnd*hxRespSampRate,2) - ncs2Min)/(ncs2Max-ncs2Min)));

hxThNorm = -1 + (2.*((hxRespTh(tStart*hxRespSampRate:tEnd*hxRespSampRate) - hxThMin)/(hxThMax-hxThMin)));
hxAbdNorm = -1 + (2.*((hxRespAbd(tStart*hxRespSampRate:tEnd*hxRespSampRate) - hxAbdMin)/(hxAbdMax-hxAbdMin)));

figure('Units', 'pixels','Position', [100 100 400 350]);
nFig = 2;
ax1(1) = subplot(nFig,1,1);
% yyaxis left
plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    hxThNorm,':',...
    'color','k','LineWidth',2);
hold on
plot(ax1(1),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    hxAbdNorm,...
    'color',[0.3,0.75,0.84],'LineWidth',2);

plotCute1('Time (s)','Hx (a.u.)',ax1(1),[],{'Thorax','Abdomen'},1,'Horizontal');
ylim([-2.5,1.5])
xlim([tStart-tStart,tEnd-tStart])

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    ncs1Norm,':',...
    'color','k','LineWidth',2);
hold on 
plot(ax1(2),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    ncs2Norm,...
    'color',[0.3,0.75,0.84],'LineWidth',2);
plot(ax1(2),tResp(tStart*hxRespSampRate:tEnd*hxRespSampRate)-tStart,...
    -1+(-0.25.*detIsoVol(tStart*hxRespSampRate:tEnd*hxRespSampRate)),...
    'color',[0.9,0.75,0.2],'LineWidth',2);
plotCute1('Time (s)','NCS (a.u.)',ax1(2),[],{'Thorax','Abdomen','Isovolumetric'},1,'Horizontal');

ylim([-2.5,1.5])
xlim([tStart-tStart,tEnd-tStart])
linkaxes(ax1(:),'x')

%%
figure('Units', 'pixels', ...
    'Position', [100 100 600 700]);
% tHxResp = 0:1/hx
nFig = 4;
ax1(1) = subplot(nFig,1,1);
yyaxis left
plot(ax1(1),tResp,hxRespTh,':','color','k','LineWidth',2); %./max(hxRespTh)
plotCute2([],'Thorax (ml)',ax1(1),[],[],0);
yyaxis right
plot(ax1(1),tResp,hxRespAbd,'-','color',[0.3,0.75,0.93]); %./max(hxRespAbd)
plotCute1('Time (s)','Abdomen (ml)',ax1(1),...
    [],{'Thorax',...
    'Abdomen'},1,'Horizontal');

ax1(2) = subplot(nFig,1,2);
plot(ax1(2),tResp,ncsResp(:,1),'color','k'); 
plotCute2('Time (s)','a.u.',ax1(2),[],[],0);

ax1(3) = subplot(nFig,1,3);
plot(ax1(3),tTV,hxTV,'--','color','k','LineWidth',2);
hold on
plot(tTV,ncsTVAmp,'color',[0.3,0.75,0.93],'LineWidth',2);
xLabel = 'Time (s)';
yLabel = 'Respiratory Volume (mL)';
plotLegend = {'Hx','NCS'};
plotCute2(xLabel,yLabel,ax1(3),[],plotLegend,1,'Horizontal');

ax1(4) = subplot(nFig,1,4);
plot(ax1(4),tBR,hxBR,'--','color','k','LineWidth',2);
hold on
plot(ax1(4),tBR,ncsBR(:,1),'color',[0.3,0.75,0.93],'LineWidth',2);
xLabel = 'Time (s)';
yLabel = 'Breath rate (BPM)';
plotLegend = {'Hx','NCS'};
plotCute2(xLabel,yLabel,ax1(4),[],plotLegend,1,'Horizontal');

linkaxes(ax1,'x')

%%
ncsHeartIdx = heartAmpMinMax(heartAmpMinMax(:,2) == 0,1);
nFig = 2;
figure('Units', 'pixels', 'Position', [100 100 600 600]);

tStart = 356;
tEnd = 400;

ax1(1) = subplot(nFig,1,1);
plot(ax1(1),tHR,ncsHR(:,1),'color',[0.3,0.75,0.93],'LineWidth',2);
plotCute1('Time (s)','Heart Rate (BPM)',ax1(1),[],[],0);
    
ax1(2) = subplot(nFig,1,2);
plot(tHeart,ncsHeart(:,1),'color',[0.3,0.75,0.93],'LineWidth',1.5);

% hold on
% plot(tHeart(ncsHeartIdx),ncsHeart(ncsHeartIdx,1),'o');
plotCute1('Time (s)','Heartbeat (a.u.)',ax1(2),...
    [],[],0);
linkaxes(ax1,'x');
xlim(ax1(2),[tStart,tEnd]);


%%
ncsHeartIdx = heartAmpMinMax(heartAmpMinMax(:,2) == 0,1);
nFig = 2;
figure('Units', 'pixels', 'Position', [100 100 600 600]);

tStart = 356;
tEnd = 400;

% tStart = 20;
% tEnd = 120;

ax1(1) = subplot(nFig,1,1);
plot(ax1(1),tHR,hxHR,':','color','k','LineWidth',2);
hold on
plot(ax1(1),tHR,ncsHR(:,1),'color',[0.3,0.75,0.93],'LineWidth',2);
plotCute1('Time (s)','HR (BPM)',ax1(1),...
    [],{'ECG','NCS'},1);
hold off
    
ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(tHeart,hxEcg,'-.','color','k','LineWidth',0.5);
plotCute1('Time (s)','Hx ECG (mV)',ax1(2),[],[],0);

ylim([-0.8,1.2])
yyaxis right
plot(tHeart,ncsHeart(:,1),'color',[0.3,0.75,0.93],'LineWidth',1.5);

hold on
plot(tHeart(ncsHeartIdx),ncsHeart(ncsHeartIdx,1),'o');
plotCute1('Time (s)','Heartbeat (a.u.)',ax1(2),...
    [],{'ECG','NCS','NCS Peak'},1,'Horizontal');
hold off
linkaxes(ax1,'x');
xlim(ax1(2),[tStart,tEnd]);

%%
ncsHeartIdx = heartPhMinMax(heartPhMinMax(:,2) == 0,1);
nFig = 6;
figure('Units', 'pixels', 'Position', [100 100 400 800]);

tStart = 60;
tEnd = 200;

ax1(1) = subplot(nFig,1,1);
plot(ax1(1),tHR(tStart:tEnd),hxHR(tStart:tEnd),'color',[0.3,0.75,0.93],'LineWidth',2);
hold on
plot(ax1(1),tHR(tStart:tEnd),ncsHR(tStart:tEnd,2),':','color','k','LineWidth',1);
plotCute1([],'HR (BPM)',ax1(1),...
    [],{'ECG','NCS'},1);
hold off
    
ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(tHeart,hxEcg,'color',[0.3,0.75,0.93],'LineWidth',1);
plotCute1([],'Hx ECG (mV)',ax1(2),[],[],0);

ylim([-0.8,1.2])
yyaxis right
plot(tHeart,ncsHeart(:,2),'-.','color','k','LineWidth',1);

hold on
plot(tHeart(ncsHeartIdx),ncsHeart(ncsHeartIdx,2),'o');
plotCute1([],'Heartbeat (a.u.)',ax1(2),...
    [],{'ECG','NCS','NCS Peak'},1);
xlim(ax1(2),[115,130]);
hold off

ax1(3) = subplot(nFig,1,3);
% First plot contains Hx Abdomen and Thoracic waveforms
yyaxis left
plot(tResp,hxRespTh,':','color','k','LineWidth',2);
plotCute1([],'Thorax (ml)',ax1(3),[],[],0);
ylim(ax1(3),[2.675e4,2.683e4])

yyaxis right
plot(tResp,hxRespAbd,'color',[0.3,0.75,0.93],'LineWidth',1);
plotCute1([],'Abdomen (ml)',ax1(3),...
    [],{'Thorax',...
    'Abdomen'},1);
ylim(ax1(3), [1.708e4,1.7135e4])

% Second plot NCS respiration waveforms derived from amplitude and phase 
ax1(4) = subplot(nFig,1,4);
plot(ax1(4),tResp,ncsResp(:,1),'color','k','LineWidth',1);
plotCute1([],'NCS (a.u.)',ax1(4),...
    [],[],0);
xlim([60,200])

% Third plot Breath Rate
ax1(5) = subplot(nFig,1,5);
plot(ax1(5),tBR,hxBR,':','color','k','LineWidth',2);
hold on 
plot(ax1(5),tBR,ncsBR(:,1),'color',[0.3,0.75,0.93],'LineWidth',1);
plotCute1([],'BR (BPM)',ax1(5),...
    [],{'BR_H_x',...
    'BR_N_C_S'},1);
hold off

% Fourth plot Tidal Volume
ax1(6) = subplot(nFig,1,6);
plot(ax1(6),tTV,hxTV,':','color','k','LineWidth',2);
hold on 
plot(ax1(6),tTV,ncsTVAmpPhSum,'color',[0.3,0.75,0.93],'LineWidth',1);
plotCute1('Time (s)','TV (ml)',ax1(6),...
    [],{'TV_H_x',...
    'TV_N_C_S'},1);
hold off

xlim(ax1(1),[60,200])
xlim(ax1(3:end),[60,200])

%
% linkaxes([ax1(1),ax1(3:end)],'x')
% axis(ax1,'tight')

%% Plot NCS and Hx Resp waveforms and tidal volume: for isovolumetric
degSign = char(0176);

nFig = 3;
figure('Units', 'pixels', 'Position', [100 100 400 400]);
ax2(1) = subplot(nFig,1,1);

% First plot contains Hx Abdomen and Thoracic waveforms
yyaxis left
plot(tResp,hxRespTh,':','color','k','LineWidth',1);
plotCute1([],'Thorax (ml)',ax2(1),[],[],0);
ylim([2.65e4,2.665e4])

yyaxis right
plot(tResp,hxRespAbd,'color',[0.3,0.75,0.93],'LineWidth',1);
plotCute1([],'Abdomen (ml)',ax2(1),...
    [],{'Thorax',...
    'Abdomen'},1);

% Second plot NCS respiration waveforms derived from amplitude and phase 
ax2(2) = subplot(nFig,1,2);
plot(ax2(2),tResp,ncsResp(:,1),'color','k','LineWidth',1);
ylim([1e-3,11e-3])
plotCute1([],'NCS (a.u.)',ax2(2),...
    [],[],0);

ax2(3) = subplot(nFig,1,3);
plot(ax2(3),tTV,hxTV,':','color','k','LineWidth',1.5);
hold on
plot(ax2(3),tTV,ncsTVAmpPhSum,'color',[0.3,0.75,0.93],'LineWidth',1);
% ylim([0,600])
hold off
plotCute1('Time (s)','TV (ml)',ax2(3),...
    [],{'Hexoskin','NCS'},1);
ylim(ax2(3),[0,1000])
linkaxes(ax2,'x')
% axis(ax2,'tight')

xlim(ax2,[95,140])

%%
% Plot NCS original Amp and phase and corresponding resp and heartbeat
degSign = char(0176);

tSync = 0:1/ncsSampRate:(length(ncsSync) - 1)/ncsSampRate;

tStartSec = 2;
tEndSec = 18;
tStart = tStartSec*ncsSampRate+1; % 1
tEnd = tEndSec*ncsSampRate+1; % length(tSync)

nFig = 3;
figure('Units', 'pixels', 'Position', [100 100 400 400]);

ax4(1) = subplot(nFig,1,1);
yyaxis left
plot(ax4(1), tSync(tStart:tEnd)-tStartSec, ncsSync(tStart:tEnd,1),'color',[0.3,0.75,0.93],'LineWidth',1.5);
plotCute1([],'Amplitude (V)',ax4(1),[],[],0);
% ylim(ax4(1),[5e-3, 10e-3])

yyaxis right
plot(ax4(1), tSync(tStart:tEnd)-tStartSec, ncsSync(tStart:tEnd,2),':','color','k','LineWidth',1);
plotCute1([],['Phase (',degSign,')'],ax4(1),[],{'NCS Amplitude', 'NCS Phase'},1);

ax4(2) = subplot(nFig,1,2);
plot(ax4(2),tSync(tStart:tEnd)-tStartSec,ncsRespProcessed(tStart:tEnd,2),':','color','k','LineWidth',2);
plotCute1([],'Respiration (a.u.)',ax4(2),[],'Respiration',1);

ax4(3) = subplot(nFig,1,3);
plot(ax4(3),tSync(tStart:tEnd)-tStartSec,ncsHeartProcessed(tStart:tEnd,1),'color',[0.3,0.75,0.93]);
plotCute1('Time (s)','Heartbeat (a.u.)',ax4(3),[],'Heartbeat',1);

linkaxes(ax4,'x')
