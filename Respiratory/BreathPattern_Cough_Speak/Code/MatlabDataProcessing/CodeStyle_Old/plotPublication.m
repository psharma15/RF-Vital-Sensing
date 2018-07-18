% Result plot: This code plots the final figures for the paper. Results are
% obtained from ncsHxCompare.m

%% 
figure('Units', 'pixels', ...
    'Position', [100 100 1200 700]);
% tHxResp = 0:1/hx
nFig = 3;
ax1(1) = subplot(nFig,1,1);
plot(ax1(1),tResp,hxRespTh./max(hxRespTh),':','color','k','LineWidth',2); %./max(hxRespTh)
hold on
plot(ax1(1),tResp,hxRespAbd./max(hxRespAbd),'-','color',[0.3,0.75,0.93]); %./max(hxRespAbd)
plotCute1('Time (sec)','a.u.',ax1(1),...
    'Hexoskin Respiration (thorax & abdomen)',{'Thorax (a.u.)',...
    'Abdomen (a.u.)'},1);

ax1(2) = subplot(nFig,1,2);
yyaxis left
plot(ax1(2),tResp,ncsResp(:,1)./max(ncsResp(:,1))); 
plotCute1([],'a.u.',ax1(2),[],[],0);
% ylim(ax1(2),[0.58,1.01])

yyaxis right
plot(ax1(2),tResp,ncsResp(:,2)./max(ncsResp(:,2))); 
plotCute1('Time (sec)','a.u.',ax1(2),...
    'NCS Respiration (amplitude & phase)',{'NCS Amp (a.u.)','NCS Ph (a.u.)'},1);
% ylim([0.6,1.001])


ax1(3) = subplot(nFig,1,3);
plot(ax1(3),tTV,hxTV,'color',[0 0.9 0],'LineWidth',2);
hold on
plot(tTV,ncsTVPh,'--','color',[0.9 0.2 0.2]);

hold off
grid on

xLabel = 'Time (sec)';
yLabel = 'Tidal Volume (mL)';
plotTitle = 'Estimated TV from Hexoskin and calibrated TV from NCS';

plotLegend = {'Hx TV','NCS TV: D*ph'};
plotCute1(xLabel,yLabel,ax1(3),plotTitle,plotLegend,1);
axis(ax1(3),'tight')

linkaxes(ax1,'x')
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

tSync = 0:1/ncsHighSampRate:(length(ncsSync) - 1)/ncsHighSampRate;

nFig = 3;
figure('Units', 'pixels', 'Position', [100 100 400 400]);

ax4(1) = subplot(nFig,1,1);
yyaxis left
plot(ax4(1), tSync, ncsSync(:,1),'color',[0.3,0.75,0.93],'LineWidth',1.5);
plotCute1([],'Amplitude (V)',ax4(1),[],[],0);
ylim(ax4(1),[5e-3, 10e-3])

yyaxis right
plot(ax4(1), tSync, ncsSync(:,2),':','color','k','LineWidth',2);
plotCute1([],['Phase (',degSign,')'],ax4(1),[],{'NCS Amplitude', 'NCS Phase'},1);

ax4(2) = subplot(nFig,1,2);
plot(ax4(2),tResp,ncsResp(:,2),':','color','k','LineWidth',2);
plotCute1([],'Respiration (a.u.)',ax4(2),[],[],0);

ax4(3) = subplot(nFig,1,3);
plot(ax4(3),tHeart,ncsHeart(:,1));
plotCute1('Time (s)','Heartbeat (a.u.)',ax4(3),[],[],0);
