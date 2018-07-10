% This code performs wavelet decomposition and does further reconstruction
% in different possible options.

%% Run postProcess3.m before this. Displays reconstructed waveforms
[ncs, ~, t] = postProcess3(0.5,1,4,5);

%% Multi-resolution wavelet transform
sampStart = 1;
nSample = length(ncs);
sampEnd = sampStart + nSample - 1;
ncsSample = ncs(sampStart:sampEnd); %Desired signal - See if negative sign needed or not

%% Wavelet decomposition
nCoeff = 9;
wNameNcs = 'bior3.5';
[Cncs,Lncs] = wavedec(ncsSample,nCoeff,wNameNcs);

%% Wavelet reconstruction, by leaving out detailed coefficients starting from d1
figure
maxSubPlot = 5;
axNum = 1;
figNum = 1;
ax(axNum) = subplot(maxSubPlot,1,figNum);
plot(ncsSample); grid on;
title('Filtered Normalized NCS','FontSize',12)
nLevel = 1;
tmpNcs = [0; cumsum(Lncs)];
CncsNew = Cncs;

for i = 1:maxSubPlot-1
    figNum = figNum + 1;
    axNum = axNum + 1;
    ax(axNum) = subplot(maxSubPlot,1,figNum);
    CncsNew(tmpNcs(end-(i-1)-2)+1:tmpNcs(end-(i-1)-1)) = 0;
    recSignal = waverec(CncsNew,Lncs,wNameNcs); 
    plot(recSignal); grid(ax(axNum),'on');
    title(['Reconstructed NCS, upto d',num2str(i),' removed'],'FontSize',12)

    nLevel = nLevel + 1;
end

figure
figNum = 0;
for i = 1:nCoeff  - maxSubPlot + 1
    figNum = figNum + 1;
    axNum = axNum + 1;
    ax(axNum) = subplot(nCoeff - maxSubPlot + 1,1,figNum);
    CncsNew(tmpNcs(end-(i-1)-2-maxSubPlot+1)+1:tmpNcs(end-(i-1)-1-maxSubPlot+1)) = 0;
    recSignal = waverec(CncsNew,Lncs,wNameNcs); 
    plot(recSignal); grid(ax(axNum),'on');
    title(['Reconstructed NCS, upto d',num2str(i + maxSubPlot - 1),' removed'],'FontSize',12)

    nLevel = nLevel + 1;
end

linkaxes(ax);
axis('tight')
ylim(ax,[-10,10])
%% Looking at instantaneous frequency of wavelet reconstruction with upto d7 removed. 
CncsNew = Cncs;
CncsNew(tmpNcs(end-2-7+1)+1:tmpNcs(end-1))=0;

ncsRec = waverec(CncsNew,Lncs,wNameNcs);

figure
ax1(1) = subplot(2,1,1);
plot(ncsSample)
grid on
legend('NCS')
title('FontSize',12)
ax1(2) = subplot(2,1,2);
plot(ncsRec)
grid on
legend('Reconstructed NCS')
title('Wavelet Reconstruction - upto d7 removed','FontSize',12)
linkaxes(ax1,'x')

fs = 512;
zNcs = hilbert(ncsRec);
instfreqNcs = fs/(2*pi)*diff(unwrap(angle(zNcs)));
avgHRncs = mean(instfreqNcs(nSample*0.1:nSample*0.9));

figure
plot([instfreqNcs,ones(nSample,1)*avgHRncs])
title('Frequency of NCS waveform')
legend('Instantaneous freq NCS','Avg freq NCS')

%% Wavelet: Viewing just reconstruction using detailed coefficients, starting from d1
figure
maxSubPlot = 5;
axNum = 1;
figNum = 1;
ax(axNum) = subplot(maxSubPlot,1,figNum);
plot(ncsSample)
title('Filtered Normalized NCS','FontSize',12)

nLevel = 1;
tmpNcs = [0; cumsum(Lncs)];

for i = 1:maxSubPlot-1
    CncsNew = zeros(length(Cncs),1);

    figNum = figNum + 1;
    axNum = axNum + 1;
    ax(axNum) = subplot(maxSubPlot,1,figNum);
    CncsNew(tmpNcs(i)+1:tmpNcs(i+1)) = Cncs(tmpNcs(i)+1:tmpNcs(i+1));
    recSignal = waverec(CncsNew,Lncs,wNameNcs); 
    plot(recSignal);
    if i==1
        title(['Ncs A',num2str(nCoeff)],'FontSize',12)
    else
        title(['Ncs D',num2str(nCoeff-i+2)],'FontSize',12)
    end

    nLevel = nLevel + 1;
end


figure
figNum = 0;
for i = 1:nCoeff  - maxSubPlot + 1
    CncsNew = zeros(length(Cncs),1);

    figNum = figNum + 1;
    axNum = axNum + 1;
    ax(axNum) = subplot(nCoeff - maxSubPlot + 1,1,figNum);
    CncsNew(tmpNcs(i+maxSubPlot-1):tmpNcs(i+maxSubPlot-1+1)) = Cncs(tmpNcs(i+maxSubPlot-1):tmpNcs(i+maxSubPlot-1+1));
    recSignal = waverec(CncsNew,Lncs,wNameNcs); 
    plot(recSignal);
    title(['Ncs D',num2str(nCoeff-i+2-maxSubPlot+1)])

    nLevel = nLevel + 1;
end

linkaxes(ax,'x');
axis('tight')
ylim(ax,[-2.5,2.5])

%% Looking at NCS D8
i = 3;
CncsNew = zeros(length(Cncs),1);
CncsNew(tmpNcs(i):tmpNcs(i+1)) = Cncs(tmpNcs(i):tmpNcs(i+1));
recSignal_Ncsd8 = waverec(CncsNew,Lncs,wNameNcs); 

figure
% subplot(2,1,1)
plot(t,ncs,t,recSignal_Ncsd8);
title(['Ncs D',num2str(nCoeff-i+2)],'FontSize',12)
legend('NCS','NCS D8')
xlabel('Samples in time')
grid on
axis 'tight'

nSampT30 = find(t>=30,1); %Index at 30sec
[r,p,RL,RU] = corrcoef(ncs(1:nSampT30),recSignal_Ncsd8(1:nSampT30));
pValue = p(1,2);
fprintf('Correlation at rest: %f',r(1,2));
fprintf('\nP-value at rest: %e\n',pValue);

% Xncsd8 = fftshift(fft(recSignal_Ncsd8));
% subplot(2,1,2)
% fs = 512;
% df = fs/nSample;
% f = -fs/2:df:fs/2 - df;
% plot(f,abs(Xncsd8)/nSample);
% title('Frequency spectrum of NCS d8','FontSize',12)
% xlabel('Frequency(Hz)')
