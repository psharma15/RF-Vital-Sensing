function [ncsReconstr,pCorr,pValue] = wavedecrec1(ncs,fs,wName,nCoeff,dispOn)
%%WAVEDECREC1 function performs wavelet decomposition and reconstruction of
%%ncs. Input needs ncs waveform, time, wavelet name, and which detailed
%%coefficient needs to be reconstructed. Output is reconstructed NCS
%%corresponding to reconstruction coefficient, Pearson's correlation and
%%P-value.
% Example:
%   [recNcsd8,pCorr,pValue] = wavedecrec1(ncs,t,'bior5.5',8,0);

[Cncs,Lncs] = wavedec(ncs,nCoeff,wName);
tmpNcs = [0; cumsum(Lncs)];

% Reconstructing detailed coefficient 'reconstructionCoeff'
i = 2; % First is approximation coeff cA_<last_decomp_level>
CncsNew = zeros(length(Cncs),1);
% Second is last detail coeff cA_<last_decomp_level>. Only keeping this for
% reconstruction.
CncsNew(tmpNcs(i):tmpNcs(i+1)) = Cncs(tmpNcs(i):tmpNcs(i+1));
ncsReconstr = waverec(CncsNew,Lncs,wName); 

t = 0:1/fs:(length(ncs)-1)/fs;
tStart = t(1);

nSampT60 = find(t>=60+tStart,1); %Index at t=30 sec
[rRest,pRest,~,~] = corrcoef(ncs(1:nSampT60),ncsReconstr(1:nSampT60));

% nSampT65 = find(t>=65+tStart,1); %Index at t=65 sec
% nSampT95 = find(t>=95+tStart,1); 
% [rMotion,pMotion,~,~] = corrcoef(ncs(nSampT65:nSampT95),reconstructedNcs(nSampT65:nSampT95));

% pCorr = [rRest(1,2), rMotion(1,2)]; %Pearson's correlation at rest and motion
% pValue = [pRest(1,2), pMotion(1,2)]; %P-value at rest and motion

pCorr = rRest(1,2); %Pearson's correlation at rest and motion
pValue = pRest(1,2); %P-value at rest and motion

if dispOn ==1
    figure
    plot(t,ncs,t,ncsReconstr);
    title(['Ncs D',num2str(nCoeff)],'FontSize',12)
    legend('Input NCS',['NCS D',num2str(nCoeff)])
    xlabel('Samples in time')
    grid on
    axis 'tight'
    
    fprintf('Correlation for first 60 sec: %f',rRest(1,2));
    fprintf('\nP-value for first 60 sec: %e\n', pRest(1,2));
end

end