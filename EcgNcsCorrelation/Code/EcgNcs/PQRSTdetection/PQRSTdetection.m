% Source: https://www.mathworks.com/matlabcentral/fileexchange/66098-ecg-p-qrs-t-wave-detecting-matlab-code

sig=load('ecg_60hz_200.dat');
fs=200;

% dataPath = 'D:\Research\SummerFall17Spring18\CnC\NCS\EcgNcsCorrelation\CodeAndData\Data\Mar03';
% addpath(dataPath);
% load('freq2G2.mat');
% rmpath(dataPath);
% ecgData = freq2G2(:,3);
% sig = ecgData;
% fs = 512;

N=length(sig);
t=[0:N-1]/fs;
figure(1);subplot(4,2,1);plot(sig)
title('Original Signal')

%%
     %           Low Pass Filter

b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];
sigL=filter(b,a,sig);
subplot(4,2,3);plot(sigL)
title('Low Pass Filter')
subplot(4,2,4);zplane(b,a)

%%
     %           High Pass Filter

b=[-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32];
a=[1 -1];
sigH=filter(b,a,sigL);
subplot(4,2,5);plot(sigH)
title('High Pass Filter')
subplot(4,2,6);zplane(b,a)

%%
     %          Derivative Base Filter

b=[1/4 1/8 0 -1/8 -1/4];
a=[1];
sigD=filter(b,a,sigH);
subplot(4,2,7);plot(sigD)
title('Derivative Base Filter')
subplot(4,2,8);zplane(b,a)

%%
     %      be tavane 2 miresanim
sigD2=sigD.^2;

%%
     %      normalization
signorm=sigD2/max(abs(sigD2));

%%
  
h=ones(1,31)/31;
sigAV=conv(signorm,h);
sigAV=sigAV(15+[1:N]);
sigAV=sigAV/max(abs(sigAV));
figure(2);plot(sigAV)
title('Moving Average filter')

%%
treshold=mean(sigAV);
P_G= (sigAV>0.01);
figure(3);plot(P_G)
title('treshold Signal')
figure;plot(sigL)
%%
difsig=diff(P_G);
left=find(difsig==1);
raight=find(difsig==-1);

%%
     %      run cancel delay
     %      6 sample delay because of LowPass filtering
     %      16 sample delay because of HighPass filtering
left=left-(6+16);
raight=raight-(6+16);

%%
    % P-QRS-t
for i=1:length(left)
   
    [R_A(i) R_t(i)]=max(sigL(left(i):raight(i)));
    R_t(i)=R_t(i)-1+left(i) %add offset
   
    [Q_A(i) Q_t(i)]=min(sigL(left(i):R_t(i)));
    Q_t(i)=Q_t(i)-1+left(i)
  
    [S_A(i) S_t(i)]=min(sigL(left(i):raight(i)));
    S_t(i)=S_t(i)-1+left(i)
    
    [P_A(i) P_t(i)]=max(sigL(left(i):Q_t(i)));
    P_t(i)=P_t(i)-1+left(i)
    
    [T_A(i) T_t(i)]=max(sigL(S_t(i):raight(i)));
    T_t(i)=T_t(i)-1+left(i)+47
    
   
end

%%



figure;plot(t,sigL,t(Q_t),Q_A,'*g',t(S_t),S_A,'^k',t(R_t),R_A,'ob',t(P_t),P_A,'+b',t(T_t),T_A,'+r');
for i=1:((length(P_t))-1)
    
    HRV=P_t(i+1)-P_t(i)
end

    