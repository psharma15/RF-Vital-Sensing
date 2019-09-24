function [R_i,R_amp,S_i,S_amp,T_i,T_amp,Q_i,Q_amp,heart_rate,buffer_plot]=peakdetect2(ecg,fs,view)
%% function [R_i,R_amp,S_i,S_amp,T_i,T_amp]=peakdetect(ecg,fs,view)
%% Online Adaptive QRS detector
%% Description
% QRS detection
% Detects Q , R and S waves,T Waves
% Uses the state-machine logic to determine different peaks in an ECG
% signal. It has the ability to confront noise by canceling out the noise
% by high pass filtering and baseline wander by low pass. Besides, check
% out criterion to stop detection of spikes.
% The code is written in a way for future online implementation.
%% Inputs
% ecg : raw ecg vector
% fs : sampling frequency
% view : display span of the signal e.g. 8 seconds (8 seconds is the default)

%% Outputs
% indexes and amplitudes of R_i, R_amp, etc
% heart_rate computed heart rate
% buffer_plot : processed signal
%% how to use
% for example after loading the the ecg mat files in matlab call the
% function as below ;
% [R_i,R_amp,S_i,S_amp,T_i,T_amp]=peakdetect(EKG1,250,10);


%% Author : Hooman Sedghamiz    contact :hoose792@student.liu.se , hooman650@yahoo.com
% Dont forget to reference if you found this script usefull

%% Changes
% Updated by: Pragya Sharma March 16, 2018
% Removing windowing to smooth the signal, as final peaks are a bit
% shifted, when applied on the Raw ECG waveform.

%%
close all ; % close the active plots
if nargin < 3
    view = 8; % on default the first 8 seconds are viewed
    if nargin <2
       fs = 250; %default Sampling frequency
       if nargin < 1
           [FileName,PathName] = uigetfile('*.mat');
           localdir = dir;
           cd(PathName);
           load(FileName);
           cd(localdir);
           ecg = EKG1; % on default the program uses EKG
       end
    end
    
end

%% initialize
R_i = [];%save index of R wave
R_amp = []; %save amp of R wave
S_i = [];%save index of S wave
S_amp = []; %save amp of S wave
T_i = [];%save index of T wave
T_amp = [];%save amp of T wave
thres_p =[]; %for plotting adaptive threshold
buffer_plot =[];
buffer_long=[]; % buffer for online processing
state = 0 ; % determines the state of the machine in the algorithm
c = 0; % counter to determine that the state-machine doesnt get stock in T wave detection wave
T_on = 0; % counter showing for how many samples the signal stayed above T wave threshold
T_on1=0; % counter to make sure its the real onset of T wave
S_on = 0; % counter to make sure its the real onset of S wave
sleep = 0; % counter that avoids the detection of several R waves in a short time
S_amp1 = []; % buffer to set the adaptive T wave onset
buffer_base=[]; %buffer to determine online adaptive mean of the signal
dum = 0; %counter for detecting the exact R wave
% window = round(fs/25); % averaging window size
window = 0; % averaging window size
weight = 1.8; %initial value of the weigth
co = 0; % T wave counter to come out of state after a certain time
thres2_p = []; %T wave threshold indices
thres_p_i = []; %to save indices of main thres
S_amp1_i = []; %to save indices of S thres
thres2_p_i = []; %to save indices of T threshold
Q_i = []; % vectors to store Q wave
Q_amp =[]; %vectors to store Q wave
%% preprocess
ecg = ecg (:); % make sure its a vector
ecg_raw =ecg; %take the raw signal for plotting later
time_scale = length(ecg_raw)/fs; % total time;
%Noise cancelation(Filtering)
f1=0.5; %cuttoff low frequency to get rid of baseline wander
f2=45; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg = filtfilt(a,b,ecg);

%% define two buffers

buffer_mean=mean(abs(ecg(1:2*fs)-mean(ecg(1:2*fs)))); % adaptive threshold DC corrected (baseline removed)
buffer_T = mean(ecg(1:2*fs));%second adaptive threshold to be used for T wave detection
%% start online inference (Assuming the signal is being acquired online)
for i = 1 : length(ecg)
    
buffer_long = [buffer_long ecg(i)] ; % save the upcoming new samples
buffer_base = [buffer_base ecg(i)] ; % save the baseline samples

%% Renew the mean and adapt it to the signal after 1 second of processing
if length(buffer_base) >= 2*fs
    buffer_mean = mean(abs(buffer_base(1:2*fs)-mean(buffer_base(1:2*fs))));
    buffer_T = mean(buffer_base(1:2*fs));
    buffer_base =[];
end

%% smooth the signal by taking the average of 15 samples and add the new upcoming samples
  if length(buffer_long)>= window % take a window with length 15 samples for averaging
      mean_online = mean(buffer_long);  % take the mean
      buffer_plot =[buffer_plot mean_online]; % save the processed signal
      
      
    %% Enter state 1(putative R wave) as soon as that the mean exceeds the double time of threshold  
    if state == 0  
     if length(buffer_plot) >= 3   %added to handle bugg for now
      if mean_online > buffer_mean*weight && buffer_plot(i-1-window) > buffer_plot(i-window)    %2.4*buffer_mean   
          state = 1; % entered R peak detection mode
          currentmax = buffer_plot(i-1-window);
          ind = i-1-window;
          thres_p = [thres_p buffer_mean*weight];
          thres_p_i = [thres_p_i ind];
      else     
          state = 0;
      end
     end
    end
    
    %% Locate the R wave location by finding the highest local maxima
      if state == 1 % look for the highest peak
          
            if  currentmax > buffer_plot(i-window)
                dum = dum + 1;
                if dum > 4 
                R_i = [R_i ind];%save index
                R_amp = [R_amp buffer_plot(ind)]; %save index
                % Locate Q wave
                [Q_tamp Q_ti] = min(buffer_plot(ind-round(0.040*fs):(ind)));
                Q_ti = ind-round(0.040*fs) + Q_ti -1;
                Q_i = [Q_i Q_ti];
                Q_amp = [Q_amp Q_tamp];
                
                
                if length(R_amp) > 8
                weight = 0.30*mean(R_amp(end-7:end)); %calculate the 35% of the last 8 R waves
                weight = weight/buffer_mean;
                end
                state = 2; % enter S detection mode state 2
                dum = 0;
                end
            else
                dum = 0;
                state = 0;
            end 
            
      end
      
    %% check weather the signal drops below the threshold to look for S wave
      if state == 2 
        if  mean_online <= buffer_mean     % check the threshold
             state = 3;   %enter S detection           
        end
      end
      
      %% Enter S wave detection state3 (S detection)
          if state == 3
            co = co + 1; 
            
          if co < round(0.200*fs)
            if buffer_plot(i-window-1) <= buffer_plot(i-window) % see when the slope changes
             S_on = S_on + 1; % set a counter to see if its a real change or just noise
             if S_on >= round(0.0120*fs)
             S_i = [S_i i-window-4];%save index of S wave
             S_amp = [S_amp buffer_plot(i-window-4)];%save index
             S_amp1 = [S_amp1  buffer_plot(i-window-4)]; %ecg(i-4)
             S_amp1_i = [S_amp1_i ind]; %index of S_amp1_i
             state = 4; % enter T detection mode
             S_on = 0;
             co = 0;
             end
            end
          else
              state = 4;
              co = 0;
          end
          end
      
       %% enter state 4 possible T wave detection
       if state == 4    
         if mean_online < buffer_mean % see if the signal drops below mean 
           state = 6; % confirm
         end
       end
       %% Enter state 6 which is T wave possible detection
       if state ==6   
         c = c + 1; % set a counter to exit the state if no T wave detected after 0.3 second
         if c <= 0.7*fs  
             % set a double threshold based on the last detected S wave and
             % baseline of the signal and look for T wave in between these
             % two threshold
             thres2 = ((abs(abs(buffer_T)-abs(S_amp1(end))))*3/4 + S_amp1(end)); 
             thres2_p =[thres2_p thres2];
             thres2_p_i =[thres2_p_i ind];
             if mean_online > thres2
              T_on = T_on +1; % make sure it stays on for at least 3 samples
              if T_on >= round(0.0120*fs)
               if buffer_plot(i-window-1)>= buffer_plot(i-window)
                   T_on1 = T_on1+1; % make sure its a real slope change
                  if T_on1 > round(0.0320*fs) 
                   T_i = [T_i i-window-11];%save index of T wave
                   T_amp = [T_amp  buffer_plot(i-window-11)];%save index
                   state = 5; % enter sleep mode
                   T_on = 0;
                   T_on1 = 0;
                  end
                                      
               end
              end
             end
             
            
          
         else
             state= 5; % enter Sleep mode
         end
         
       end  
      
        
       %% this state is for avoiding the detection of a highly variate noise or another peak
       % this avoids detection of two peaks R waves less than half a second
       if state==5
           sleep =sleep+c+1;
           c = 0;
           if sleep/fs >= 0.400
               state = 0;
               sleep = 0;%look for the next peak
           end  
       end
      
      % update the online buffer by removing the oldest sample
      buffer_long(1)=[];
      
      
  end
  
end

%% conditions
R_R = diff(R_i); % calculate the distance between each R wave
heart_rate=length(R_i)/(time_scale/60); % calculate heart rate
% msgbox(strcat('Heart-rate is = ',mat2str(heart_rate)));

% compute the min max R-R wave
max_R_interval = max(R_R);
min_R_interval = min(R_R);
% detect arythmia if there is any irregularity
% if (max_R_interval/fs)-(min_R_interval/fs) > 0.16
%     msgbox('Irregular Rhythm','Arrythmia Detected');
% end


%% plottings 
time = 1/fs:1/fs:view;
R = find(R_i <= view*fs); % determine the length for plotting vectors
S = find(S_i <= view*fs); % determine the length for plotting vectors
T = find(T_i <= view*fs); % determine the length for plotting vectors
Q = find(Q_i <= view*fs); % determine the length for plotting vectors
L1 = find(thres_p_i <= view*fs);
L2 = find(S_amp1_i <= view*fs);
L3 = find(thres2_p_i <= view*fs);
if view*fs > length(buffer_plot)
ax(1) = subplot(211);plot(time(1:length(buffer_plot)),buffer_plot(1:end));   
else
ax(1) = subplot(211);plot(time,buffer_plot(1:(view*fs)));
end
hold on,scatter(R_i(1:R(end))./fs,R_amp(1:R(end)),'r');
hold on,scatter(S_i(1:S(end))./fs,S_amp(1:S(end)),'g');
hold on,scatter(T_i(1:T(end))./fs,T_amp(1:T(end)),'k');
hold on,scatter(Q_i(1:Q(end))./fs,Q_amp(1:Q(end)),'m');
hold on,plot(thres_p_i(1:L1(end))./fs,thres_p(1:L1(end)),'LineStyle','-.','color','r',...
    'LineWidth',2.5);
hold on,plot(S_amp1_i(1:L2(end))./fs,S_amp1(1:L2(end)),'LineStyle','--','color','c',...
    'LineWidth',2.5);
hold on,plot(thres2_p_i(1:L3(end))./fs,thres2_p(1:L3(end)),'-k','LineWidth',2);
legend('Raw ECG Signal','R wave','S wave','T wave','R adaptive thres','Latest S wave','T wave adaptive threshold threshold','Location','NorthOutside','Orientation','horizontal'); 
xlabel('Time(sec)'),ylabel('V');
axis tight; title('Zoom in to see both signal details overlaied');title('Filtered, smoothed and processed signal');
ax(2) =subplot(212);plot(time,ecg_raw(1:(view*fs)));title('Raw ECG')
xlabel('Time(sec)'),ylabel('V');
legend(); 
linkaxes(ax,'x');
zoom on;
% close all;


    
