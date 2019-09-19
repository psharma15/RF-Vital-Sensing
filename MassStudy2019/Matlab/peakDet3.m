% Pragya Sharma, ps847@cornell.edu
% Modified: 19th May 2019 (Created: Feb 13th 2019)

function [minimaIdx, maximaIdx] = peakDet3(y, fs,opts)
% Performing peak detection.
% Based on: W. Lu, "A semi-automatic method for peak and valley
% detection in free-breathing respiratory waveforms."


%% Checking options
if ~isfield(opts,'tWin')
    opts.tWin = 2; % In seconds
    fprintf('The default window for peak detection is %d.\n',opts.tWin);
end
if ~isfield(opts,'minInterceptDist')
    opts.minInterceptDist = 0.05; % minimum time (s) between two intercepts
    fprintf('The default minimum intercept distance is %3.2f. \n',opts.minInterceptDist);
end

%% ------------------------------------------------------------------------
t = (0:(length(y)-1))./fs;

tWin = opts.tWin; % Window time length
mac = zeros(length(y),1);
nSampWin = tWin * fs;
%fLocal = [];
stopIdx = length(y);
T = [];
nSampAvg = 3; % For detecting up and down intercepts, use average of samples of length nAvgWin
tWait = 2*fs; % FFT is costly, so only do it every few seconds?

idxUpDownInt = zeros(stopIdx,2); % Up Intercept: +1, Down Intercept: -1
countInt = 1;
for i = 1:stopIdx
   
    ti = t(i);
    
    % Step 1: Find T but only every tWait
    if mod(i-1,tWait) == 0
        if i <= nSampWin
            yWin = y(1:nSampWin);
        elseif (i > nSampWin) &&(i < stopIdx-nSampWin)
            yWin = y(i-nSampWin+1:i);
        else
            yWin = y(stopIdx-nSampWin+1:stopIdx);
        end
        Y = fft(yWin);
        P2 = abs(Y/length(yWin));
        P1 = P2(1:length(yWin)/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:(length(yWin)/2))/length(yWin);
        [~,idxPk] = max(P1);
        %fLocal = [fLocal; f(idxPk)];
        if( f(idxPk) < 1.2) && (f(idxPk) >= 0.08)
            T = [T; 1/f(idxPk)]; % In seconds
        elseif (f(idxPk) < 0.08)
            T = [T; tWin];
        else
            T = [T; 1/1.2];
        end
    else
        T = [T; T(i-1)];
        %fLocal = [fLocal; fLocal(i-1)];
    end
    nSampT = ceil(T(i) * fs);

    % Step 2: Calculation of moving average
    if (0<=ti) && (ti<=T(i))
        mac(i) = mean(y(1:2*nSampT));
    elseif (T(i)<ti) && (ti<=t(stopIdx)-T(i))
        mac(i) = mean(y(i-nSampT+1:i+nSampT));
    else
        mac(i) = mean(y(stopIdx-2*nSampT+1:stopIdx));
    end
    % Step 3: Calculation of up and down intercepts
    if (i>nSampAvg) && (i<stopIdx-nSampAvg)
        if (mean(y(i-nSampAvg:i-1)) <= mac(i-1)) && (mean(y(i:i+nSampAvg-1)) >= mac(i))
            % Up intercept
            idxUpDownInt(countInt,:) =  [i,+1]; countInt = countInt + 1;
        elseif (mean(y(i-nSampAvg:i-1)) >= mac(i-1)) && (mean(y(i:i+nSampAvg-1)) <= mac(i))
            % Down intercept
            idxUpDownInt(countInt,:) = [i,-1]; countInt = countInt + 1;
        end    
    elseif i>1
        if (y(i-1) <= mac(i-1)) && (y(i) >= mac(i))
            % Up intercept 
            idxUpDownInt(countInt,:) = [i,+1]; countInt = countInt + 1;
        elseif (y(i-1) >= mac(i-1)) && (y(i) <= mac(i))
            idxUpDownInt(countInt,:) = [i,-1]; countInt = countInt + 1;
        end
    end        
    
end

countInt = countInt - 1;
idxUpDownInt = idxUpDownInt(1:countInt,:);

upIntIdx = idxUpDownInt(idxUpDownInt(:,2)==1,1); % Up intercept index
downIntIdx = idxUpDownInt(idxUpDownInt(:,2)==-1,1); % Down intercept index

fig1 = figure;
ax(1) = subplot(2,1,1);
plot(t,y);
hold on
plot(t,mac);
plot(t(upIntIdx),y(upIntIdx),'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(t(downIntIdx ),y(downIntIdx),'v','MarkerEdgeColor','r','MarkerFaceColor','r');

%% ------------------------------------------------------------------------
% Step 4: Correcting the intercept and finding peaks between them

minimaIdx = [];
maximaIdx = [];
counter = 1;
%fig2 = figure;
while counter < length(idxUpDownInt)
    if idxUpDownInt(counter,2) == 1
        % Up intercept
        if idxUpDownInt(counter+1,2) == 1
            % If 2 consecutive up intercept, skip.
            counter = counter + 1; 
        else
            % If next is down intercept, check if there is minimum number
            % of points in between two. Otherwise skip.
            if idxUpDownInt(counter+1,1)-idxUpDownInt(counter,1) >= (opts.minInterceptDist*fs)
                % Valid intercepts. Find maxima between them
                [~,idx] = max(y(idxUpDownInt(counter,1):idxUpDownInt(counter+1,1)));
                maximaIdx = [maximaIdx; idx + idxUpDownInt(counter,1) - 1]; %#ok<AGROW
            end
            counter = counter + 1;
        end
    else
        % Down intercept
        if idxUpDownInt(counter+1,2) == 1
            % If next is up intercept, check if there is minimum number of
            % points in between two. Otherwise skip.
            if idxUpDownInt(counter+1,1)-idxUpDownInt(counter,1) >= (opts.minInterceptDist*fs)
                % Valid intercepts. Find minima between them
                %figure(fig2);plot(y(idxUpDownInt(counter,1):idxUpDownInt(counter+1,1)));
                [~,idx] = min(y(idxUpDownInt(counter,1):idxUpDownInt(counter+1,1)));
                minimaIdx = [minimaIdx; idx + idxUpDownInt(counter,1) - 1]; %#ok<AGROW
            end
            counter = counter + 1;
        else
            % If next is also down intercept, skip.
            counter = counter + 1;
        end
    end
end
                
figure(fig1)
ax(2) = subplot(2,1,2);
plot(t,y);
hold on
plot(t,mac);
plot(t(maximaIdx),y(maximaIdx),'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(t(minimaIdx),y(minimaIdx),'v','MarkerEdgeColor','r','MarkerFaceColor','r');
% ax(2) = subplot(2,1,2);
% plot(t,fLocal)
linkaxes(ax)

                
    

