% This code estimates zero-crossing points of a waveform 

function zcPosNeg = zeroCrossDet(data,fs,opts)

%%
% Gives indices of zero-crossing. This WILL give error if exact 0 exists in
% the data. In that case two successive indices will be marked as
% zero-cross points. Ignore the first one.
zci = @(v) find(v(:).*circshift(v(:), -1) <= 0); 

%%
if ~isfield(opts,'minTime')
    % Minimum time between two zerocrossing points
    opts.minTime = 0;
end

zc = zci(data);
t = (0:(length(data)-1))/fs;

zcPosNeg = zeros(length(zc),2);
zcPosNeg(:,1) = zc;
zcIdx = 1;
for i = 1:length(zc)
    if zc(i) > 1
        zcPosNeg(zcIdx,2) = sign(data(zc(i)) - data(zc(i)-1));
    else
        zcPosNeg(zcIdx,2) = sign(data(zc(i)+1) - data(zc(i)));
    end
    
    if zcIdx>1
       if zcPosNeg(zcIdx,2) == zcPosNeg(zcIdx-1,2)
           % If two are same, keep the second one, throw away the first.
           
           zcPosNeg = zcPosNeg([1:zcIdx-2,zcIdx:end],:); 
           zcIdx = zcIdx-1;
           
       end
       if zcIdx > 1 % It may be decreased by previous step
           if ((zcPosNeg(zcIdx,1)-zcPosNeg(zcIdx-1,1))/fs) < opts.minTime
               
               zcPosNeg = zcPosNeg([1:zcIdx-2,zcIdx+1:end],:);
               zcIdx = zcIdx-2;
               
           end
       end
    end
    zcIdx = zcIdx + 1;
    
end

% % Plot figure
% figure
% plot(t,data); hold on
% idxPos = zcPosNeg(:,2) > 0;
% idxNeg = zcPosNeg(:,2) <= 0;
% plot (t(zcPosNeg(idxPos,1)),data(zcPosNeg(idxPos,1)),'+');
% plot (t(zcPosNeg(idxNeg,1)),data(zcPosNeg(idxNeg,1)),'o');
% legend('airflow','Inspire Start','Expire Start')
end