% This function estimates fractional inspiratory time (Ti/Tt, ratio of
% inspiratory to total breath time) using Hexoskin Inspiration and
% Expiration data.
% Inhalation: Defined as end of inhalation
% Exhalation: Defined as end of exhalation
% April 24, 2018
% Pragya Sharma, ps847@cornell.edu

function hxTiTt = hxEstTiTt(tHxInsp,tHxExp,hxInsp,hxExp)
% Input:
%       hxInsp: Inspiration amplitude, most likely from the respiratory
%       waveform calculated using ponderate sum of thoracic and abdominal.
%       hxExp: Expiration amplitude, most likely from the respiratory
%       waveform calculated using ponderate sum of thoracic and abdominal.
%       tAbsHxInsp: Date Time vector corresponding to each Inspiration. Its
%       columns are [yyyy mm dd HH MM SS.ffff]
%       tAbsHxExp: Date Time vector corresponding to each Expiration. Its
%       columns are [yyyy mm dd HH MM SS.ffff]
% Output:
%       hxTiTt: First coumn is inspiration time. Second is corresponding 
%       calculated Ti/Tt value

%% ------------------------------------------------------------------------
% This is Hexoskin data. Assuming, this has one inspiration and one
% expiration in one breath cycle. Correcting such that it starts with an 
% inspiration (expiration minima) and ends with an expiration (expiration
% minima).

if hxInsp(1) < hxExp(1)
    fprintf(['Hexoskin end of inhalation amplitude is less than end of ',...
        'exhalation amplitude. Check if issue with data or need to ',...
        'exchange inhalation and exhalation data. Code stopping...\n']);
    return
end

if ((length(tHxInsp) > length(tHxExp) + 1) || (length(tHxInsp) < length(tHxExp) - 1))
    fprintf(['Hexoskin data has incorrect number of inspiration and ',...
        'expiration detection. Code stopping... \n']);
    return
end

if tHxInsp(1) < tHxExp(1)
    % First point is beginning of expiration. Ignore this point.
    tHxInsp = tHxInsp(2:end);
end
if length(tHxInsp) >= length(tHxExp)
    % End point is beginning of expiration. Ignore this point.
    tHxInsp = tHxInsp(1:(length(tHxExp)-1));
end

tI = tHxInsp - tHxExp(1:end-1);
tT = diff(tHxExp);

hxTiTt = [tHxInsp, tI./tT];

