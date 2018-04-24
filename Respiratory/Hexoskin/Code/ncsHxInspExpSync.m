% This function syncs inspiration and expiration data of NCS, that is not
% uniformly sampled like the raw data.

function [hxInhale,hxExhale]...
    = ncsHxInspExpSync(dataPath,hxFolder,hxDataNum,ncsDataNum,...
              manualTimeOffset,dataDuration,ncsTstart)
% Input:
%       dataPath: Path to experiment data folder
%       hxFolder: Relative to dataPath
%       hxDataNum: This is fixed for inspiration and expiration. Refer to
%                  readHxData.m for the folder numbers.
%       ncsDataNum: NCS data number. Refer to ncsFileInfo for the data
%                   number corresponding to the '.mat' file.
%       manualTimeOffset: Manual time offset needed in seconds.
%       dataDuration: Total data duration in seconds.
%       ncsTstart: This gives offset referred to NCS starting time
% Output:
%       hxInhale: Truncated Inhalation data arranged in columns with time
%                 information in seconds - [time, inhalation]. Inhalation 
%                 data is in arbitrary units. 
%       hxTimeTrunc: Time in seconds, corresponding to inhalation

%% ------------------------------------------------------------------------
% 