% This function automates calling all the files and performs analysis for
% each file (that is saved automatically in one place)
% Pragya Sharma, ps847@cornell.edu
% 25 Sept 2019

%% ------------------------------------------------------------------------
% Caliing ncsBioProcessFunc.m
% Reading all the saved file names

% set(gcf,'Visible','off');
% set(0,'DefaultFigureVisible','off');

[caseOrder,routineName,fileName,~,~] = case7_30Info();
numFiles = length(caseOrder);

for i = 1:numFiles
    caseNum = caseOrder(i);
    ncsBioProcessFunc(caseNum,fileName{i,1},fileName{i,2},fileName{i,3},fileName{i,4});
end
 
%%
% set(gcf,'Visible','on')
% set(0,'DefaultFigureVisible','on');
