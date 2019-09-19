
function [minimaIdx, maximaIdx] = peakDet2(y,~, fs)
% peakDet2 uses some function to find peaks. 

t = (0:(length(y)-1))./fs;
x = 1:length(y);

fig1 = figure;
plot(t,y);

pkfinderPath = 'E:\NCS\HumanStudyMarch2019\Matlab\PeakFinder';
addpath(pkfinderPath);

slopeThreshold = .0001;
ampThreshold = 0.01;
smoothWidth = 200;
peakgroup = 0;
smoothType = 1;

maximaIdx = findpeaks(x,y,slopeThreshold,ampThreshold,smoothWidth,peakgroup,smoothType);
maximaIdx = round(maximaIdx(:,2));
minimaIdx = findvalleys(x,y,slopeThreshold,ampThreshold,smoothWidth,peakgroup,smoothType);
minimaIdx = round(minimaIdx(:,2));

figure(fig1);hold on;
plot(t(maximaIdx),y(maximaIdx),'^');
plot(t(minimaIdx),y(minimaIdx),'v');

rmpath(pkfinderPath);
end