function [shiftLocNormRates, shiftLocShiftRates, shiftLocPcInd] = cueShiftInhibPlace(cueShiftStruc, toPlot);


% Clay 2019
% Look at effect of shifted cue on place cells normally at that location
% 1. find cells with PF in location to which cue is shifted
% 2. look at mean rate in normal vs. shift laps

refLapType = findRefLapType(cueShiftStruc);
pc = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1);

normPosRates = cueShiftStruc.PCLappedSessCell{refLapType}.posRates(pc,:);

if refLapType==1
    shiftLapType=2;
else
    shiftLapType=1;
end

shiftPosRates = cueShiftStruc.PCLappedSessCell{shiftLapType}.posRates(pc,:);

[maxVal, maxInd] = max(normPosRates');
% [newInd, oldInd] = sort(maxInd);
% sortInd = oldInd;

% look at peaks of 25:40

shiftLocInd = find(maxInd>=25 & maxInd<=40); % cells (pc) with PF around shift location

shiftLocNormRates = normPosRates(shiftLocInd,:);
shiftLocShiftRates = shiftPosRates(shiftLocInd,:);

shiftLocPcInd = pc(shiftLocInd);

if toPlot
figure; 
plot(mean(shiftLocNormRates,1)); 
hold on; 
plot(mean(shiftLocShiftRates,1),'g');
plot(mean(shiftLocShiftRates,1)-mean(shiftLocNormRates,1),'r'); 
yl = ylim;
%line([25 25], yl);
end
