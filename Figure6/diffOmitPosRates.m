function [dPFratePC, pkPosPC, relPFratePC] = diffOmitPosRates(cueShiftStruc, toPlot)

%% USAGE: [dPFratePC, pkPosPC, relPFratePC] = diffOmitPosRates(cueShiftStruc, toPlot)
% NOTE: posRates1 should always be normal laps, and posRates2 shift
% mod from: xcorrPosRates
%
% Clay 2019
% Compare posRates for all cells for normal cue vs. omit laps
% To make plot of how much cells are affected by cue omission, by position


filename = cueShiftStruc.filename;
[refLapType] = findRefLapType(cueShiftStruc);
posRates1 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates; % normal laps
posRates2 = cueShiftStruc.PCLappedSessCell{end}.posRates; % omit laps should always be last


for i = 1:size(posRates1,1) % for all cells
    %% center pk to middle (mainly helps w edge cells)
    % (prob not neces for same session?)
    % shift shiftCue laps the same amount as normal ones
    
    [val, pkPosSeg] = max(posRates1(i,:)); % find peak rate pos
    pkPos(i) = pkPosSeg;
    if pkPosSeg<50
        posRates1s = circshift(posRates1(i,:), 50-pkPosSeg); % and shift to middle
        posRates2s = circshift(posRates2(i,:), 50-pkPosSeg);
    else
        posRates1s = circshift(posRates1(i,:), -(pkPosSeg-50));
        posRates2s = circshift(posRates2(i,:), -(pkPosSeg-50));
    end
    
%     figure; hold on;
%     plot(posRates1(i,:)); plot(posRates2(i,:),'r');
%     plot(posRates1s,'c'); plot(posRates2s,'m');
    
    %% diff in centered posRates (bet omit and normal)
    dPFrate(i) = mean(posRates1s(40:60)-posRates2s(40:60));  % diff PF rate centered, norm - omit
    relPFrate(i) = mean(posRates2s(40:60))/mean(posRates1s(40:60)); %
    
end

pc = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1);
pkPosPC = pkPos(pc); 
dPFratePC = dPFrate(pc); % just take PCs
relPFratePC = relPFrate(pc);

if toPlot==1
       
    figure; plot(pkPosPC, dPFratePC,'.');
    try title(filename); catch; end
    xlabel('pkPos'); ylabel('dPFrate');
    
    
    nBins = 20;
    [N, edges, bins] = histcounts(pkPosPC,nBins);
    for j=1:nBins
        binDiff = dPFratePC(bins==j);
        avDiff(j) = nanmean(binDiff);
        semDiff(j) = std(binDiff)/sqrt(length(avDiff));
    end
    figure; bar(avDiff); 
    xlabel('posBin'); ylabel('xcorr pk shift (bins out of 100)');
    try title(filename); catch; end
end
