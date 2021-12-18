function [pkShiftPC, pkPosPC] = xcorrPosRates(cueShiftStruc, toPlot)

%% USAGE: [pkShift, pkPosSeg] = xcorrPosRates(cueShiftStruc, toPlot);
% NOTE: posRates1 should always be normal laps, and posRates2 shift

filename = cueShiftStruc.filename;
[refLapType] = findRefLapType(cueShiftStruc);
posRates1 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
if refLapType==1; shiftLap=2; else shiftLap=1; end
posRates2 = cueShiftStruc.PCLappedSessCell{shiftLap}.posRates;


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
    
    %% xcorr centered posRates
    xc = xcorr(posRates1s, posRates2s);  % xcorr centered
    
    [val, offset] = max(xc); % xcorr peak
    
    pkShift(i) = offset - 100;  % calc pkDiff from xcorr offset
    
end

pkShift(abs(pkShift)>50)=NaN; % sometimes shift doesn't exist so NaN it

pc = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1);
pkPosPC = pkPos(pc); pkShiftPC = pkShift(pc); % just take PCs

if toPlot==1
       
    figure; plot(pkPosPC, pkShiftPC,'.');
    try title(filename); catch; end
    xlabel('pkPos'); ylabel('pkShift');
    
    
    nBins = 20;
    [N, edges, bins] = histcounts(pkPosPC,nBins);
    for j=1:nBins
        binPkShift = pkShiftPC(bins==j);
        avShift(j) = nanmean(binPkShift);
        semShift(j) = std(binPkShift)/sqrt(length(avShift));
    end
    figure; bar(avShift); 
    xlabel('posBin'); ylabel('xcorr pk shift (bins out of 100)');
    try title(filename); catch; end
end
