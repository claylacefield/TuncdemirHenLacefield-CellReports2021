function [posRatesRm1, posRatesRm2] = roomSwitchTuning()

% Clay July 2020
% Use files in /Backup20TB/clay/DGdata/cuePaper/roomSwitch
% to find all registered cells' tuning

[path] = uigetdir();
cd(path);
currDir = dir;

posRatesRm1 = [];
posRatesRm2 = [];

for i=1:4
    mouseName = currDir(i+2).name;
    cd(mouseName);
    
    refLapType = 2;
    
    if contains(mouseName, 'IR-519-P4')
        sessInd = [2 3];
        load(findLatestFilename('cueOmitShiftMint-'));
        posRates1 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
        pc1 = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC);
        load(findLatestFilename('cueOmitShiftMintRoom-'));
        posRates2 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
        pc2 = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC);
    else
        sessInd = [1 2];
        load(findLatestFilename('cueShiftOmitStatRewOlfIAA-'));
        posRates1 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
        pc1 = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC);
        load(findLatestFilename('cueShiftOmitStatRewOlfIAARoom-'));
        posRates2 = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
        pc2 = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC);
    end
    
    
    load(findLatestFilename('egist')); % cell-registered-struc
%     load(findLatestFilename('ultSessSeg')); % multSessSegStruc
%     toPlot = 0;
%     [sameCellTuningStruc] = sameCellTuning2P(multSessSegStruc, cell_registered_struct, toPlot);
%     cellsInAll = sameCellTuningStruc.cellsInAll;
    
    mapInd = cell_registered_struct.cell_to_index_map;
    mapInd = mapInd(:,sessInd);
    cellRegIndInAll = find(min(mapInd, [], 2));
    cellsInAll = mapInd(cellRegIndInAll,:);
    
    % to only take PCs?
    % XXX
    
    posRates1b = posRates1(cellsInAll(:,1),:);
    posRates2b = posRates2(cellsInAll(:,2),:);
    
%     figure;
%     subplot(2,1,1);
%     [sortInd] = plotPosRates(posRates1b, 0, 1);
%     subplot(2,1,2);
%     imagesc(posRates2b(sortInd,:));
    
    posRatesRm1 = [posRatesRm1; posRates1b];
    posRatesRm2 = [posRatesRm2; posRates2b];
    
    cd ..;
    
end

% look at corrcoefs between room1,2, by position
% from xcorrPosRates.m

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
    cc = corrcoef(posRates1s, posRates2s);
    
    [val, offset] = max(xc); % xcorr peak
    
    pkShift(i) = offset - 100;  % calc pkDiff from xcorr offset
    maxXc(i) = val;
    corCoef(i) = cc(2); % correlation coef using corrcoef
    xc0(i) = nanmean(xc(90:110));   % middle range of xcorr
    
end

pkShift(abs(pkShift)>50)=NaN; % sometimes shift doesn't exist so NaN it


%%


figure;
subplot(2,1,1);
[sortInd, pfPos] = plotPosRates(posRatesRm1, 0, 1);
subplot(2,1,2);
imagesc(posRatesRm2(sortInd,:));

%figure;
binSize = 5;
binByPos(maxXc, pkPos, binSize);
binByPos(corCoef, pkPos, binSize);