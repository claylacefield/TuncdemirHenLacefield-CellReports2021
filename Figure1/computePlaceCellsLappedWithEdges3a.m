function PCLappedSess = computePlaceCellsLappedWithEdges3(spikes, treadPos, T, shuffN)
%function [PCLappedSess, PCSimple] = computePlaceCellsLappedWithEdges2(spikes, treadPos, T, shuffN);
%written by Andres Grosmark, towards the end of 2017, and also the
%beginning of 2018.
%spikes = rows are samples, and columns are neurons

if size(spikes,2)>size(spikes,1)
    spikes = spikes';
end

%minimum percentage of laps in which cell must fire above its mean rate, also
%minumum number of laps that a bin must be occupied in order to be included
%in the analyisis
lapPerc = 10; 
lapPerc = lapPerc/100;

%absolute minimum number of laps (for sessions with few laps)
lapMin = 2;

%calculate which lap is which:
[lapVec, ~] = calcLaps1(treadPos, T);

%run the place field (with 0 shuffles) once to find the initial valid
%time-points
PCSimple = computePlaceTransVectorLapCircShuffWithEdges4(spikes, treadPos, T, lapVec, 0);

nanPos = PCSimple.nanedPos;
lapVec(isnan(nanPos)) = NaN;
lapU = unique(lapVec(~isnan(lapVec)));


%calculate the occupancy for each lap
posRateByLap = zeros(size(spikes, 2), 100, length(lapU));
posRateRawByLap = zeros(size(spikes, 2), 100, length(lapU));
rawOccuByLap = [];
OccuByLap = [];
for i = 1:length(lapU)
    notLap = lapVec ~= lapU(i);
    pc = computePlaceTransVectorLapCircShuffWithEdges4(spikes, treadPos, T, lapVec, 0, notLap);
    posRateByLap(:, :, i) = pc.posRates;
    posRateRawByLap(:, :, i) = pc.rawPosRates;
    rawOccuByLap = [rawOccuByLap; pc.rawOccupancy];
    OccuByLap = [OccuByLap; pc.Occupancy];
end

%calculate the minumum lap number:
nLap = size(OccuByLap, 1);
minL = round(nLap*lapPerc);
if minL < lapMin
    minL = lapMin;
end
badBins = find(sum(rawOccuByLap > 0, 1) < minL);
posBins = linspace(0, 1, 101);
posBins(end) = posBins(end) + 0.01;
[~, whichBin] = histc(treadPos, posBins);
badPos = ismember(whichBin, badBins);

%run the main place code itself, this time excluding 'bad' bins (bins for
%which the occupancy is not > 0 for at least minL of laps)
[PCSimple, ratesAllShuff] = computePlaceTransVectorLapCircShuffWithEdges4(spikes, treadPos, T, lapVec, shuffN, badPos);

PCSimple.ByLap.posRateByLap = posRateByLap;
PCSimple.ByLap.posRateRawByLap = posRateRawByLap;
PCSimple.ByLap.rawOccuByLap = rawOccuByLap;
PCSimple.ByLap.OccuByLap = OccuByLap;
PCSimple.badBins = badBins;
bins = 1:100;
PCSimple.goodBins = bins(~ismember(bins, badBins));
goodBins = PCSimple.goodBins;

pfm = mean(PCSimple.Shuff.isPC);
occu = nansum(PCSimple.Occupancy);

%%%%Info (bits/spike and bits/second) calcuations:
nCells = size(PCSimple.posRates, 1);
PCSimple.InfoPerSpk = NaN(nCells, 1);
PCSimple.InfoPerSec = NaN(nCells, 1);
PCSimple.Shuff.InfoPerSpk = NaN(nCells, size(ratesAllShuff, 3));
PCSimple.Shuff.InfoPerSec = NaN(nCells, size(ratesAllShuff, 3));
PCSimple.Shuff.InfoPerSpkZ = NaN(nCells, 1);
PCSimple.Shuff.InfoPerSecZ = NaN(nCells, 1);
PCSimple.Shuff.InfoPerSpkP = NaN(nCells, 1);
PCSimple.Shuff.InfoPerSecP = NaN(nCells, 1);
for i = 1:nCells
    [infoSp, infoSec] = makeInfoPerSpikeAndSec(PCSimple.posRates(i, goodBins), PCSimple.Occupancy(goodBins));
    PCSimple.InfoPerSpk(i) =  infoSp;
    PCSimple.InfoPerSec(i) = infoSec;
    for sh = 1:size(ratesAllShuff, 3)
        [infoSp, infoSec] = makeInfoPerSpikeAndSec(ratesAllShuff(i, goodBins, sh), PCSimple.Occupancy(goodBins));
        PCSimple.Shuff.InfoPerSpk(i, sh) =  infoSp;
        PCSimple.Shuff.InfoPerSec(i, sh) = infoSec;
    end
    PCSimple.Shuff.InfoPerSpkZ(i) = (PCSimple.InfoPerSpk(i) - nanmean(PCSimple.Shuff.InfoPerSpk(i, :)))./nanstd(PCSimple.Shuff.InfoPerSpk(i, :));
    %Z score
    PCSimple.Shuff.InfoPerSecZ(i) = (PCSimple.InfoPerSec(i) - nanmean(PCSimple.Shuff.InfoPerSec(i, :)))./nanstd(PCSimple.Shuff.InfoPerSec(i, :));
    PCSimple.Shuff.InfoPerSpkP(i) = nanmean(PCSimple.Shuff.InfoPerSpk(i, :) >= PCSimple.InfoPerSpk(i));
    %%P value
    PCSimple.Shuff.InfoPerSecP(i) = nanmean(PCSimple.Shuff.InfoPerSec(i, :) >= PCSimple.InfoPerSec(i));
    
end
%%%%



pfsByLapCell = {};
pfsByLapMat = [];
pfsByLapMat = nan([length(PCSimple.Shuff.isPC), 2]);
isPC = find(PCSimple.Shuff.isPC);
for i = 1:length(isPC)
    lapRate = squeeze(PCSimple.ByLap.posRateByLap(isPC(i), :, :))';
    if size(lapRate, 2) == 1
        lapRate = lapRate';
    end
    lapRate = lapRate./repmat(nanmean(lapRate, 2), 1, size(lapRate, 2));
    pf = PCSimple.Shuff.PFInPos{isPC(i)};
    m1 = nanmean(lapRate(:, pf), 2);
    pfsByLapCell{isPC(i)} = m1;
    pfsByLapMat(isPC(i), :) = [mean(m1 > 1), sum(m1 > 1)];
end


PCLappedSess = PCSimple;
nLap = size(PCSimple.ByLap.OccuByLap, 1);
minL = round(nLap*lapPerc);
if minL < lapMin
    minL = lapMin;
end
PCLappedSess.minLapsN = minL;
pos = 1:100;
isPC = find(PCSimple.Shuff.isPC);
PCLappedSess.Shuff.isPC(:) = 0;
PCLappedSess.Shuff.PFPeakPos(:) = NaN;
PCLappedSess.Shuff.BestField = NaN(size(PCLappedSess.posRates, 1), 1);
PCLappedSess.Shuff.PeakRate(:) = NaN;
for i = 1:size(PCLappedSess.posRates, 1)
    PCLappedSess.Shuff.PFInAllPos{i} = {};
    PCLappedSess.Shuff.LapRelRate{i} = {};
    PCLappedSess.Shuff.PFInPos{i} = [];
end
for i = 1:length(isPC)
    lapRate = squeeze(PCSimple.ByLap.posRateRawByLap(isPC(i), :, :))';
    if size(lapRate, 2) == 1
        lapRate = lapRate';
    end
    lapRate = lapRate./repmat(nanmean(lapRate, 2), 1, size(lapRate, 2));
    kFields = [];
    for ii = 1:length(PCSimple.Shuff.PFInAllPos{isPC(i)})
        pf = PCSimple.Shuff.PFInAllPos{isPC(i)}{ii};
        m1 = nanmean(lapRate(:, pf), 2);
        PCLappedSess.Shuff.LapRelRate{i}{end + 1} = m1;
        if sum(m1 > 1) >= minL
            kFields = [kFields; 1];
        else
            kFields = [kFields; 0];
        end
    end
    if sum(kFields > 0)
        PCLappedSess.Shuff.isPC(isPC(i)) = 1;
        PCLappedSess.Shuff.PFInAllPos{isPC(i)} = PCSimple.Shuff.PFInAllPos{isPC(i)}(kFields > 0);
        PCLappedSess.Shuff.LapRelRate{i} = PCLappedSess.Shuff.LapRelRate{i}(kFields > 0);
        fieldRates = [];
        peakPosAll = [];
        for ii = 1:length(PCLappedSess.Shuff.PFInAllPos{isPC(i)})
            pf = PCLappedSess.Shuff.PFInAllPos{isPC(i)}{ii};
            posIn = pos(pf);
            ratesIn = PCLappedSess.posRates(isPC(i), pf);
            [peakRate, peakPos] = max(ratesIn);
            peakPos = posIn(peakPos);
            fieldRates = [fieldRates; peakRate];
            peakPosAll = [peakPosAll; peakPos];
        end
        [~, bestField] = max(fieldRates);
        PCLappedSess.Shuff.BestField(isPC(i)) = bestField;
        PCLappedSess.Shuff.PeakRate(isPC(i)) = fieldRates(bestField);
        PCLappedSess.Shuff.PFPeakPos(isPC(i)) = peakPosAll(bestField);
        PCLappedSess.Shuff.PFInPos{isPC(i)} = PCLappedSess.Shuff.PFInAllPos{isPC(i)}{bestField};
    end
end
