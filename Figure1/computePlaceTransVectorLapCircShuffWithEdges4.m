function [out, ratesAllShuff] = computePlaceTransVectorLapCircShuffWithEdges4(activity, treadPos, T, lapVec, shuffN, varargin)
%function out = computePlaceTransVectorLapCircShuffWithEdges4(activity, treadPos, T, lapVec, shuffN, excludeVec(opt))

excludeVec = [];
if ~isempty(varargin)
    excludeVec = varargin{1};
    if length(varargin)==2
        minVel = varargin{2};
    else
        minVel = 5;
    end
else
    minVel = 5;
end


if size(activity, 1) > size(activity, 2)
    activity = activity';
end
%randomizes the random number generator
rng('shuffle');

%set parameters
RatePerc = 99; %perc of shuffles that cell's pc has to be above 
edgeRateMultiple = 2; %doesn't matter a ton, additional criteria for edges of place field
trimRunStarts = 0.25; %to eliminate possible artifacts during run starts
trimRunEnds = 0.25; %at the beginning and end of a run epoch, can set to 0
minRunTime = 2; %min duration of a running bout
minPFBins = 5; %min # of spatial bins a plce field needs to be above its shuffle value for it to count as a place cell
%minVel = 5; %if -5 calculate vel<5 cm/s (now put this in as a varargin or
%is set automatically above)

out = [];
out.Params.RatePerc = 99;
out.Params.edgeRateMultiple = edgeRateMultiple;
out.Params.trimRunStarts = trimRunStarts;
out.Params.trimRunEnds = trimRunEnds;
out.Params.minRunTime = minRunTime;
out.Params.minPFBins = minPFBins;
out.Params.minVel = minVel;
out.Params.shuffN = shuffN;


pos = treadPos;

dT = median(diff(T));

g = fspecial('Gaussian',[10, 1], 5);
h1 = linspace(0, 1, 101);
h1(end) = 1.0001;

%determine the valid run epochs (and trim them)
if minVel > 0
    runTimes = calcMovEpochs1(treadPos, T, minVel);
    runTimes(:, 2) = runTimes(:, 2) + 0.00001; %b/o histc works it will exclude the last sample unless you add this
    runTimes(:, 1) = runTimes(:, 1) + trimRunStarts;
    runTimes(:, 2) = runTimes(:, 2) - trimRunEnds;
    runTimes = runTimes(diff(runTimes, [], 2) >= minRunTime, :);
else
    if minVel == 0
        runTimes = [min(T), max(T)];
    end
    if minVel < 0
        runTimes = calcStillEpochs1(treadPos, T, -1*minVel, minRunTime);
        runTimes(:, 2) = runTimes(:, 2) + 0.00001;
        runTimes(:, 1) = runTimes(:, 1) + trimRunStarts;
        runTimes(:, 2) = runTimes(:, 2) - trimRunEnds;
        runTimes = runTimes(diff(runTimes, [], 2) >= minRunTime, :);
    end
end


kRun = inInterval(runTimes, T);
kRun(~(lapVec > 0)) = 0;
if ~isempty(excludeVec)
    kRun(excludeVec > 0) = 0;
end

out.runTimes = runTimes;

%set non-valid time points to NaN
pos(~(kRun > 0)) = NaN;
%kRun are the samples are the valid running times
actNaN = mean(isnan(activity), 1) == 1;
pos(actNaN) = NaN;

whichLap = lapVec(~isnan(pos));
out.nanedPos = pos;

%keep only the positions and activity correspending to valid samples
activity = activity(:, ~isnan(pos));
pos = pos(~isnan(pos));

%calculate the occupancy
[occ1, whichPlace] = histc(pos, h1);
occ1 = occ1(1:(end - 1));
occ1 = occ1*dT;

out.rawOccupancy = occ1;
rawSums = zeros(size(activity, 1), 100);

%calculate the number of events at each spatial bin for each cell
for i = 1:100
    k = whichPlace == i;
    if sum(k) > 0
        rawSums(:, i) = nansum(activity(:, k), 2);
    end
end
out.posSums = rawSums;

%divide the summed by the occupancy to get rates by position
rawPosRates = rawSums./repmat(occ1, size(rawSums, 1), 1);
out.rawPosRates = rawPosRates;

rawPosRates(isnan(rawPosRates)) = 0;

%smooth the occupancy and the rates
out.Occupancy = convolve2(occ1', g, 'wrap')';
out.posRates = convolve2(rawPosRates', g, 'wrap')';

%for each lap circularly shuffle all the activity associated wt that lap
%shuffle the position vector within each lap for valid positions, keep activity
%matrix  the same, compare one wt other, what is my expectation of the
%cell's activity based on firing rate, what's expectation based on
%occupancy
if shuffN > 1
    rng('shuffle');
    out.Shuff.RatePerc = RatePerc;
    posReal = pos;
    ratesAllShuff = zeros(size(activity, 1), 100, shuffN);
    rawRatesAllShuff = zeros(size(activity, 1), 100, shuffN);
    rawSumsAll = rawRatesAllShuff;
    out.Shuff.whichLap = int16(whichLap);
    uLaps = unique(whichLap);
    circShuffLapN = zeros(length(uLaps), shuffN);
    for i = 1:length(uLaps)
        circShuffLapN(i, :)  = randi(sum(whichLap == uLaps(i)), [1, shuffN]);
    end
    out.Shuff.circShuffLapN = circShuffLapN;
    for sh = 1:shuffN
        pos = posReal;
        %shuffle the position vector
        for i = 1:length(uLaps)
            pos(whichLap == uLaps(i)) = circshift(pos(whichLap == uLaps(i)), [0, circShuffLapN(i, sh)]);
        end
        
        [~, whichPlace] = histc(pos, h1);
        
        rawSums = zeros(size(activity, 1), 100);
        
        for i = 1:100
            k = whichPlace == i;
            if sum(k) > 0
                rawSums(:, i) = nansum(activity(:, k), 2);
            end
        end
        rawSumsAll(:, :, sh) = rawSums;
        
        rawPosRates = rawSums./repmat(out.rawOccupancy, size(rawSums, 1), 1);
        rawRatesAllShuff(:, :, sh) = rawPosRates;
        
        rawPosRates(isnan(rawPosRates)) = 0;
        
        posRates = convolve2(rawPosRates', g, 'wrap')';
        %these are all the shuffled rates by position for each cell for each
        %of n shuffles
        ratesAllShuff(:, :, sh) = posRates;
    end
    
    ratesR = out.posRates;
    %these concatanations are to account for the stuff around the edges by
    %putting them back to back can get continous set of bins
    ratesR = [ratesR, ratesR];
    pos = 1:100;
    pos = [pos, pos];
    
    %perc of shuffles to create the threshold
    percRate = prctile(ratesAllShuff, RatePerc, 3);
    out.Shuff.shuffMeanRate = nanmean(ratesAllShuff, 3);
    out.Shuff.ThreshRate = percRate;
    percRate = [percRate, percRate];
    
    sigRate = double(ratesR > percRate);
       %PFinPos give you all of the bins that are within PF, pstn in PF,
       %can get the width of place field from this calc
       %PFPEakPos if its a PC where is its peak? 
    
    pC = find(sum(sigRate, 2) >= minPFBins);
    out.Shuff.isPC = zeros(size(ratesR, 1), 1);
    out.Shuff.PeakRate = NaN(size(ratesR, 1), 1);
    out.Shuff.PFInPos = {};
    out.Shuff.PFInAllPos = {};
    out.Shuff.PFPeakPos = NaN(size(ratesR, 1), 1);
    for i = 1:size(out.posRates, 1)
        out.Shuff.PFInAllPos{i} = {};
        out.Shuff.PFInPos{i} = [];
    end
    
    %next, it iterates over cells to determine whether they are place cells 
    %a lot of the logic is to get around the circular nature of the data
    %a) it finds within the concatanated set, all of the bins that are above the 
    %threshold and asks are any of these last for at least 5 bins (minPFBins), 
    %b) if you have more than one place field don't count it twice, start
    % from largest PF.
    %ST- not sure which line corresponds to a or b
 
    for i = 1:length(pC)
        [sigBins, sigBinsL] = suprathresh(sigRate(pC(i), :), 0.5);
        sigBins = sigBins(sigBinsL >= minPFBins, :);
        sigBinsAll = [];        
        maxR = [];
        
        if ~isempty(sigBins)
            sigBinsAll = [];
            for ii = 1:size(sigBins, 1)
                sigBinsAll = [sigBinsAll, sigBins(ii, 1):sigBins(ii, 2)];
            end
            
            sBins = ratesR(pC(i), :) > nanmean(ratesR(pC(i), :))*edgeRateMultiple;
            sBins = sBins | sigRate(pC(i), :);
            [pfS, pfL] = suprathresh(double(sBins), 0.5);
            pfS = pfS(pfL >= minPFBins, :);
            pfL = pfL(pfL >= minPFBins, :);
            
            k = [];
            for ii = 1:size(pfS, 1)
                if sum(ismember(sigBinsAll, pfS(ii, 1):pfS(ii, 2))) >= minPFBins
                    k = [k; 1];
                else
                    k = [k; 0];
                end
            end
            pfS = pfS(k > 0, :);
            pfL = pfL(k > 0);
            
            k = sum(pfS > 100, 2) < 2;
            pfS = pfS(k, :);
            pfL = pfL(k, :);
            
            [~, s] = sort(pfL, 'descend');
            pfL = pfL(s, :);
            pfS = pfS(s, :);
            
            inPFAll = zeros(1, 100);
            for ii = 1:length(pfL)
                inPFIn = zeros(1, 100);
                if pfS(ii, 2) > 100
                    inPFIn(1:mod(pfS(ii, 2), 100)) = 1;
                    inPFIn(pfS(ii, 1):end) = 1;
                else
                    inPFIn(pfS(ii, 1):pfS(ii, 2)) =1;
                end
                if sum(inPFAll(inPFIn > 0)) == 0
                    inPFAll(inPFIn > 0) = 1;
                    out.Shuff.PFInAllPos{pC(i)}{end + 1} = find(inPFIn > 0);
                    maxR = [maxR; max(ratesR(pC(i), out.Shuff.PFInAllPos{pC(i)}{end}))];
                end
            end
            out.Shuff.isPC(pC(i)) = 1;
            [out.Shuff.PeakRate(pC(i)), keepF] = max(maxR);
            out.Shuff.PFInPos{pC(i)} = out.Shuff.PFInAllPos{pC(i)}{keepF};
            posIn = pos(out.Shuff.PFInAllPos{pC(i)}{keepF});
            ratesIn = ratesR(pC(i), out.Shuff.PFInAllPos{pC(i)}{keepF});
            
            [~, maxBin] = max(ratesIn);
            out.Shuff.PFPeakPos(pC(i)) = posIn(maxBin);
        end
    end
end
