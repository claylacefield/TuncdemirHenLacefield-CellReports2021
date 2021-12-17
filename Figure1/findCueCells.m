function [cueCellStruc] = findCueCells(cueShiftStruc, eventName, segDictCode, toPlot)

%% USAGE: findCueCells(cueShiftStruc);
% This function finds the identity of putative start and middle/variable
% cue cells based upon:
% 1.) is a place cell (i.e. fires significantly more at cue location than
% elsewhere)
% 2.) position: firing within certain range of the cue
% 3.) (for middle cue) omitCue responses: if responses go away when cue not
% present.

% Issues:
%
% 073119:
%   - exactly what range do we use for cue location (start and middle)
%       - static location, or based upon TDML cue location (pin or
%       context?)
%   - omitCue response: should I take from avgCueTrigSig? or posRates at
%   middle cue location?

if segDictCode~=0
    if isnumeric(segDictCode)
    segDictName = findLatestFilename('segDict', 'goodSeg');
    else
        segDictName = segDictCode;
    end
else
    segDictName = uigetfile('*.mat', 'Select segDict/seg2P file');
end


%% Find "normal" lap type (middle cue in typical location)
lapTypeArr = cueShiftStruc.lapCueStruc.lapTypeArr;
lapTypeArr(lapTypeArr==0) = max(lapTypeArr)+1; % make omitCue last
for i=1:length(cueShiftStruc.pksCellCell)
    numLapType(i) = length(find(lapTypeArr==i));
end
[val, refLapType] = max(numLapType); % use ref lap from one with most laps 

% extract posRates for normal/reference laps
posRatesRef = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;

%% Criteria 1: place cell
% take PCs from normal laps
pc = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1);

%% Criteria 2: position
% find the peak firing position of each cell
% (or should we use centerOfMass?)
[maxVal, maxInd] = max(posRatesRef');
pcPkPos = maxInd(pc); % just PCs

%% Start cue cells
% find start cue cells just based upon position (ind w.re. to all seg)
startCueCellInd = pc(find(pcPkPos>=90 | pcPkPos<=10));

% extract posRates for omit laps (all cells)
posRatesOmit = cueShiftStruc.PCLappedSessCell{end}.posRates;
%posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;

% plot to check
if toPlot
figure('Position',[0,50,800,800]);

subplot(2,2,1);
[sortInd] = plotUnitsByTuning(posRatesRef(startCueCellInd,:), 0, 1);
cl = caxis;
title('start cue cells');

% plot to check
subplot(2,2,3); 
colormap(jet); imagesc(posRatesOmit(startCueCellInd(sortInd),:));  caxis(cl);
title('omit laps');

subplot(2,2,2); 
%figure;
plot(mean(posRatesRef(startCueCellInd,:),1), 'b');
hold on; 
plot(mean(posRatesOmit(startCueCellInd,:),1), 'r');
title('avgs');
xlabel('pos');
ylabel('mean rate (Hz)');
legend('cue laps', 'omit laps');
end

%% Middle PC ind
midCellInd = pc(find(pcPkPos>=40 & pcPkPos<=65)); % pc(find(pcPkPos>=45 & pcPkPos<=65));
% NOTE that 45 might not include some predictive cells that fire just
% before the cue (but may not omit/shift)

% plot to check
if toPlot
%maxVal = max(max(posRatesRef(midCellInd,:)));
figure('Position',[0,50,800,800]);
subplot(2,2,1);
[sortInd] = plotUnitsByTuning(posRatesRef(midCellInd,:), 0, 1);
cl = caxis;
title('midCellInd');

% extract posRates for omit laps
%posRatesOmit = cueShiftStruc.PCLappedSessCell{end}.posRates;

% plot to check
subplot(2,2,3); 
colormap(jet); imagesc(posRatesOmit(midCellInd(sortInd),:));  caxis(cl);
title('omit laps');

subplot(2,2,2); 
%figure;
plot(mean(posRatesRef(midCellInd,:),1), 'b');
hold on; 
plot(mean(posRatesOmit(midCellInd,:),1), 'r');
title('avgs');
xlabel('pos');
ylabel('mean rate (Hz)');
legend('cue laps', 'omit laps');
end

%% Criteria 3: lack of omitCue response
try  % TRY/CATCH all sections based on omit laps, just in case they don't exist
% find mean firing rate for area around peak (+/-5)
for i=1:length(midCellInd)
    midFieldRateRef(i) = mean(posRatesRef(midCellInd(i),maxInd(midCellInd(i))-5:maxInd(midCellInd(i))+5),2);
    midFieldRateOmit(i) = mean(posRatesOmit(midCellInd(i),maxInd(midCellInd(i))-5:maxInd(midCellInd(i))+5),2);
end

%% Method #1: "cue cells" are middle cells with >= 2x posRate for cue trials vs. omit
% quick and dirty for now: midCueCellInd = twice as big response to cue as omit
midCueCellInd = midCellInd(find(midFieldRateRef./midFieldRateOmit>2));
% NOTE: some are Inf, which is still>2, but these points don't plot below

% take out Inf's for ratios (max at 25)
refOmitRatio = midFieldRateRef./midFieldRateOmit;
for i=1:length(refOmitRatio)
    if refOmitRatio(i)==Inf
        refOmitRatio(i)=25;
    end
end

nonCueCellInd = setdiff(midCellInd, midCueCellInd); % nonCueCells are middle cells that still have omit activ

% plots
if toPlot
    try
    figure;
    subplot(2,2,1);
    plot(maxInd(midCellInd),refOmitRatio,'.');
    hold on;
    line([45 65], [1 1]);
    title('spatial profile of cue/omit');
    xlabel('pos');
    
    %figure;
    subplot(2,2,2);
    plot(midFieldRateRef,midFieldRateOmit,'.');
    hold on;
    line([0 0.5], [0 0.5]);
    title('cue vs. omit');
    xlabel('cue'); ylabel('omit');
    
    %figure;
    subplot(2,2,3);
    hist(refOmitRatio,20);
    xlabel('midFieldRateRef./midFieldRateOmit');
    
    figure('Position', [50,100,800,800]);
    subplot(2,2,1);
    [sortInd] = plotUnitsByTuning(posRatesRef(midCueCellInd,:), 0, 1);
    cl = caxis;
    title('midCueCell (2x PF rate) cueLaps');
    subplot(2,2,3);
    colormap(jet); imagesc(posRatesOmit(midCueCellInd(sortInd),:)); caxis(cl);
    title('midCueCell omitLaps');
    if length(numLapType)==3
        posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;
        subplot(2,2,4);
        colormap(jet); imagesc(posRatesShift(midCueCellInd(sortInd),:)); caxis(cl);
        title('midCueCell ShiftLaps');
    end
    
    subplot(2,2,2);
    plot(mean(posRatesRef(midCueCellInd,:),1), 'b');
    hold on;
    plot(mean(posRatesOmit(midCueCellInd,:),1), 'r');
    title('avgs');
    xlabel('pos');
    ylabel('mean rate (Hz)');
    legend('cue laps', 'omit laps');
    
    figure('Position', [50,100,800,800]);
    subplot(2,2,1);
    [sortInd] = plotUnitsByTuning(posRatesRef(nonCueCellInd,:), 0, 1);
    cl = caxis;
    title('nonCueCell (<2x PF rate) cueLaps');
    subplot(2,2,3);
    colormap(jet); imagesc(posRatesOmit(nonCueCellInd(sortInd),:)); caxis(cl);
    title('nonCueCell omitLaps');
    if length(numLapType)==3
        %posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;
        subplot(2,2,4);
        colormap(jet); imagesc(posRatesShift(nonCueCellInd(sortInd),:)); caxis(cl);
        title('nonCueCell ShiftLaps');
    end
    
    subplot(2,2,2);
    plot(mean(posRatesRef(nonCueCellInd,:),1), 'b');
    hold on;
    plot(mean(posRatesOmit(nonCueCellInd,:),1), 'r');
    title('avgs');
    xlabel('pos');
    ylabel('mean rate (Hz)');
    legend('cue laps', 'omit laps');
    catch
    end
end

%% scratch

% playing with shuffle significance test for cueVsOmit response
% NOTE: may also t-test on amplitude
% pseudo:
% 1.) find cue-triggered responses for each unit, find peak/auc amplitude
% of each response
% 2.) calc diff in mean response
% 3.) put together cue and omit and randomly resample cue and omit laps
% 4.) calc diff in mean response

%% Method #2: "cue cells" are middle cells w. ttest2<0.05 for event amp
% This is based upon ttest2 bet cue and omit (NOTE: event max amp)
%eventName = 'tact';

% if segDictCode~=0
%     if isnumeric(segDictCode)
%     segDictName = findLatestFilename('segDict', 'goodSeg');
%     else
%         segDictName = segDictCode;
%     end
% else
%     segDictName = uigetfile('*.mat', 'Select segDict/seg2P file');
% end

midCueCellInd2 = [];
for i = 1:length(midCellInd)
    try
    % calc avgCueTrigSig for each middle cell
%     if segDictCode~=0
        [cueTrigSigStruc] = avgCueTrigSig(midCellInd(i), eventName, 0, segDictName);
%     else
%         try
%             [cueTrigSigStruc] = avgCueTrigSig(midCellInd(i), eventName, 0, cueTrigSigStruc.segDictName);
%         catch
%             [cueTrigSigStruc] = avgCueTrigSig(midCellInd(i), eventName, 0);
%         end
%     end
    
    omitCueSig = cueTrigSigStruc.omitCueSig;
    midCueSig = cueTrigSigStruc.midCueSig;
    
    % find max event amplitude following cue (minus baseline) for norm laps
    for j=1:size(midCueSig,2)
        midCueAmp(j) = max(midCueSig(30:130,j)-midCueSig(30,j));
    end
    
    % and omit laps
    for j=1:size(omitCueSig,2)
        omitCueAmp(j) = max(omitCueSig(30:130,j)-omitCueSig(30,j));
    end
    
    % ttest2 on event amplitudes
    [h,p,ci,stats] = ttest2(midCueAmp, omitCueAmp);
    
    % if cue event amplitudes signif > omit, then add cell to list
    if p<=0.05 && mean(midCueAmp)>mean(omitCueAmp)%h==1
        midCueCellInd2 = [midCueCellInd2 midCellInd(i)];
    end
    catch
    end
end

% plot
if toPlot
    try
        figure('Position', [100,150,800,800]);
        subplot(2,2,1);
        [sortInd] = plotUnitsByTuning(posRatesRef(midCueCellInd2,:), 0, 1);
        cl = caxis;
        title('midCueCell (ttest) cueLaps');
        subplot(2,2,3);
        colormap(jet); imagesc(posRatesOmit(midCueCellInd2(sortInd),:)); caxis(cl);
        title('midCueCell omitLaps');
        if length(numLapType)==3
            %posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;
            subplot(2,2,4);
            colormap(jet); imagesc(posRatesShift(midCueCellInd2(sortInd),:)); caxis(cl);
            title('midCueCell ShiftLaps');
        end
        
        subplot(2,2,2);
        plot(mean(posRatesRef(midCueCellInd2,:),1), 'b');
        hold on;
        plot(mean(posRatesOmit(midCueCellInd2,:),1), 'r');
        title('avgs');
        xlabel('pos');
        ylabel('mean rate (Hz)');
        legend('cue laps', 'omit laps');
    catch
    end
end

%% Method #3: "cue cells" have posRate diff cue-omit > 95% lap shuffled
% now trying shuffle with laps
midCellRateDiff = midFieldRateRef-midFieldRateOmit; % observed diff in middle rate for every middle cell
posRateLapRef = cueShiftStruc.PCLappedSessCell{refLapType}.ByLap.posRateByLap;
posRateLapOmit = cueShiftStruc.PCLappedSessCell{end}.ByLap.posRateByLap;
numRefLaps = size(posRateLapRef,3);
numOmitLaps = size(posRateLapOmit,3);

posRateLapRefOmit = cat(3, posRateLapRef, posRateLapOmit);
    
midCueCellInd3 = [];
for i = 1:length(midCellInd)
    try
    for j=1:100
        
        % resample laps
        refLapRes = randsample(numRefLaps+numOmitLaps, numRefLaps);
        omitLapRes = setdiff(1:(numRefLaps+numOmitLaps),refLapRes);
        
        posRateRefRes = squeeze(mean(posRateLapRefOmit(midCellInd(i),:,refLapRes),3));
        posRateOmitRes = squeeze(mean(posRateLapRefOmit(midCellInd(i),:,omitLapRes),3));
        
        midFieldRateRefRes = mean(posRateRefRes(maxInd(midCellInd(i))-5:maxInd(midCellInd(i))+5));
        midFieldRateOmitRes = mean(posRateOmitRes(maxInd(midCellInd(i))-5:maxInd(midCellInd(i))+5));
        resampMeanDiff(j) = midFieldRateRefRes-midFieldRateOmitRes;
    end
    
    numGreater = length(find(resampMeanDiff>=midCellRateDiff(i)));
    
    if numGreater<=5
        midCueCellInd3 = [midCueCellInd3 midCellInd(i)];
    end
    catch
    end
    
end

if toPlot
    try
        figure('Position', [150,200,800,800]);
        subplot(2,2,1);
        [sortInd] = plotUnitsByTuning(posRatesRef(midCueCellInd3,:), 0, 1);
        cl = caxis;
        title('midCueCell (shuff) cueLaps');
        subplot(2,2,3);
        colormap(jet); imagesc(posRatesOmit(midCueCellInd3(sortInd),:)); caxis(cl);
        title('midCueCell omitLaps');
        if length(numLapType)==3
            %posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;
            subplot(2,2,4);
            colormap(jet); imagesc(posRatesShift(midCueCellInd3(sortInd),:)); caxis(cl);
            title('midCueCell ShiftLaps');
        end
        
        subplot(2,2,2);
        plot(mean(posRatesRef(midCueCellInd3,:),1), 'b');
        hold on;
        plot(mean(posRatesOmit(midCueCellInd3,:),1), 'r');
        title('avgs');
        xlabel('pos');
        ylabel('mean rate (Hz)');
        legend('cue laps', 'omit laps');
    catch
    end
end

catch
    disp('No omit laps');
end

%% shift
% testing on 190628_B4, seg 46, 87 (lower on shift)
% a.) calc rate at normal location vs shift location
% OR
% b.) shuffle avgCueTrigSig; [cueTrigSigStruc] = avgCueTrigSig(segNum, eventName, toPlot)
% c.) also look at shift in cue cells by omit, and all midCellInd (i.e. are
% cells that don't show omit more likely to shift?

% posRates shift vs normal
try
    %if length(numLapType)==3 %3
    
    %posRatesShift = cueShiftStruc.PCLappedSessCell{1}.posRates;
    
    lapTypeCuePos = cueShiftStruc.lapCueStruc.lapTypeCuePos;
    
    for i=1:length(midCellInd)
        try
        midFieldRateRef(i) = mean(posRatesRef(midCellInd(i),20:35),2);
        midFieldRateShift(i) = mean(posRatesShift(midCellInd(i),20:35),2);
        catch
        end
    end
    
    % now trying shuffle with laps
    midCellRateDiff = midFieldRateShift-midFieldRateRef;
    posRateLapRef = cueShiftStruc.PCLappedSessCell{refLapType}.ByLap.posRateByLap;
    posRateLapShift = cueShiftStruc.PCLappedSessCell{1}.ByLap.posRateByLap;
    numRefLaps = size(posRateLapRef,3);
    numShiftLaps = size(posRateLapShift,3);
    
    posRateLapRefShift = cat(3, posRateLapRef, posRateLapShift); % concatenate ref and shift lap posRates for resampling
    
    midShiftCellInd3 = [];
    for i = 1:length(midCellInd)
        try
        for j=1:100 % shuffle iterations
            
            % resample laps
            refLapRes = randsample(numRefLaps+numShiftLaps, numRefLaps);
            shiftLapRes = setdiff(1:(numRefLaps+numShiftLaps),refLapRes);
            
            % find rates based upon resampled laps
            posRateRefRes = squeeze(mean(posRateLapRefShift(midCellInd(i),:,refLapRes),3));
            posRateShiftRes = squeeze(mean(posRateLapRefShift(midCellInd(i),:,shiftLapRes),3));
            
            midFieldRateRefRes = mean(posRateRefRes(20:35)); % mean of shift range (~25) for ref laps (resamp)
            midFieldRateShiftRes = mean(posRateShiftRes(20:35)); % mean of shift range (~25) for shift laps (resamp)
            resampMeanDiff(j) = midFieldRateShiftRes-midFieldRateRefRes; % diff bet ref/mid and shift
        end
        
        numGreater = length(find(resampMeanDiff>=midCellRateDiff(i))); % # shuff where shift pos posRates on shift laps>ref laps
        
        if numGreater<=5
            midShiftCellInd3 = [midShiftCellInd3 midCellInd(i)];
        end
        catch
        end
        
    end
    
    %end
    
    if toPlot
        try
            figure('Position', [200,250,800,800]);
            subplot(2,2,1);
            [sortInd] = plotUnitsByTuning(posRatesRef(midShiftCellInd3,:), 0, 1);
            cl = caxis;
            title('midShiftCell (posRates shuff) cueLaps');
            subplot(2,2,3);
            colormap(jet); imagesc(posRatesShift(midShiftCellInd3(sortInd),:)); caxis(cl);
            title('midShiftCell shiftLaps');
            subplot(2,2,4);
            colormap(jet); imagesc(posRatesOmit(midShiftCellInd3(sortInd),:)); caxis(cl);
            title('midShiftCell omitLaps');
            
            subplot(2,2,2);
            plot(mean(posRatesRef(midShiftCellInd3,:),1), 'b');
            hold on;
            plot(mean(posRatesShift(midShiftCellInd3,:),1), 'r');
            title('avgs');
            xlabel('pos');
            ylabel('mean rate (Hz)');
            legend('cue laps', 'shift laps');
        catch
        end
    end
catch
end
% end % end IF numLapType>=3

%% now shuffle on shift amplitudes
try
inds = midCueCellInd; %3; % only look amongst midCueCells
catch
    disp('No cueCells/omitCue, so using all for shift');
    inds = midCellInd;
end
% midCueCellInd3 = posRate diff cue-omit > 95% lap shuffled
% midCueCellInd2 =  middle cells w. ttest2<0.05 for event amp
% midCueCellInd = 2x posRate at cue vs. omit
try
    %if length(numLapType)==3 % if there are 3 lap types (thus shift)
    midShiftCellInd2 = [];
    for i = 1:length(inds)
        try
        % calc avgCueTrigSig for each middle cell
%         try
            toZ = 1;
            [cueTrigSigStruc] = avgCueTrigSig(inds(i), eventName, 0, segDictName, toZ); % cueTrigSigStruc.
%         catch
%             [cueTrigSigStruc] = avgCueTrigSig(inds(i), eventName, 0);
%         end
        
        shiftCueSig = cueTrigSigStruc.shiftCueSig;
        midCueSig = cueTrigSigStruc.midCueSig;
        
        % find max event amplitude following cue (minus baseline) for norm laps
        for j=1:size(midCueSig,2)
            midCueAmp(j) = max(midCueSig(30:130,j)-midCueSig(30,j)); % or sum?
        end
        
        % and Shift laps
        for j=1:size(shiftCueSig,2)
            shiftCueAmp(j) = max(shiftCueSig(30:130,j)-shiftCueSig(30,j)); % or sum?
        end
        
        %         % ttest2 on event amplitudes
        %         [h,p,ci,stats] = ttest2(midCueAmp, shiftCueAmp);
        
        % shuffle
        allAmps = [midCueAmp shiftCueAmp];
        for j = 1:100
            % resample laps
            
            refLapRes = randsample(length(midCueAmp)+length(shiftCueAmp), length(midCueAmp));
            shiftLapRes = setdiff(1:(length(midCueAmp)+length(shiftCueAmp)),refLapRes);
            
            avMidCueAmpRes(j) = mean(allAmps(refLapRes));
            avShiftCueAmpRes(j) = mean(allAmps(shiftLapRes));
            
        end
        
        % if cue event amplitudes signif > omit, then add cell to list
        avMidCueAmp(i) = mean(midCueAmp); avShiftCueAmp(i) = mean(shiftCueAmp);
        if length(find(abs(avMidCueAmpRes-avShiftCueAmpRes)>=abs(avMidCueAmp(i)-avShiftCueAmp(i))))<=5
            midShiftCellInd2 = [midShiftCellInd2 inds(i)];
        end
        catch
        end
        
    end
    
    % plot
    if toPlot
        try
        figure('Position', [100,150,800,800]);
        subplot(2,2,1);
        [sortInd] = plotUnitsByTuning(posRatesRef(midShiftCellInd2,:), 0, 1);
        cl = caxis;
        title('midShiftCell (evAmp shuff) cueLaps');
        subplot(2,2,3);
        colormap(jet); imagesc(posRatesShift(midShiftCellInd2(sortInd),:)); caxis(cl);
        title('midShiftCell ShiftLaps');
        subplot(2,2,4);
        colormap(jet); imagesc(posRatesOmit(midShiftCellInd2(sortInd),:)); caxis(cl);
        title('midShiftCell omitLaps');
        
        subplot(2,2,2);
        plot(mean(posRatesRef(midShiftCellInd2,:),1), 'b');
        hold on;
        plot(mean(posRatesOmit(midShiftCellInd2,:),1), 'r');
        title('avgs');
        xlabel('pos');
        ylabel('mean rate (Hz)');
        legend('cue laps', 'Shift laps');
        catch
        end
    end
catch
    disp('Prob with shift calc #2');
end

%end % end IF numLapType>=3

%% pack some stuff in output structure
cueCellStruc.path = cueShiftStruc.path;
cueCellStruc.cueShiftName = cueShiftStruc.filename;
cueCellStruc.posRatesRef = posRatesRef; % posRates for reference laps
cueCellStruc.startCueCellInd = startCueCellInd; 
cueCellStruc.placeCellInd = setdiff(pc, [midCellInd; startCueCellInd]); 
cueCellStruc.midCellInd = midCellInd; % just cells in middle
try % in case no omit laps
    cueCellStruc.midCueCellInd = midCueCellInd; % cells with 2x PF rate vs omit
    cueCellStruc.nonCueCellInd = nonCueCellInd;
    try % in case nothing significant
        cueCellStruc.midCueCellInd3 = midCueCellInd3;   % cue cells by lap shuffle event amplitude
        cueCellStruc.midCueCellInd2 = midCueCellInd2;   % cue cells by ttest
    catch
    end
catch
end

try % in case no shift laps
    cueCellStruc.midShiftCellInd3 = midShiftCellInd3;
    cueCellStruc.midShiftCellInd2 = midShiftCellInd2;
    try
        cueCellStruc.allMidCueCellInd = unique([midCueCellInd2 midCueCellInd3 midShiftCellInd2 midShiftCellInd3]);
    catch
        cueCellStruc.allMidCueCellInd = unique([midShiftCellInd2 midShiftCellInd3]);
    end
    
    cueCellStruc.avMidCueAmp = avMidCueAmp;
    cueCellStruc.avShiftCueAmp = avShiftCueAmp;
    
catch
    try
    cueCellStruc.allMidCueCellInd = unique([midCueCellInd2 midCueCellInd3]);
    catch
    end
    
    try
    cueCellStruc.avMidCueAmp = avMidCueAmp;
    cueCellStruc.avShiftCueAmp = avShiftCueAmp;
    catch
    end
    
end


if toPlot
    try
        figure('Position', [100,150,800,800]);
        subplot(2,2,1);
        [sortInd] = plotUnitsByTuning(posRatesRef(cueCellStruc.allMidCueCellInd,:), 0, 1);
        cl = caxis;
        title('allMidCueCells (by omit and shift) cueLaps');
        subplot(2,2,3);
        colormap(jet);
        try
            imagesc(posRatesShift(cueCellStruc.allMidCueCellInd(sortInd),:)); caxis(cl);
        catch
        end
        title('allMidCueCells ShiftLaps');
        subplot(2,2,4);
        colormap(jet);
        try
            imagesc(posRatesOmit(cueCellStruc.allMidCueCellInd(sortInd),:)); caxis(cl);
        catch
        end
        title('allMidCueCells omitLaps');
        
        subplot(2,2,2);
        plot(mean(posRatesRef(cueCellStruc.allMidCueCellInd,:),1), 'b');
        hold on;
        try
            plot(mean(posRatesOmit(cueCellStruc.allMidCueCellInd,:),1), 'r');
        catch
        end
        try
            plot(mean(posRatesShift(cueCellStruc.allMidCueCellInd,:),1), 'g');
        catch
        end
        
        title('avgs');
        xlabel('pos');
        ylabel('mean rate (Hz)');
        legend('cue laps', 'Shift laps');
    catch
    end
end
