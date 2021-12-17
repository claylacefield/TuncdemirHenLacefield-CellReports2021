function [posBinFrac, posInfo, pcRatesBlanked, pcOmitRatesBlanked, pfOnlyRates, pfOnlyRatesOmit, pcRatesSorted, omitRatesSorted] = cuePosInhib(cueShiftStruc, goodSeg, refLapType, toPlot);

% formerly popPosInfo

%refLapType = 2; 
posRates = cueShiftStruc.PCLappedSessCell{refLapType}.posRates;
if goodSeg==0
pc = find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1);
elseif goodSeg==1
    pc = 1:length(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC);
else
    pc = goodSeg;
end
pcRates = posRates(pc,:);

[maxs, inds] = max(pcRates'); % find bin of peak firing rate for PCs
[maxs, inds2] = max(posRates(find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1),:)');

numBins = 20; %10;
[counts, edges, binInd] = histcounts(inds2, numBins);

if toPlot
figure; 
subplot(2,1,1);
bar(counts);
title('numUnits tuned to pos');
end

posBinFrac = counts/length(inds2);

ips = cueShiftStruc.PCLappedSessCell{refLapType}.InfoPerSec(find(cueShiftStruc.PCLappedSessCell{refLapType}.Shuff.isPC==1)); % for goodSeg

% avg infoPerSpk/Sec for units with similar tuning
for i = 1:numBins
    unitInds = find(binInd==i);
    posInfo(i) = mean(ips(unitInds));
end

if toPlot
subplot(2,1,2); 
bar(posInfo);
title('pos mean info per spk');
end

%% look at cue inhibition

[sorted, sortInds] = sort(inds);
pcRatesSorted = pcRates(sortInds,:);
%figure; imagesc(pcRatesSorted);
posRates2 = cueShiftStruc.PCLappedSessCell{end}.posRates; % for Omit laps
pcRates2 = posRates2(pc,:);
omitRatesSorted = pcRates2(sortInds,:);

% go through posRates, and blank out time around peak
j=0;
for i = 1:size(pcRatesSorted,1) % for all cells (actually not just pc now)
    rates = pcRatesSorted(i,:);
    rates2 = omitRatesSorted(i,:);
    pfRates = pcRatesSorted(i,:);
    pfRates2 = omitRatesSorted(i,:);
    pkPos = sorted(i); % pos bin of place field peak
    if pkPos <= 10 % wraparound for units at beginning
        rates(1:pkPos+9) = NaN; % rates(1:pkPos+10) = NaN;
        rates(100-(9-pkPos):100) = NaN; % rates(100-10-pkPos:100) = NaN;
        rates2(1:pkPos+9) = NaN; % rates2(1:pkPos+10) = NaN;
        rates2(100-(9-pkPos):100) = NaN; % rates2(100-10-pkPos:100) = NaN;
        
        % place field rates only 
        pfRates(pkPos+10:100-(10-pkPos)) = NaN;
        pfRates2(pkPos+10:100-(10-pkPos)) = NaN;
        
    elseif pkPos>90 % wraparound for units at end
        rates(pkPos-9:100) = NaN; % changing all 10 to 9 for pfBlanked
        rates(1:9-(100-pkPos)) = NaN;
        rates2(pkPos-9:100) = NaN;
        rates2(1:9-(100-pkPos)) = NaN;
        
        % place field rates only 
        pfRates(10-(100-pkPos):pkPos-10) = NaN;
        pfRates2(10-(100-pkPos):pkPos-10) = NaN;
    else % units in middle
        rates(pkPos-9:pkPos+9) = NaN;
        rates2(pkPos-9:pkPos+9) = NaN;
        
        % place field rates only 
        pfRates(1:pkPos-10) = NaN;
        pfRates(pkPos+10:end) = NaN;
        pfRates2(1:pkPos-10) = NaN;
        pfRates2(pkPos+10:end) = NaN;
        
        if pkPos>45 && pkPos<55 % for potential middle-cue cells (PF in pos 40-60)
            j=j+1;
            cuePosRate(j) = mean(pcRatesSorted(i,45:55),2); % mean rate in cue trials
            omitPosRate(j) = mean(omitRatesSorted(i,45:55),2);  % and omit trials
            cuePkPos(j) = pkPos;
        end
        
    end
    pcRatesBlanked(i,:) = rates;
    pcOmitRatesBlanked(i,:) = rates2;
    pfOnlyRates(i,:) = pfRates;
    pfOnlyRatesOmit(i,:) = pfRates2;   
end

% cuePosRate = cuePosRate(cuePosRate~=0); % just elim zeros
% omitPosRate = omitPosRate(omitPosRate~=0);

if toPlot
figure; 
colormap(jet); 
subplot(2,2,1);
imagesc(pcRatesBlanked);
title('pcRatesBlanked');
subplot(2,2,2);
imagesc(pcOmitRatesBlanked);
title('pcOmitRatesBlanked');
subplot(2,2,3);
plot(nanmean(pcRatesBlanked,1));
hold on;
plot(nanmean(pcOmitRatesBlanked,1),'r');
legend('refLaps', 'omitLaps');
subplot(2,2,4);
%plot(cuePosRate,omitPosRate,'x');
plot(cuePkPos,cuePosRate-omitPosRate,'x');
%hold on; line([0 0.2], [0 0.2]);
title('middle cell mean ref lap rate - omit');
end

nonEp1 = 11:50; %11:44;
nonEp2 = 61:89; %56:89;
startEp1 = 1:10; 
startEp2 = 90:100;
midEp = 51:60; %45:55;

if toPlot
figure;
%bar([mean(mean(pcRates(:,40:60),2),1) mean(mean(pcRates2(:,40:60),2),1)]);
bar([mean([nanmean(mean(pcRatesBlanked(:,nonEp1),2),1) nanmean(mean(pcRatesBlanked(:,nonEp2),2),1)]) mean([nanmean(mean(pcRatesBlanked(:,startEp1),2),1) nanmean(mean(pcRatesBlanked(:,startEp2),2),1)]) nanmean(mean(pcRatesBlanked(:,midEp),2),1) nanmean(mean(pcOmitRatesBlanked(:,midEp),2),1)]);
title('pkPos blanked non-cue, startCue, middleCue, omitCue');
%legend('startCueBlanked', 'middleCueBlanked', 'middleOmitBlanked');
ylabel('mean rate');

figure; 
subplot(2,1,2);
plot(nanmean(pfOnlyRates,1));
subplot(2,1,1);
imagesc(pfOnlyRates);
title('pfOnly rates');

end
