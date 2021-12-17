function [cueCellInd, relAmp] = cueCompetAnalysis()




% 1. find cue cells (vs. omit) for each cue
% 2. run avgCueTrigSigNew
% 3. compare avg amplitudes

% criteria for cue cells
% - bigger than pre-cue
% - bigger than omit
% - 

% then 2D plot of ratios of single cue to paired

load(findLatestFilename('cueShiftStruc'));
segDictName = findLatestFilename('segDict');
load(segDictName);
[refLapType] = findRefLapType(cueShiftStruc);

cueTypeCell = {'olf', 'led'};
% cueType1 = 'olf';
% cueType2 = 'led';
% cuePos1 = 500;
% cuePos2 = 1150;


%% now shuffle on shift amplitudes

%eventName = cueType1;
cueCellInd = [];
cueCellInd2 = [];
numCell = 0;
for k=1:length(cueTypeCell)
    tic;
    eventName = cueTypeCell{k};
    try

        %cueCellInd = [];
        for i = 1:length(cueShiftStruc.pksCell)
            cueAmp = []; omitCueAmp = [];
            try
                if mod(i,20)==0
                    disp(['Processing unit ' num2str(i) ' of ' num2str(length(cueShiftStruc.pksCell))]);
                end
                
                ca = C(i,:);
                maxCa = max(ca)-min(ca);
                
                % calc avgCueTrigSig for each cell
                [cueTrigSigStruc] = avgCueTrigSigNew(i, eventName, 0, segDictName); % cueTrigSigStruc.
                % NOTE: need to fix this to do omit pos correctly
                
                evTrigSigCell = cueTrigSigStruc.evTrigSigCell;
                refLapEvTrigSig = evTrigSigCell{refLapType}; % refLapType
                omitTrigSig = cueTrigSigStruc.omitTrigSig;
                
                % find max event amplitude following cue (minus baseline) for norm laps
                for j=1:size(refLapEvTrigSig,2)
                    try
                        cueAmp(j) = max(refLapEvTrigSig(30:110,j))-min(refLapEvTrigSig(1:30,j)); % findpeaks(refLapEvTrigSig(30:100,j),'MinPeakDistance', 69)-refLapEvTrigSig(30,j); % max(refLapEvTrigSig(30:130,j)-refLapEvTrigSig(30,j)); % or sum?
                    catch
                        cueAmp(j) = 0;
                    end
                end
                
                % and omit laps
                for j=1:size(omitTrigSig,2)
                    try
                        omitCueAmp(j) = max(omitTrigSig(30:110,j))-min(omitTrigSig(1:30,j)); % findpeaks(omitTrigSig(30:100,j),'MinPeakDistance',69)-omitTrigSig(30,j); %  % or sum?
                    catch
                        omitCueAmp(j) = 0;
                    end
                end
                
                % ttest2 on event amplitudes
                [h,p,ci,stats] = ttest2(cueAmp, omitCueAmp);
                [h2,p2,ci,stats] = ttest2(mean(refLapEvTrigSig(50:110,:),1), mean(refLapEvTrigSig(1:30,:),1));
                
                % having a hard time getting these criteria right
                % Criteria
                % 1.) pretty good pval cue/omit ttest2, 2.) sig pval cue vs pre
                % cue, 3.) cue > omit, 4.) cue > pre cue
                if p<0.14 && h2==1 && mean(cueAmp)>mean(omitCueAmp) && mean(mean(refLapEvTrigSig(50:120,:),1))>mean(mean(refLapEvTrigSig(1:30,:),1)) && max(cueAmp)>maxCa/2
                    numCell = numCell+1;
                    cueCellInd = [cueCellInd i];
                    disp(['Found ' cueTypeCell{k} ' cell #' num2str(i)]);
                    %cueCellAmps(numCell,1) = mean(cueAmp);
                    relAmp(:,numCell) = [max(mean(evTrigSigCell{1},2))/max(mean(evTrigSigCell{2},2)), max(mean(evTrigSigCell{3},2))/max(mean(evTrigSigCell{2},2))];
                end
                
%                         % shuffle
%                         allAmps = [cueAmp omitCueAmp];
%                         for j = 1:100
%                             % resample laps
%                 
%                             refLapRes = randsample(length(cueAmp)+length(omitCueAmp), length(cueAmp));
%                             omitLapRes = setdiff(1:(length(cueAmp)+length(omitCueAmp)),refLapRes);
%                 
%                             avCueAmpRes(j) = mean(allAmps(refLapRes));
%                             avOmitCueAmpRes(j) = mean(allAmps(omitLapRes));
%                 
%                         end
%                 
%                         % if cue event amplitudes signif > omit, then add cell to list
%                         avCueAmp(i) = mean(cueAmp); avOmitCueAmp(i) = mean(omitCueAmp);
%                         if length(find(abs(avCueAmpRes-avOmitCueAmpRes)>=abs(avCueAmp(i)-avOmitCueAmp(i))))<=5
%                             cueCellInd2 = [cueCellInd2 inds(i)];
%                         end
            catch
                disp(['Prob w cell #' num2str(i)]);
            end 
        end
        
    catch
        disp(['Prob w cueType ' num2str(k)]);
    end
    toc;
end


%figure; 
%plot();
    
%     % plot
%     if toPlot
%         try
%         figure('Position', [100,150,800,800]);
%         subplot(2,2,1);
%         [sortInd] = plotUnitsByTuning(posRatesRef(cueCellInd2,:), 0, 1);
%         cl = caxis;
%         title('midShiftCell (evAmp shuff) cueLaps');
%         subplot(2,2,3);
%         colormap(jet); imagesc(posRatesShift(cueCellInd2(sortInd),:)); caxis(cl);
%         title('midShiftCell ShiftLaps');
%         subplot(2,2,4);
%         colormap(jet); imagesc(posRatesOmit(cueCellInd2(sortInd),:)); caxis(cl);
%         title('midShiftCell omitLaps');
%         
%         subplot(2,2,2);
%         plot(mean(posRatesRef(cueCellInd2,:),1), 'b');
%         hold on;
%         plot(mean(posRatesOmit(cueCellInd2,:),1), 'r');
%         title('avgs');
%         xlabel('pos');
%         ylabel('mean rate (Hz)');
%         legend('cue laps', 'Shift laps');
%         catch
%         end
%     end
