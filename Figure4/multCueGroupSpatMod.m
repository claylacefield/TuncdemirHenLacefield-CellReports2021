function [multCueSpatModStruc] = multCueGroupSpatMod(groupCueStrucStruc, toPlot)

% Clay 2020
% Take array of groupCueStrucs from multCue and calc spatMod by day
%
% Use on data Backup20TB/clay/DGdata/190517-22
%
%


% 1. go through each day's groupCueStruc
% 2. run findCueCellsMulti on each cueShiftStruc
% 3. take cue cells (from both first and second cue) based upon omit for
% that cue (see if there are more from first or second?)
% and calc spatMod (rate pref/rate other)
% - also compare other cue response for non-cue cells and cue cells (2x
% rate or shuff?)
% 4. avg spat mod for each day



for d = 1:length(groupCueStrucStruc) % for each day
    
    groupCueStruc = groupCueStrucStruc(d).groupCueStruc;
    
    for m = 1:length(groupCueStruc) % for each mouse in that day
        
        cueShiftStruc = groupCueStruc(m);
        multCueSpatModStruc(d).mouse(m).segDictName = cueShiftStruc.filename;
        multCueSpatModStruc(d).mouse(m).path = cueShiftStruc.path;
        
        % cue cells for first cue
        [cueCellStruc] = findCueCellsMulti(cueShiftStruc, [20,40],2,0);
        cueCellInd1a = cueCellStruc.midCueCellInd; % 2x omit rate cells
        cueCellInd1b = cueCellStruc.midCueCellInd3; % cells signif by shuffle
        nonCueCellInd1 = cueCellStruc.nonCueCellInd; % <2x omit
        
        % cue cells for second cue
        [cueCellStruc] = findCueCellsMulti(cueShiftStruc, [70,90],3,0);
        cueCellInd2a = cueCellStruc.midCueCellInd;
        cueCellInd2b = cueCellStruc.midCueCellInd3;
        nonCueCellInd2 = cueCellStruc.nonCueCellInd;
        
        posRatesAll = cueShiftStruc.PCLappedSessCell{1}.posRates;
        
        fieldNamesCell = {'cueCellInd1a', 'cueCellInd1b', 'nonCueCellInd1', 'cueCellInd2a', 'cueCellInd2b', 'nonCueCellInd2'};
        
        % for each type of cells (cue1/2, noncue)
        if toPlot
            figure; p=0;
        end
        
        for f=1:length(fieldNamesCell)
%             try
                spatMod = [];
                field = fieldNamesCell{f};
                eval(['cellInds = ' field ';']);
                for c = 1:length(cellInds)
                    posRatesN = posRatesAll(cellInds(c),:);
                    cue1rate = max(posRatesN(20:40));
                    cue2rate = max(posRatesN(70:90));
                    if strfind(field,'1')
                        spatMod(c) = cue2rate/cue1rate; % non-preferred/preferred cue loc max rate
                    else
                        spatMod(c) = cue1rate/cue2rate;
                    end
                end
                posRatesCell = posRatesAll(cellInds,:);
                posRateLapCell = squeeze(nanmean(cueShiftStruc.PCLappedSessCell{1}.ByLap.posRateByLap(cellInds,:,:),1));
                posRateLapCell2 = cueShiftStruc.PCLappedSessCell{1}.ByLap.posRateByLap(cellInds,:,:);
                
                % save to output
                multCueSpatModStruc(d).mouse(m).cellType(f).name = fieldNamesCell{f};
                multCueSpatModStruc(d).mouse(m).cellType(f).inds = cellInds;
                multCueSpatModStruc(d).mouse(m).cellType(f).spatMod = spatMod;
                multCueSpatModStruc(d).mouse(m).cellType(f).posRates = posRatesCell;
                multCueSpatModStruc(d).mouse(m).cellType(f).posRatesLap = posRateLapCell;
                multCueSpatModStruc(d).mouse(m).cellType(f).posRatesLapCell = posRateLapCell2;
                
                
                if toPlot
                    p=p+1;
                    subplot(2,6,p); imagesc(squeeze(posRateLapCell)');title(field);
                    p=p+1;
                    if strfind(field,'1') % change depending upon preferred cue
                        r =max(posRateLapCell(70:90,:))./max(posRateLapCell(20:40,:)); subplot(2,6,p); plot(r); title(field);
                    else
                        r =max(posRateLapCell(20:40,:))./max(posRateLapCell(70:90,:)); subplot(2,6,p); plot(r); title(field);
                    end
                end
                
%             catch
%                 disp(['Problem with ' field]);
%             end
%             
        end % end FOR f=each cell type
        
    end  % end FOR m=each mouse
    
end % end FOR d=each day