function [sameCellTuningStruc] = sameCellTuning(multSessSegStruc, cell_registered_struct, toPlot);

%% USAGE: [sameCellTuningStruc] = sameCellTuning(multSessTuningStruc, cell_registered_struct);
%
% Questions:
% - are place cells more likely to be present in all sessions?
% - is there a change in average tuning strength over sessions? (on avg,
% for strong place cells?)
% - i.e. how do cells transition between groups?

% ziv array of indices of registered cells in each sessions (with respect to goodSegs,
% usually; thus must reference to goodSeg array to get orig C/A indices)
mapInd = cell_registered_struct.cell_to_index_map;
% NOTE: placeCellStruc is for greatSegs,
% whereas cellReg is for goodSegs

% spatial footprints of all ziv array cells
unitSpatCell = cell_registered_struct.spatial_footprints_corrected;
zivCentroids = cell_registered_struct.registered_cells_centroids;

%%
% for all sessions
for i = 1:size(mapInd,2)
    
    % find index of cellReg cells with respect to original C/A
    goodSeg = multSessSegStruc(i).goodSeg;
    cellRegArr = mapInd(:,i); % mapInd/cellRegArr are in terms of goodSeg indices
    
    % find orig CA index of registered cells 
    for j = 1:length(cellRegArr) % for all cells in cellReg array from this session (0's for cells not in session)
        if cellRegArr(j)~=0  % if this cell is present in this session
            mapInd2(j,i) = goodSeg(cellRegArr(j));
        else
            mapInd2(j,i) = 0; % NOTE: mapInd2 is now ziv cellReg array with orig C/A indices
        end
    end
    
    % find ind of place cells w. re. to orig C/A
%     greatSeg = multSessTuningStruc(i).goodSegPosPkStruc.greatSeg;
%     placeCellOrigInd{i} = greatSeg(multSessTuningStruc(i).placeCellStruc.goodRay);
    
    % find place cell info for cells in this session (from goodSegs)
    disp('Identifying place fields, moving epochs'); tic;
    [out, PCLappedSess] = wrapAndresPlaceFieldsClay(multSessSegStruc(i).pksCell(multSessSegStruc(i).goodSeg), 0, multSessSegStruc(i).treadBehStruc);
    placeCellOrigInd{i} = goodSeg(find(PCLappedSess.Shuff.isPC));
    multSessSegStruc(i).outPC = out;
    multSessSegStruc(i).PCLapSess = PCLappedSess;
    toc;
    
%     disp('And non-movement epochs (and rewCells)'); tic;
%     [outNon, PCLappedSess] = wrapAndresPlaceFieldsClay(multSessSegStruc(i).pksCell(multSessSegStruc(i).goodSeg), 0, multSessSegStruc(i).treadBehStruc, -3);
%     multSessSegStruc(i).outNonPC = outNon;
%     RewLoc =[(round(mean((multSessSegStruc(i).treadBehStruc.rewZoneStartPos)/20)):1:(round(mean((multSessSegStruc(i).treadBehStruc.rewZoneStopPos)/20))))];
%     multSessSegStruc(i).RewLoc = RewLoc;
%     rewCellOrigInd{i} = goodSeg(find(nanmean(outNon.posRates(:,RewLoc),2)>0));
%     toc;
end

%%

disp('Finding cells (and their tuning) in diff sessions'); tic;

% find ziv array cells present in all sessions
cellRegIndInAll = find(min(mapInd2, [], 2)); % [1,1,1,...] in all col (i.e. none have zeros)
cellsInAllOrig = mapInd2(cellRegIndInAll,:); % orig C/A index of all ziv registered cells present in all sessions
cellsInAll = mapInd(cellRegIndInAll,:);

%%
% see if cells present in all sessions are place cells in all
for i = 1:size(cellsInAllOrig,1) % for all cells 
    for j = 1:size(cellsInAllOrig,2)    % for each session from that cell
        if find(placeCellOrigInd{j}== cellsInAllOrig(i,j)) % see if it's a place cell
            sameCellPlaceBool(i,j) = 1;
        else
            sameCellPlaceBool(i,j) = 0;
        end
    end
end

%%
% cellRegInd(cellsInAll) for cells present in all sessions, that are place cells in all
placeCellsInAll = find(min(sameCellPlaceBool, [], 2)); % index in array of only place cells present in all sessions
placeCellsInNone = find(~max(sameCellPlaceBool, [], 2)); % or cells present in all sessions that are place cells in none
placeCellsInAny = find(max(sameCellPlaceBool,[],2));% or cells present in all sessions that are place cells in at least one
%%
% tuning in placeCellsInAll
placeCellAllGoodSegInd = cellsInAll(placeCellsInAll,:);
placeCellAllOrigInd = cellsInAllOrig(placeCellsInAll,:);
placeCellInNoneOrigInd = cellsInAllOrig(placeCellsInNone,:);
placeCellInAnyOrigInd = cellsInAllOrig(placeCellsInAny,:);
%= cellRegIndInAll(placeCellsInAll,:);


%% Save useful vars to output struc

sameCellTuningStruc.multSessSegStruc = multSessSegStruc; % just save orig struc (not too huge)
sameCellTuningStruc.unitSpatCell = unitSpatCell;  % cell array of spatial profiles of ziv cells
sameCellTuningStruc.zivCentroids = zivCentroids;    % centroids of these cells
sameCellTuningStruc.placeCellOrigInd = placeCellOrigInd;  % ind of place cells (goodRay) w. re. to orig C/A
sameCellTuningStruc.cellsInAll = cellsInAll;
sameCellTuningStruc.cellsInAllOrig = cellsInAllOrig; % orig C/A index of all ziv registered cells present in all sessions
sameCellTuningStruc.placeCellAllGoodSegInd = placeCellAllGoodSegInd;
sameCellTuningStruc.placeCellAllOrigInd = placeCellAllOrigInd; % orig C/A index of all cells that are place cells in all sessions
sameCellTuningStruc.placeCellInNoneOrigInd = placeCellInNoneOrigInd; % orig C/A index of all cells that are not place cells in all sessions
sameCellTuningStruc.placeCellInAnyOrigInd = placeCellInAnyOrigInd; % orig C/A index of all cells that are place cells in at least one session
sameCellTuningStruc.regMapOrigInd = mapInd2; % orig C/A ind for all cells in ziv mat
sameCellTuningStruc.regMapGoodSegInd = mapInd;
sameCellTuningStruc.sameCellPlaceBool = sameCellPlaceBool;  % boolean for this mat

%% Plotting
% plot tuning of cells that are place cells (outPC.Shuff.isPC) in all
% sessions

figure; color = {'r' 'g' 'b'};
%title('Place cells in all sessions');
% for all cells
dim = ceil(sqrt(length(placeCellsInAll)));
for i = 1:length(placeCellsInAll)
    subplot(dim,dim,i);
    hold on;
    for j = 1:size(placeCellAllOrigInd,2)
        origInd = placeCellAllOrigInd(i,j); % orig C/A index of this cell
        goodSegInd = find(multSessSegStruc(j).goodSeg == origInd);
        posPks = multSessSegStruc(j).PCLapSess.posRates(goodSegInd,:);
        samePlacePosPks(i,j,:) = posPks;
        plot(posPks, color{j});
    end
end
title('Place cells in all sessions');

%%
% %save spatial profiles of cells from all sessions, aligned with ziv
% for i = 1:length(unitSpatCell)
%     allCell = mean(unitSpatCell{i},1);
%     rgb(:,:,i) = allCell/max(allCell(:));
% end
% 
% %imwrite(rgb, 'cellRegRGB.tif');
% 
% %cellRegIndInAll
% 
% % mark centroids of place cells in all
% try
% rgb2 = insertMarker(10*rgb, zivCentroids(cellRegIndInAll(placeCellsInAll),:));
% %rgb2 = insertMarker(10*rgb, zivCentroids);
% catch
%     rgb2=rgb;
% end
% 
% figure;
% imshow(rgb2);

% also, want to look at whether there are clusters or domains of tuning
% similarity
% 1.) find all place cells, 2.) find tuning center for each place cell, 3.)
% 'color', color, from array of positions, 4.) 

%% for all units in multiple sessions (not just all sessions), find tuning
% similarity; go ahead and get posPks for all ziv cells, to look at tuning of cells
% that appear or disappear in different sessions

% PLOT posRates for all sameCells over sessions
% for all cells
if toPlot == 2
for i = 1:size(mapInd2,1)
    for j = 1:size(mapInd2,2)
        try
        origInd = mapInd2(i,j); % orig C/A index of this cell
        if origInd ~= 0
        goodSegInd = find(multSessSegStruc(j).goodSeg == origInd);
        posPks = multSessSegStruc(j).outPC.posRates(goodSegInd,:);
        else
           posPks = zeros(40,1); 
        end
        allPlacePosPks(i,j,:) = posPks;
        catch
        end
    end
end

% PLOT only place cells in all sesions
% posRates over sessions for placeCellsInAll
color = {'r' 'g' 'b'};
%figure;
for i = 1:size(allPlacePosPks,1)
    %ceil(size(allPlacePosPks,1)/25)
    if mod(i,25)-1==0
        figure;
        n = 1;
    else
        n=n+1;
    end
    subplot(5,5,n);
   for j = 1:size(allPlacePosPks,2)
       hold on;
        plot(squeeze(allPlacePosPks(i,j,:)), color{j});
   end
end
end

%% 
% and from this I see that some cells seem to have consistent apparent place fields
% between sessions, but aren't considered place cells in some, so maybe I should check and see if cells are place cells
% when considering all sessions from cellReg registered cells

%figure; pie([229 35 8 42]);
