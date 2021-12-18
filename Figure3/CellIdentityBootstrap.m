function [BootstrapStruc] = CellIdentityBootstrap (sameCellCueShiftTuningStruc)
% shuffle commented out on 8/19/19 to calculate only overlapping cells
% (denominator)-to correct error 
%%
mapInd = sameCellCueShiftTuningStruc.regMapInd;
placeCellInd = sameCellCueShiftTuningStruc.placeCellInd;
MidCueCellInd = sameCellCueShiftTuningStruc.MidCueCellInd;
EdgeCueCellInd = sameCellCueShiftTuningStruc.EdgeCueCellInd;
nonCueCellInd = sameCellCueShiftTuningStruc.nonCueCellInd;
cellsInAll = sameCellCueShiftTuningStruc.cellsInAll;
TotalCueCellInd = {};
for i = 1:length(MidCueCellInd)
    TotalCueCellInd{i} = unique([MidCueCellInd{i}; EdgeCueCellInd{i}]);
end
%to compare 2-3 only
placeCellIndLastTwo ={}; placeCellIndLastTwo{1,1} = placeCellInd {1,2}; placeCellIndLastTwo{1,2} = placeCellInd {1,3};
MidCueCellIndLastTwo ={}; MidCueCellIndLastTwo{1,1} = MidCueCellInd {1,2}; MidCueCellIndLastTwo{1,2} = MidCueCellInd {1,3};
EdgeCueCellIndLastTwo ={}; EdgeCueCellIndLastTwo{1,1} = EdgeCueCellInd {1,2}; EdgeCueCellIndLastTwo{1,2} = EdgeCueCellInd {1,3};
TotalCueCellIndLastTwo ={}; TotalCueCellIndLastTwo{1,1} =TotalCueCellInd {1,2}; TotalCueCellIndLastTwo{1,2} = TotalCueCellInd {1,3};
nonCueCellIndLastTwo ={}; nonCueCellIndLastTwo{1,1} = nonCueCellInd {1,2}; nonCueCellIndLastTwo{1,2} = nonCueCellInd {1,3};

% find ziv array cells present in sessions ref session is 2
cellRegIndInAll = find(min(mapInd, [], 2)); % [1,1,1,...] in all col (i.e. none have zeros)
cellsInAll = mapInd(cellRegIndInAll,:);
cellRegIndFirstTwo=find(min(mapInd(:,1:2), [], 2));
cellsInFirstTwo = mapInd(cellRegIndFirstTwo,1:2);
cellRegIndLastTwo=find(min(mapInd(:,2:3), [], 2));
cellsInLastTwo = mapInd(cellRegIndLastTwo,2:3);

BootstrapStruc =[];
tic
%% Place Cells
% frac of pcs in sess2 overlapping with sess1 and sess3 pcs
% numPCInFirstTwoAllAnim = [];

sameCell12placeBool=[];
placeCellsFirstTwoInd =[];
PCfrac12 =[];
for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
    for j = 1:2    % for each session from that cell
        if find(placeCellInd{j}== cellsInFirstTwo(i,j)) % see if it's a place cell
            sameCell12placeBool(i,j) = 1;
        else
            sameCell12placeBool(i,j) = 0;
        end
    end
end
% cellRegInd(cellsInAll) for cells present in all sessions, that are place cells in all
placeCellsInFirstTwo= find(min(sameCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
placeCellsFirstTwoInd = cellsInFirstTwo(placeCellsInFirstTwo,:);
PC1only = find((sameCell12placeBool(:, 1) == 1 & sameCell12placeBool(:, 2) == 0));
PC2onlyin1= find((sameCell12placeBool(:, 1) == 0 & sameCell12placeBool(:, 2) == 1));
PC1= find((sameCell12placeBool(:, 1)==1)); 
PC2in1= find((sameCell12placeBool(:, 2) == 1));
PCfrac12 = length(placeCellsFirstTwoInd)/max(length(PC1), length(PC2in1));

rng('shuffle')
RandCell12placeBool =[];
placeCellsInFirstTwoShuff=[];
PCfrac12Shuff = [];
for sh = 1:10000
    sameCell12placeRand = cellsInFirstTwo;
    sameCell12placeRand(:, 1) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 1);
    sameCell12placeRand(:, 2) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 2);
    for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(placeCellInd{j}== sameCell12placeRand(i,j)) % see if it's a place cell
                RandCell12placeBool(i,j) = 1;
            else
                RandCell12placeBool(i,j) = 0;
            end
        end
    end
    placeCellsInFirstTwoShuff = find(min(RandCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(placeCellsInFirstTwoShuff)
        placeCellsInFirstTwoIndShuff = cellsInFirstTwo(placeCellsInFirstTwoShuff,:);
        PCfrac12Shuff(sh) = length(placeCellsInFirstTwoIndShuff)/max(length(PC1), length(PC2in1));
    else
        placeCellsInFirstTwoIndShuff = [];
        PCfrac12Shuff(sh) = 0;
    end
end

sameCell23placeBool = [];
placeCellsLastTwoInd=[];
PCfrac23 =[];
for i = 1:size(cellsInLastTwo,1) % for all cells tracked across last two sessions
    for j = 1:2
        if find(placeCellIndLastTwo{j}== cellsInLastTwo(i,j)) % see if it's a place cell
            sameCell23placeBool(i,j) = 1;
        else
            sameCell23placeBool(i,j) = 0;
        end
    end
end
placeCellsInLastTwo= find(min(sameCell23placeBool, [], 2)); % index in array of only place cells present in all sessions
placeCellsLastTwoInd = cellsInLastTwo(placeCellsInLastTwo,:);
PC2onlyin3 = find((sameCell23placeBool(:, 1) == 1 & sameCell23placeBool(:, 2) == 0));
PC3only= find((sameCell23placeBool(:, 1) == 0 & sameCell23placeBool(:, 2) == 1));
PC2in3= find((sameCell23placeBool(:, 1)==1)); 
PC3= find((sameCell23placeBool(:, 2) == 1));
PCfrac23 = length(placeCellsLastTwoInd)/max(length(PC2in3), length(PC3));

RandCell23placeBool = [];
placeCellsInLastTwoShuff=[];
PCfrac23Shuff = [];
for sh = 1:10000
    sameCell23PlaceRand = cellsInLastTwo;
    sameCell23PlaceRand(:, 1) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 1);
    sameCell23PlaceRand(:, 2) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 2);
    for i = 1:size(cellsInLastTwo,1)
        for j = 1:2
            if find(placeCellIndLastTwo{j}== sameCell23PlaceRand(i,j))
                RandCell23placeBool(i,j) = 1;
            else
                RandCell23placeBool(i,j) = 0;
            end
        end
    end
    placeCellsInLastTwoShuff = find(min(RandCell23placeBool, [], 2));
    if ~isempty(placeCellsInLastTwoShuff)
        placeCellsInLastTwoIndShuff = cellsInLastTwo(placeCellsInLastTwoShuff,:);
        PCfrac23Shuff(sh) = length(placeCellsInLastTwoIndShuff)/max(length(PC2in3), length(PC3));
    else
        placeCellsInLastTwoIndShuff = [];
        PCfrac23Shuff(sh) = 0;
    end
end


figure;
subplot (1,2,1); hist(PCfrac12Shuff, 20);
mean(PCfrac12Shuff < PCfrac12);
hold on; plot([PCfrac12, PCfrac12], [0, 1000], '--r');
bounds = prctile(PCfrac12Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('Place Cell Overlap Fraction in Sess12');
PCfrac12P= 1 - mean(PCfrac12Shuff < PCfrac12);
subplot (1,2,2); hist(PCfrac23Shuff, 20);
mean(PCfrac23Shuff < PCfrac23);
hold on; plot([PCfrac23, PCfrac23], [0, 1000], '--r');
bounds = prctile(PCfrac23Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('Place Cell Overlap Fraction in Sess23');
PCfrac23P= 1 - mean(PCfrac23Shuff < PCfrac23);
%% Mid Cue Cells
% frac of MCs in sess2 overlapping with sess1 and sess3 MCs
sameCell12placeBool = []; sameCell23placeBool = [];
MidCueCellsFirstTwoInd =[]; MidCueCellsLastTwoInd =[];
MCfrac12 =[]; MCfrac23 =[];
for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
    for j = 1:2    % for each session from that cell
        if find(MidCueCellInd{j}== cellsInFirstTwo(i,j))
            sameCell12placeBool(i,j) = 1;
        else
            sameCell12placeBool(i,j) = 0;
        end
    end
end
for i = 1:size(cellsInLastTwo,1) % for all cells
    for j = 1:2
        if find(MidCueCellIndLastTwo{j}== cellsInLastTwo(i,j))
            sameCell23placeBool(i,j) = 1;
        else
            sameCell23placeBool(i,j) = 0;
        end
    end
end
MidCueCellsInFirstTwo= find(min(sameCell12placeBool, [], 2));
MidCueCellsFirstTwoInd = cellsInFirstTwo(MidCueCellsInFirstTwo,:);
MidCueCellsInLastTwo= find(min(sameCell23placeBool, [], 2));
MidCueCellsLastTwoInd = cellsInLastTwo(MidCueCellsInLastTwo,:);
MC1only = find((sameCell12placeBool(:, 1) == 1 & sameCell12placeBool(:, 2) == 0));
MC2onlyin1= find((sameCell12placeBool(:, 1) == 0 & sameCell12placeBool(:, 2) == 1));
MC1= find((sameCell12placeBool(:, 1)==1)); 
MC2in1= find((sameCell12placeBool(:, 2) == 1));
MCfrac12 = length(MidCueCellsFirstTwoInd)/max(length(MC1), length(MC2in1));
MC2onlyin3 = find((sameCell23placeBool(:, 1) == 1 & sameCell23placeBool(:, 2) == 0));
MC3only= find((sameCell23placeBool(:, 1) == 0 & sameCell23placeBool(:, 2) == 1));
MC2in3= find((sameCell23placeBool(:, 1)==1)); 
MC3= find((sameCell23placeBool(:, 2) == 1));
MCfrac23 = length(MidCueCellsLastTwoInd)/max(length(MC2in3), length(MC3));


% now shuffle
RandCell12placeBool =[]; RandCell23placeBool =[];
MidCueCellsInFirstTwoIndShuff=[]; MidCueCellsInLastTwoIndShuff=[];
MCfrac12Shuff = []; MCfrac23Shuff = [];
for sh = 1:10000
    sameCell12MCRand = cellsInFirstTwo;
    sameCell12MCRand(:, 1) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 1);
    sameCell12MCRand(:, 2) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 2);
    for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(MidCueCellInd{j}== sameCell12MCRand(i,j)) % see if it's a place cell
                RandCell12placeBool(i,j) = 1;
            else
                RandCell12placeBool(i,j) = 0;
            end
        end
    end
    MidCueCellsInFirstTwoShuff = find(min(RandCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(MidCueCellsInFirstTwoShuff)
        MidCueCellsInFirstTwoIndShuff = cellsInFirstTwo(MidCueCellsInFirstTwoShuff,:);
        MCfrac12Shuff(sh) = length(MidCueCellsInFirstTwoIndShuff)/max(length(MC1), length(MC2in1));
    else
        MidCueCellsInFirstTwoIndShuff=[];
        MCfrac12Shuff(sh)=0;
    end
end

for sh = 1:10000
    sameCell23MCRand = cellsInLastTwo;
    sameCell23MCRand(:, 1) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 1);
    sameCell23MCRand(:, 2) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 2);
    for i = 1:size(cellsInLastTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(MidCueCellIndLastTwo{j}== sameCell23MCRand(i,j)) % see if it's a place cell
                RandCell23placeBool(i,j) = 1;
            else
                RandCell23placeBool(i,j) = 0;
            end
        end
    end
    MidCueCellsInLastTwoShuff = find(min(RandCell23placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(MidCueCellsInLastTwoShuff)
        MidCueCellsInLastTwoIndShuff = cellsInLastTwo(MidCueCellsInLastTwoShuff,:);
        MCfrac23Shuff(sh) = length(MidCueCellsInLastTwoIndShuff)/max(length(MC2in3), length(MC3));
    else
        MidCueCellsInLastTwoIndShuff=[];
        MCfrac23Shuff(sh) = 0;
    end
end

figure; subplot (1,2,1); hist(MCfrac12Shuff, 20);
mean(MCfrac12Shuff < MCfrac12);
hold on; plot([MCfrac12, MCfrac12], [0, 1000], '--r');
bounds = prctile(MCfrac12Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('MidCue Cell Overlap Fraction in Sess12');
MCfrac12P= 1 - mean(MCfrac12Shuff < MCfrac12);
subplot (1,2,2); hist(MCfrac23Shuff, 20);
mean(MCfrac23Shuff < MCfrac23);
hold on; plot([MCfrac23, MCfrac23], [0, 1000], '--r');
bounds = prctile(MCfrac23Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('MidCue Cell Overlap Fraction in Sess23');
MCfrac23P= 1 - mean(MCfrac23Shuff < MCfrac23);
%% Edge Cue Cells
% frac of ECs in sess2 overlapping with sess1 and sess3 ECs
sameCell12placeBool = []; sameCell23placeBool = [];
EdgeCueCellsFirstTwoInd =[]; EdgeCueCellsLastTwoInd =[];
ECfrac12 =[]; ECfrac23 =[];
for i = 1:size(cellsInFirstTwo,1)
    for j = 1:2
        if find(EdgeCueCellInd{j}== cellsInFirstTwo(i,j)) % see if it's an edge cell
            sameCell12placeBool(i,j) = 1;
        else
            sameCell12placeBool(i,j) = 0;
        end
    end
end
for i = 1:size(cellsInLastTwo,1) % for all cells
    for j = 1:2
        if find(EdgeCueCellIndLastTwo{j}== cellsInLastTwo(i,j))
            sameCell23placeBool(i,j) = 1;
        else
            sameCell23placeBool(i,j) = 0;
        end
    end
end
EdgeCueCellsInFirstTwo= find(min(sameCell12placeBool, [], 2));
EdgeCueCellsFirstTwoInd = cellsInFirstTwo(EdgeCueCellsInFirstTwo,:);
EdgeCueCellsInLastTwo= find(min(sameCell23placeBool, [], 2));
EdgeCueCellsLastTwoInd = cellsInLastTwo(EdgeCueCellsInLastTwo,:);
EC1only = find((sameCell12placeBool(:, 1) == 1 & sameCell12placeBool(:, 2) == 0));
EC2onlyin1= find((sameCell12placeBool(:, 1) == 0 & sameCell12placeBool(:, 2) == 1));
EC1= find((sameCell12placeBool(:, 1)==1)); 
EC2in1= find((sameCell12placeBool(:, 2) == 1));
EC2onlyin3 = find((sameCell23placeBool(:, 1) == 1 & sameCell23placeBool(:, 2) == 0));
EC3only= find((sameCell23placeBool(:, 1) == 0 & sameCell23placeBool(:, 2) == 1));
EC2in3= find((sameCell23placeBool(:, 1)==1)); 
EC3= find((sameCell23placeBool(:, 2) == 1));
ECfrac12 = length(EdgeCueCellsFirstTwoInd)/max(length(EC1), length(EC2in1));
ECfrac23 = length(EdgeCueCellsLastTwoInd)/max(length(EC2in3), length(EC3));

% now shuffle
RandCell12placeBool =[]; RandCell23placeBool =[];
EdgeCueCellsInFirstTwoIndShuff=[]; EdgeCueCellsInLastTwoIndShuff=[];
ECfrac12Shuff = []; ECfrac23Shuff = [];
for sh = 1:10000
    sameCell12ECRand = cellsInFirstTwo;
    sameCell12ECRand(:, 1) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 1);
    sameCell12ECRand(:, 2) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 2);
    for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(EdgeCueCellInd{j}== sameCell12ECRand(i,j)) % see if it's a place cell
                RandCell12placeBool(i,j) = 1;
            else
                RandCell12placeBool(i,j) = 0;
            end
        end
    end
    EdgeCueCellsInFirstTwoShuff = find(min(RandCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(EdgeCueCellsInFirstTwoShuff)
        EdgeCueCellsInFirstTwoIndShuff = cellsInFirstTwo(EdgeCueCellsInFirstTwoShuff,:);
        ECfrac12Shuff(sh) = length(EdgeCueCellsInFirstTwoIndShuff)/max(length(EC1), length(EC2in1));
    else
        EdgeCueCellsInFirstTwoIndShuff = [];
        ECfrac12Shuff(sh) = 0;
    end
end
for sh = 1:10000
    sameCell23ECRand = cellsInLastTwo;
    sameCell23ECRand(:, 1) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 1);
    sameCell23ECRand(:, 2) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 2);
    for i = 1:size(cellsInLastTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(EdgeCueCellIndLastTwo{j}== sameCell23ECRand(i,j)) % see if it's a place cell
                RandCell23placeBool(i,j) = 1;
            else
                RandCell23placeBool(i,j) = 0;
            end
        end
    end
    EdgeCueCellsInLastTwoShuff = find(min(RandCell23placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(EdgeCueCellsInLastTwoShuff)
        EdgeCueCellsInLastTwoIndShuff = cellsInLastTwo(EdgeCueCellsInLastTwoShuff,:);
        ECfrac23Shuff(sh) = length(EdgeCueCellsInLastTwoIndShuff)/max(length(EC2in3), length(EC3));
    else
        EdgeCueCellsInLastTwoIndShuff=[];
        ECfrac23Shuff(sh)=0;
    end
end
figure; subplot (1,2,1); hist(ECfrac12Shuff, 20);
mean(ECfrac12Shuff < ECfrac12);
hold on; plot([ECfrac12, ECfrac12], [0, 1000], '--r');
bounds = prctile(ECfrac12Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('EdgeCue Cell Overlap Fraction in Sess12');
ECfrac12P= 1 - mean(ECfrac12Shuff < ECfrac12);
subplot (1,2,2); hist(ECfrac23Shuff, 20);
mean(ECfrac23Shuff < ECfrac23);
hold on; plot([ECfrac23, ECfrac23], [0, 1000], '--r');
bounds = prctile(ECfrac23Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('EdgeCue Cell Overlap Fraction in Sess23');
ECfrac23P= 1 - mean(ECfrac23Shuff < ECfrac23);
%% non Cue Cells
% frac of nonCues in sess2 overlapping with sess1 and sess3 NCs
sameCell12placeBool = []; sameCell23placeBool = [];
nonCueCellsFirstTwoInd =[]; nonCueCellsLastTwoInd =[];
NCfrac12 =[]; NCfrac23 =[];
for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
    for j = 1:2    % for each session from that cell
        if find(nonCueCellInd{j}== cellsInFirstTwo(i,j)) % see if it's a place cell
            sameCell12placeBool(i,j) = 1;
        else
            sameCell12placeBool(i,j) = 0;
        end
    end
end
for i = 1:size(cellsInLastTwo,1) % for all cells
    for j = 1:2
        if find(nonCueCellIndLastTwo{j}== cellsInLastTwo(i,j)) % see if it's a place cell
            sameCell23placeBool(i,j) = 1;
        else
            sameCell23placeBool(i,j) = 0;
        end
    end
end
% cellRegInd(cellsInAll) for cells present in all sessions, that are place cells in all
nonCueCellsInFirstTwo= find(min(sameCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
nonCueCellsFirstTwoInd = cellsInFirstTwo(nonCueCellsInFirstTwo,:);
nonCueCellsInLastTwo= find(min(sameCell23placeBool, [], 2)); % index in array of only place cells present in all sessions
nonCueCellsLastTwoInd = cellsInLastTwo(nonCueCellsInLastTwo,:);
NC1only = find((sameCell12placeBool(:, 1) == 1 & sameCell12placeBool(:, 2) == 0));
NC2onlyin1= find((sameCell12placeBool(:, 1) == 0 & sameCell12placeBool(:, 2) == 1));
NC1= find((sameCell12placeBool(:, 1)==1)); 
NC2in1= find((sameCell12placeBool(:, 2) == 1));
NC2onlyin3 = find((sameCell23placeBool(:, 1) == 1 & sameCell23placeBool(:, 2) == 0));
NC3only= find((sameCell23placeBool(:, 1) == 0 & sameCell23placeBool(:, 2) == 1));
NC2in3= find((sameCell23placeBool(:, 1)==1)); 
NC3= find((sameCell23placeBool(:, 2) == 1));
NCfrac12 = length(nonCueCellsFirstTwoInd)/max(length(NC1), length(NC2in1));
NCfrac23 = length(nonCueCellsLastTwoInd)/max(length(NC2in3), length(NC3));

% now shuffle
RandCell12placeBool =[]; RandCell23placeBool =[];
nonCueCellsInFirstTwoIndShuff=[]; nonCueCellsInLastTwoIndShuff=[];
NCfrac12Shuff = []; NCfrac23Shuff = [];
for sh = 1:10000
    sameCell12NCRand = cellsInFirstTwo;
    sameCell12NCRand(:, 1) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 1);
    sameCell12NCRand(:, 2) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 2);
    for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(nonCueCellInd{j}== sameCell12NCRand(i,j)) % see if it's a place cell
                RandCell12placeBool(i,j) = 1;
            else
                RandCell12placeBool(i,j) = 0;
            end
        end
    end
    nonCueCellsInFirstTwoShuff = find(min(RandCell12placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(nonCueCellsInFirstTwoShuff)
        nonCueCellsInFirstTwoIndShuff = cellsInFirstTwo(nonCueCellsInFirstTwoShuff,:);
        NCfrac12Shuff(sh) = length(nonCueCellsInFirstTwoIndShuff)/max(length(NC1), length(NC2in1));
    else
        nonCueCellsInFirstTwoIndShuff=[];
        NCfrac12Shuff(sh)=0;
    end
end
for sh = 1:10000
    sameCell23NCRand = cellsInLastTwo;
    sameCell23NCRand(:, 1) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 1);
    sameCell23NCRand(:, 2) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 2);
    for i = 1:size(cellsInLastTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(nonCueCellIndLastTwo{j}== sameCell23NCRand(i,j)) % see if it's a place cell
                RandCell23placeBool(i,j) = 1;
            else
                RandCell23placeBool(i,j) = 0;
            end
        end
    end
    nonCueCellsInLastTwoShuff = find(min(RandCell23placeBool, [], 2)); % index in array of only place cells present in all sessions
    if ~isempty(nonCueCellsInLastTwoShuff)
        nonCueCellsInLastTwoIndShuff = cellsInLastTwo(nonCueCellsInLastTwoShuff,:);
        NCfrac23Shuff(sh) = length(nonCueCellsInLastTwoIndShuff)/max(length(NC2in3), length(NC3));

    else
        nonCueCellsInLastTwoIndShuff=[];
        NCfrac23Shuff(sh)=0;
    end
end
figure; subplot (1,2,1); hist(NCfrac12Shuff, 20);
mean(NCfrac12Shuff < NCfrac12);
hold on; plot([NCfrac12, NCfrac12], [0, 1000], '--r');
bounds = prctile(NCfrac12Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('nonCue Cell Overlap Fraction in Sess12');
NCfrac12P= 1 - mean(NCfrac12Shuff < NCfrac12);
subplot (1,2,2); hist(NCfrac23Shuff, 20);
mean(NCfrac23Shuff < NCfrac23);
hold on; plot([NCfrac23, NCfrac23], [0, 1000], '--r');
bounds = prctile(NCfrac23Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('nonCue Cell Overlap Fraction in Sess23');
NCfrac23P= 1 - mean(NCfrac23Shuff < NCfrac23);
%% Overlap of TotalCue Cells with NonCueCells
sameCell12CueBool = []; sameCell23CueBool = [];
sameCell12NonCueBool = []; sameCell23NonCueBool = [];
diffTypeCellsInFirstTwoInd =[]; diffTypeCellsInLastTwoInd =[];
NonandCuefrac12 =[]; NonandCuefrac23 =[];
for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
    for j = 1:2    % for each session from that cell
        if find(TotalCueCellInd{j}== cellsInFirstTwo(i,j)) % see if it's a cue cell
            sameCell12CueBool(i,j) = 1;
        else
            sameCell12CueBool(i,j) = 0;
        end
        if find(nonCueCellInd{j}== cellsInFirstTwo(i,j)) % see if it's a non cue cell
            sameCell12NonCueBool(i,j) = 1;
        else
            sameCell12NonCueBool(i,j) = 0;
        end
    end
end
for i = 1:size(cellsInLastTwo,1) %each cell in 2-3
    for j = 1:2    % for each session from that cell
        if find(TotalCueCellIndLastTwo{j}== cellsInLastTwo(i,j)) % cue cell
            sameCell23CueBool(i,j) = 1;
        else
            sameCell23CueBool(i,j) = 0;
        end
        if find(nonCueCellIndLastTwo{j}== cellsInLastTwo(i,j)) % non cue
            sameCell23NonCueBool(i,j) = 1;
        else
            sameCell23NonCueBool(i,j) = 0;
        end
    end
end

diffTypeCellsInFirstTwo= find((sameCell12CueBool(:, 1) == 1 & sameCell12NonCueBool(:, 2) == 1) | ...
    ((sameCell12CueBool(:, 2) == 1 & sameCell12NonCueBool(:, 1) == 1))); %indices that are cue cells in sess1 and non cue in 2 or vice versa
diffTypeCellsInFirstTwoInd = cellsInFirstTwo(diffTypeCellsInFirstTwo,:);
diffTypeCellsLastTwo= find((sameCell23CueBool(:, 1) == 1 & sameCell23NonCueBool(:, 2) == 1) | ...
    ((sameCell23CueBool(:, 2) == 1 & sameCell23NonCueBool(:, 1) == 1)));
diffTypeCellsInLastTwoInd = cellsInLastTwo(diffTypeCellsLastTwo,:);
NonandCuefrac12 = length(diffTypeCellsInFirstTwoInd)/max(length(PC1), length(PC2in1));
NonandCuefrac23 = length(diffTypeCellsInLastTwoInd)/max(length(PC2in3), length(PC3));
Cueonly12 = find((sameCell12CueBool(:, 1) == 1 & sameCell12NonCueBool(:, 2) == 0) | ...
    ((sameCell12CueBool(:, 2) == 1 & sameCell12NonCueBool(:, 1) == 0))); %indices that are cue cells in sess1 and in 2 
nonCueonly12 = find((sameCell12CueBool(:, 1) == 0 & sameCell12NonCueBool(:, 2) == 1) | ...
    ((sameCell12CueBool(:, 2) == 0 & sameCell12NonCueBool(:, 1) == 1))); %indices that are cue cells in sess1 and in 2 
Cueonly23 = find((sameCell23CueBool(:, 1) == 1 & sameCell23NonCueBool(:, 2) == 0) | ...
    ((sameCell23CueBool(:, 2) == 1 & sameCell23NonCueBool(:, 1) == 0))); %indices that are cue cells in sess1 and in 2 
nonCueonly23 = find((sameCell23CueBool(:, 1) == 0 & sameCell23NonCueBool(:, 2) == 1) | ...
    ((sameCell23CueBool(:, 2) == 0 & sameCell23NonCueBool(:, 1) == 1))); %indices that are cue cells in sess1 and in 2 

%shuffle
RandsameCell12CueBool = []; RandsameCell12NonCueBool=[]; RandsameCell23CueBool = []; RandsameCell23NonCueBool=[];
diffTypeCellsInFirstTwoIndShuff=[]; diffTypeCellsInLastTwoIndShuff=[];
NonandCuefrac12Shuff = []; NonandCuefrac23Shuff = [];

for sh = 1:10000
    diffCell12placeRand = cellsInFirstTwo;
    diffCell12placeRand(:, 1) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 1);
    diffCell12placeRand(:, 2) = cellsInFirstTwo(randperm(size(cellsInFirstTwo, 1)), 2);
    for i = 1:size(cellsInFirstTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(TotalCueCellInd{j}== diffCell12placeRand(i,j))
                RandsameCell12CueBool(i,j) = 1;
            else
                RandsameCell12CueBool(i,j) = 0;
            end
            if find(nonCueCellInd{j}== cellsInFirstTwo(i,j))
                RandsameCell12NonCueBool(i,j) = 1;
            else
                RandsameCell12NonCueBool(i,j) = 0;
            end
        end
        
    end
    diffTypeCellsInFirstTwoShuff= find((RandsameCell12CueBool(:, 1) == 1 & RandsameCell12NonCueBool(:, 2) == 1) | ...
        ((RandsameCell12CueBool(:, 2) == 1 & RandsameCell12NonCueBool(:, 1) == 1)));
    if ~isempty(diffTypeCellsInFirstTwoShuff)
        diffTypeCellsInFirstTwoIndShuff = cellsInFirstTwo(diffTypeCellsInFirstTwoShuff,:);
        NonandCuefrac12Shuff(sh) = length(diffTypeCellsInFirstTwoIndShuff)/max(length(PC1), length(PC2in1));
    else
        diffTypeCellsInFirstTwoIndShuff=[];
        NonandCuefrac12Shuff(sh)=0;
    end
end
for sh = 1:10000
    diffCell23placeRand = cellsInLastTwo;
    diffCell23placeRand(:, 1) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 1);
    diffCell23placeRand(:, 2) = cellsInLastTwo(randperm(size(cellsInLastTwo, 1)), 2);
    for i = 1:size(cellsInLastTwo,1) %each cell in 1-2
        for j = 1:2    % for each session from that cell
            if find(TotalCueCellIndLastTwo{j}== diffCell23placeRand(i,j))
                RandsameCell23CueBool(i,j) = 1;
            else
                RandsameCell23CueBool(i,j) = 0;
            end
            if find(nonCueCellIndLastTwo{j}== cellsInLastTwo(i,j))
                RandsameCell23NonCueBool(i,j) = 1;
            else
                RandsameCell23NonCueBool(i,j) = 0;
            end
        end
        
    end
    diffTypeCellsInLastTwoShuff= find((RandsameCell23CueBool(:, 1) == 1 & RandsameCell23NonCueBool(:, 2) == 1) | ...
        ((RandsameCell23CueBool(:, 2) == 1 & RandsameCell23NonCueBool(:, 1) == 1)));
    if ~isempty(diffTypeCellsInLastTwoShuff)
        diffTypeCellsInLastTwoIndShuff = cellsInLastTwo(diffTypeCellsInLastTwoShuff,:);
        NonandCuefrac23Shuff(sh) = length(diffTypeCellsInLastTwoShuff)/max(length(PC2in3), length(PC3));
    else
        diffTypeCellsInLastTwoIndShuff=[];
        NonandCuefrac23Shuff(sh)=0;
    end
end
figure; subplot (1,2,1); hist(NonandCuefrac12Shuff, 20);
mean(NonandCuefrac12Shuff < NonandCuefrac12);
hold on; plot([NonandCuefrac12, NonandCuefrac12], [0, 1000], '--r');
bounds = prctile(NonandCuefrac12Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('Non and Cue Overlap Fraction in Sess12');
NonandCuefrac12P= 1 - mean(NonandCuefrac12Shuff < NonandCuefrac12);


subplot (1,2,2); hist(NonandCuefrac23Shuff, 20);
mean(NonandCuefrac23Shuff < NonandCuefrac23);
hold on; plot([NonandCuefrac23, NonandCuefrac23], [0, 1000], '--r');
bounds = prctile(NonandCuefrac23Shuff, [2.5, 97.5]);
hold on; plot([bounds(1), bounds(1)], [0, 1000], '--k');
hold on; plot([bounds(2), bounds(2)], [0, 1000], '--k'); title('Non and Cue Overlap Fraction in Sess23');
NonandCuefrac23P= 1 - mean(NonandCuefrac23Shuff < NonandCuefrac23);
toc
%%
OverlapFrac = [PCfrac12, PCfrac23, MCfrac12, MCfrac23, ECfrac12, ECfrac23, NCfrac12,  NCfrac23, NonandCuefrac12, NonandCuefrac23];
Pvalues = [PCfrac12P, PCfrac23P, MCfrac12P, MCfrac23P, ECfrac12P, ECfrac23P, NCfrac12P,  NCfrac23P, NonandCuefrac12P, NonandCuefrac23P];
PC12= [max(mapInd(:,1)), max(mapInd(:,2)), length(placeCellInd{1,1}), length(placeCellInd{1,2}), length(cellsInFirstTwo), length(placeCellsFirstTwoInd), length(PC1only), length(PC2onlyin1)];
ActNumPC12 = [length(placeCellsFirstTwoInd), max(length(PC1), length(PC2in1))];
ShuffNumPC12 = PCfrac12Shuff * max(length(PC1), length(PC2in1));
PC23= [max(mapInd(:,2)), max(mapInd(:,3)), length(placeCellInd{1,2}), length(placeCellInd{1,3}), length(cellsInLastTwo), length(placeCellsLastTwoInd), length(PC2onlyin3), length(PC3only)];
ActNumPC23 = [length(placeCellsLastTwoInd), max(length(PC2in3), length(PC3))];
ShuffNumPC23 = PCfrac23Shuff * max(length(PC2in3), length(PC3));

MC12= [length(MidCueCellInd{1,1}), length(MidCueCellInd{1,2}), length(MidCueCellsFirstTwoInd), length(MC1only), length(MC2onlyin1)];
ActNumMC12 = [length(MidCueCellsFirstTwoInd), max(length(MC1), length(MC2in1))];
ShuffNumMC12 = MCfrac12Shuff * max(length(MC1), length(MC2in1));
MC23= [length(MidCueCellInd{1,2}), length(MidCueCellInd{1,3}), length(MidCueCellsLastTwoInd), length(MC2onlyin3), length(MC3only)];
ActNumMC23 = [length(MidCueCellsLastTwoInd), max(length(MC2in3), length(MC3))];
ShuffNumMC23 = MCfrac23Shuff * max(length(MC2in3), length(MC3));

EC12= [length(EdgeCueCellInd{1,1}), length(EdgeCueCellInd{1,2}), length(EdgeCueCellsFirstTwoInd), length(EC1only), length(EC2onlyin1)];
EC23= [length(EdgeCueCellInd{1,2}), length(EdgeCueCellInd{1,3}), length(EdgeCueCellsLastTwoInd), length(EC2onlyin3), length(EC3only)];
ActNumEC12 = [length(EdgeCueCellsFirstTwoInd), max(length(EC1), length(EC2in1))];
ShuffNumEC12 = ECfrac12Shuff * max(length(EC1), length(EC2in1));
ActNumEC23 = [length(EdgeCueCellsLastTwoInd), max(length(EC2in3), length(EC3))];
ShuffNumEC23 = ECfrac23Shuff * max(length(EC2in3), length(EC3));

NC12= [length(nonCueCellInd{1,1}), length(nonCueCellInd{1,2}), length(nonCueCellsFirstTwoInd), length(NC1only), length(NC2onlyin1)];
NC23= [length(nonCueCellInd{1,2}), length(nonCueCellInd{1,3}), length(nonCueCellsLastTwoInd), length(NC2onlyin3), length(NC3only)];
ActNumNC12 = [length(nonCueCellsFirstTwoInd), max(length(NC1), length(NC2in1))];
ShuffNumNC12 = NCfrac12Shuff * max(length(NC1), length(NC2in1));
ActNumNC23 = [length(nonCueCellsLastTwoInd), max(length(PC2in3), length(PC3))];
ShuffNumNC23 = NCfrac23Shuff * max(length(PC2in3), length(PC3));

ActNumNonandCue12 = [length(diffTypeCellsInFirstTwoInd), max(length(PC1), length(PC2in1))];
ShuffNumNonandCue12 = NonandCuefrac12Shuff * max(length(PC1), length(PC2in1));
ActNumNonandCue23 = [length(diffTypeCellsInLastTwoInd), max(length(PC2in3), length(PC3))];
ShuffNumNonandCue23 = NonandCuefrac23Shuff * max(length(PC2in3), length(PC3));
NumNonandCue12  =  [length(diffTypeCellsInFirstTwo), length(Cueonly12), length(nonCueonly12)];
NumNonandCue23 =  [ length(diffTypeCellsInLastTwoInd), length(Cueonly23), length(nonCueonly23)];
dirFull = pwd;
spacer = [find(dirFull == '/'), find(dirFull == '\')];
dirName = dirFull((max(spacer) + 1):end);

BootstrapStruc.MultSessID = dirName;
BootstrapStruc.OverlapFrac = OverlapFrac;
BootstrapStruc.Pvalues = Pvalues;
BootstrapStruc.ActNumPC12 = ActNumPC12;
BootstrapStruc.ShuffNumPC12 = ShuffNumPC12;
BootstrapStruc.PC12 = PC12;
BootstrapStruc.PC23 = PC23;
BootstrapStruc.ActNumPC23 = ActNumPC23;
BootstrapStruc.ShuffNumPC23 = ShuffNumPC23;
BootstrapStruc.MC12 = MC12;
BootstrapStruc.MC23 = MC23;
BootstrapStruc.ActNumMC12 = ActNumMC12;
BootstrapStruc.ShuffNumMC12 = ShuffNumMC12;
BootstrapStruc.ActNumMC23 = ActNumMC23;
BootstrapStruc.ShuffNumMC23 = ShuffNumMC23;
BootstrapStruc.EC12 = EC12;
BootstrapStruc.EC23 = EC23;
BootstrapStruc.ActNumEC12 = ActNumEC12;
BootstrapStruc.ShuffNumEC12 = ShuffNumEC12;
BootstrapStruc.ActNumEC23 = ActNumEC23;
BootstrapStruc.ShuffNumEC23 = ShuffNumEC23;
BootstrapStruc.ActNumNC12 = ActNumNC12;
BootstrapStruc.ShuffNumNC12 = ShuffNumNC12;
BootstrapStruc.ActNumNC23 = ActNumNC23;
BootstrapStruc.ShuffNumNC23 = ShuffNumNC23;
BootstrapStruc.NC12 = NC12;
BootstrapStruc.NC23 = NC23;
BootstrapStruc.ActNumNonandCue12 = ActNumNonandCue12;
BootstrapStruc.ShuffNumNonandCue12 = ShuffNumNonandCue12;
BootstrapStruc.ActNumNonandCue23 = ActNumNonandCue23;
BootstrapStruc.ShuffNumNonandCue23 = ShuffNumNonandCue23;
BootstrapStruc.NumNonandCue12 = NumNonandCue12;
BootstrapStruc.NumNonandCue23 = NumNonandCue23;

