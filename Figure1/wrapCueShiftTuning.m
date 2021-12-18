function [cueShiftStruc] = wrapCueShiftTuning(varargin) %pksCell, goodSeg, treadBehStruc,lapTypeInfo, )

% Examples
% [cueShiftStruc] = wrapCueShiftTuning(); % processes current folder, all
% segs in folder segDict (but only plots isPC?)
% [cueShiftStruc] = wrapCueShiftTuning(pksCell); % input pksCell (e.g. from
% quickTuning)
% [cueShiftStruc] = wrapCueShiftTuning(rewOmit); % for cue sessions with
% rewOmit trials
% [cueShiftStruc] = wrapCueShiftTuning(lapTypeInfo); % to input lap
% sequence or lapTypeArr
% [cueShiftStruc] = wrapCueShiftTuning(pksCell, rewOmit, lapTypeInfo); %
% (or any combinations of these?) 

% load goodSeg if present
try
filename = findLatestFilename('_goodSeg_');
load(filename); 
catch
    disp('Cant find goodSeg');
end
disp(['Calculating cue shift tuning for ' filename]);

% load treadBehStruc or create if necessary (but watch out if you have run
% guiSelectDGCs because it might create non-cue treadBehStruc, which will
% be used automatically I think
try
    load(findLatestFilename('treadBehStruc'));
catch
    disp('Couldnt find previous treadBehStruc so processing');
    [treadBehStruc] = procHen2pBehav('auto', 'cue');
end

% parse arguments (and fill in rest with defaults)
if nargin==1
    if iscell(varargin{1})  % if cell arr, then its pksCell (for quickTuning)
        pksCell = varargin{1};
        goodSeg = 1:length(pksCell);
        rewOmit = 0;
    elseif length(varargin{1})==1
        rewOmit = varargin{1};
    else
        lapTypeInfo = varargin{1};
        rewOmit = 0;
    end
elseif length(varargin)==2
    if iscell(varargin{1})  % if cell arr, then its pksCell (for quickTuning)
        pksCell = varargin{1};
        goodSeg = 1:length(pksCell);
        if length(varargin{2})==1
            rewOmit = varargin{2};
        else
            lapTypeInfo = varargin{2};
            rewOmit = 0;
        end
    else
        if length(varargin{1})==1
            rewOmit = varargin{1};
            lapTypeInfo = varargin{2};
        else
            lapTypeInfo = varargin{1};
            rewOmit = varargin{2};
        end
    end
elseif length(varargin)==3
    pksCell = varargin{1};
    goodSeg = 1:length(pksCell);
    if length(varargin{2})==1
        rewOmit = varargin{2};
        lapTypeInfo = varargin{3};
    else
        lapTypeInfo = varargin{2};
        rewOmit = varargin{3};
    end
    
else
    rewOmit=0;
    
end


%% format other stuff

% auto detect temporal downsample
T = treadBehStruc.adjFrTimes;

maxFr = [];
for numSeg = 1:length(pksCell)
    maxFr = max([max(pksCell{numSeg}) maxFr]);
end

% NOTE this is hacky and may not work sometimes
if length(T)>maxFr*1.8
    T = T(1:2:end);
end

% find lap type/epochs
if exist('lapTypeInfo')==0
[pksCellCell, posLapCell, lapCueStruc] = sepCueShiftLapSpkTimes(pksCell, goodSeg, treadBehStruc, rewOmit); %, lapTypeInfo);
else
    [pksCellCell, posLapCell, lapCueStruc] = sepCueShiftLapSpkTimes(pksCell, goodSeg, treadBehStruc, rewOmit, lapTypeInfo);
end

% posLap1 = posLapCell{1}; posLap2 = posLapCell{2};
% pksCell1 = pksCellCell{1}; pksCell2 = pksCellCell{2};
% posLap1 = posLap1/max(posLap1); posLap2 = posLap2/max(posLap2);

% build spike arrays from spk times
for typeNum = 1:length(pksCellCell)
    spikeCell{typeNum} = zeros(length(pksCellCell{typeNum}),length(posLapCell{typeNum}));
    for i = 1:length(pksCellCell{typeNum})
        spikeCell{typeNum}(i,pksCellCell{typeNum}{i})=1;
    end
end

% Calculate tuning for each lap type (and concat struc in cell array)
shuffN = 1000;
for typeNum = 1:length(pksCellCell)
    try
    spikes = spikeCell{typeNum};
    treadPos = posLapCell{typeNum}; treadPos = treadPos/max(treadPos);
    disp(['Calculating tuning for lapType ' num2str(typeNum)]); tic;
    PCLappedSessCell{typeNum} = computePlaceCellsLappedWithEdges3(spikes, treadPos, T(1:length(posLapCell{typeNum})), shuffN);
    toc;
    catch
        PCLappedSessCell{typeNum} = [];
        disp(['Prob with lapType ' num2str(typeNum)]);
    end
end

% disp('Calc lapType2 tuning'); tic;
% spikes = C2; treadPos = posLap2;
% PCLappedSess2 = computePlaceCellsLappedWithEdges3a(spikes, treadPos, T(1:length(posLap2)), shuffN);
% toc;
% 
% cueShiftStruc.PCLappedSess1 = PCLappedSess1;
% cueShiftStruc.PCLappedSess2 = PCLappedSess2;
% cueShiftStruc.pksCell1=pksCell1; cueShiftStruc.posLap1=posLap1; 
% cueShiftStruc.pksCell2=pksCell2; cueShiftStruc.posLap2=posLap2; cueShiftStruc.lapFrInds=lapFrInds;

cueShiftStruc.filename = findLatestFilename('.h5');
cueShiftStruc.pksCellCell = pksCellCell;
cueShiftStruc.posLapCell = posLapCell;
cueShiftStruc.lapCueStruc = lapCueStruc;
cueShiftStruc.PCLappedSessCell = PCLappedSessCell;

% refLapType = 2;
% plotCueShiftStruc(cueShiftStruc, refLapType);

% 
% pc = find(PCLappedSess1.Shuff.isPC==1);
% posRates1 = PCLappedSess1.posRates(pc,:);
% [maxVal, maxInd] = max(posRates1');
% [newInd, oldInd] = sort(maxInd);
% sortInd = oldInd;
% posRates1 = posRates1(sortInd,:);
% 
% posRates2 = PCLappedSess2.posRates(pc,:);
% posRates2 = posRates2(sortInd,:);
% 
% figure('Position', [0 0 800 800]);
% subplot(2,2,1);
% colormap(jet);
% imagesc(posRates1);
% xlabel('position');
% title('LapType1');
% 
% % tuning of lapType2 PCs
% subplot(2,2,2);
% colormap(jet);
% imagesc(posRates2);
% xlabel('position');
% title('LapType2');
% 
% % and mean of each
% subplot(2,2,3);
% plot(mean(posRates1,1));
% hold on;
% plot(mean(posRates2,1),'g');
% title('posRates1=b, posRates2=g');
% 
% 
