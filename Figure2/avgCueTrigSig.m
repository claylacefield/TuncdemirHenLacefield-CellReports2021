function [cueTrigSigStruc] = avgCueTrigSig(segNum, eventName, toPlot, varargin)

%% USAGE: [cueTrigSigStruc] = avgCueTrigSig(segNum, eventName, toPlot, varargin)

% for varargin, want to be able to input C, treadBehStruc
% eventName = 'led' to make 'ledTime', etc.

if length(varargin)>0 % like 'auto'
    if ~isempty(strfind(varargin{1},'seg'))
        segDictName = varargin{1};
    elseif varargin{1}==0
        segDictName = findLatestFilename('_segDict_', 'goodSeg');
    elseif varargin{1}==1
        toZ=1;
    end
    
    try
    load(segDictName);
    if ~isempty(strfind(segDictName, 'seg2P'))
        C = seg2P.C2p;
    end
    catch
    end
    
    if length(varargin)==2
        toZ = varargin{2};
    end
    
    if length(varargin)==3
        toZ = varargin{1};
        C = varargin{2};
        treadBehStruc = varargin{3};
    end
else
    segDictName = uigetfile('*.mat', 'Select segDict to use');
    %load(findLatestFilename('_treadBehStruc_'));
    try
    load(segDictName);
    if ~isempty(strfind(segDictName, 'seg2P'))
        C = seg2P.C2p;
    end
    catch
    end
end

if ~exist('treadBehStruc')
    load(findLatestFilename('_treadBehStruc_'));
end

if ~exist('C')
    segDictName = findLatestFilename('_segDict_', 'goodSeg');
    load(segDictName);
    if ~isempty(strfind(segDictName, 'seg2P'))
        C = seg2P.C2p;
    end
end

if ~exist('toZ')
    toZ = 0;
end
    

%load(findLatestFilename('goodSeg'));

% times (and positions) for all cues (slightly diff for rew)
if isempty(strfind(eventName, 'rew'))
evTimes = treadBehStruc.([eventName 'TimeStart']);
evPos = treadBehStruc.([eventName 'PosStart']);
else
    evTimes = treadBehStruc.rewZoneStartTime; %rewTime;
    evPos = treadBehStruc.rewZoneStartPos; %rewPos;
end

%[lapFrInds, lapEpochs] = findLaps(treadBehStruc.resampY(1:2:end));
[lapCueStruc] = findCueLapTypes2(0);
lapTypeArr = lapCueStruc.lapTypeArr;

y = treadBehStruc.resampY; %(1:2:end); % NOTE that lapEpochs are based upon original (non-downsampled) frames
frTimes = treadBehStruc.adjFrTimes; %(1:2:end);

%% Estimate omit times
% if there are omitCue laps, estimate a time for typical cue position each
% lap
cuePos = lapCueStruc.lapTypeCuePos;
lapEpochs = lapCueStruc.lapEpochs;
if min(lapTypeArr)==0
%     cuePos = lapCueStruc.lapTypeCuePos;
%     lapEpochs = lapCueStruc.lapEpochs;
    omitLaps = find(lapTypeArr==0);
    
    for i=1:length(omitLaps)
        try
        epochPos = y(lapEpochs(omitLaps(i),1):lapEpochs(omitLaps(i),2));
        cuePosInd = find(epochPos>cuePos(end),1);
        epochTimes = frTimes(lapEpochs(omitLaps(i),1):lapEpochs(omitLaps(i),2));
        omitCueTimes(i) = epochTimes(cuePosInd);
        catch
        end
    end
    
end

% % times for all cues
% evTimes = treadBehStruc.([eventName 'TimeStart']);

%% find times of cues at different locations, and do eventTrigSig

ca = C(segNum,:);

% zScore Ca if desired
%toZ = 1;
if toZ==1
    ca = zScoreCa(ca);
end

% evTrigSig for each cue type/pos
if length(cuePos)>1 % if multiple cue positions
    cueLapArr = lapTypeArr(find(lapTypeArr~=0)); % laps with cues
    pos1evInd = find(cueLapArr==1);
    pos2evInd = find(cueLapArr==2);
    [evTrigSig1, zeroFr] = eventTrigSig(ca, evTimes(pos1evInd), 0, [-30 120], frTimes(1:2:end));
    [evTrigSig2, zeroFr] = eventTrigSig(ca, evTimes(pos2evInd), 0, [-30 120], frTimes(1:2:end));
    if length(cuePos)>2
        pos3evInd = find(cueLapArr==3);
        [evTrigSig3, zeroFr] = eventTrigSig(ca, evTimes(pos3evInd), 0, [-30 120], frTimes(1:2:end));
        
    end
else % else if only one (probably middle cue)
    
    [evTrigSig2, zeroFr] = eventTrigSig(ca, evTimes, 0, [-30 120], frTimes(1:2:end));
end

try
    [evTrigSig0, zeroFr] = eventTrigSig(ca, omitCueTimes, 0, [-30 120], frTimes(1:2:end));
catch
end


%% save vars to output struc
cueTrigSigStruc.path = pwd;
try
cueTrigSigStruc.segDictName = segDictName;
catch
end
try
cueTrigSigStruc.omitCueSig = evTrigSig0;
catch
    %disp('No omit laps');
end
try
cueTrigSigStruc.shiftCueSig = evTrigSig1;
catch
    %disp('No shift laps');
end
try
cueTrigSigStruc.midCueSig = evTrigSig2;
catch
end

filename = findLatestFilename('.xml');
filename = filename(1:strfind(filename, '.xml')-1);

%% Plotting

if toPlot
figure; 
subplot(2,1,1); hold on;
try
plotMeanSEMshaderr(evTrigSig0, 'r',25:30);
catch
end
try
plotMeanSEMshaderr(evTrigSig1, 'g',25:30);
catch
end
try
    plotMeanSEMshaderr(evTrigSig2, 'b',25:30);
catch
end
try
    plotMeanSEMshaderr(evTrigSig3, 'c',25:30);
catch
end
yl = ylim; xl = xlim;
line([30 30], yl);
ylim(yl); xlim(xl);
title([filename ' ' eventName '-triggered avg., seg=' num2str(segNum)]);
%legend({'omit' 'pos1' 'pos2' 'pos3'});

subplot(2,1,2);
try
plot(evTrigSig0, 'r');
catch
end
hold on;
try
plot(evTrigSig1, 'g');
catch
end
% yl = ylim; xl = xlim;
% line([30 30], yl);
% ylim(yl); xlim(xl);
try
    plot(evTrigSig2, 'b');
catch
end
try
    plot(evTrigSig3, 'c');
catch
end

yl = ylim; xl = xlim;
line([30 30], yl);
ylim(yl); xlim(xl);


end


