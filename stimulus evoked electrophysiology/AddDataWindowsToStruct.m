%%% Function to add data windows to each trial of an E Struct.
%%%
%%% Inputs:
%%%     * E - E struct from a "*_Trials.mat" file
%%%     * dataFile - path to data file for this E struct
%%%         * NOTE: if current working directory is in same folder with
%%%                 data file, this is not needed
%%%     * window - [1x2] length of window [pre, post]-stimulus in ms
%%%     * chans - [1xN] electrode numbers of channels to extra windows for
%%%     * refs - [1xN] electrode numbers of channels to reference each chan to
%%%     * filtCutoffs - [1x2] bandpass filter cutoffs [lo, hi] in Hz
%%%     * pulseSelect - selector for which pulse to grab the window around
%%%         * 1 = 1st pulse
%%%         * 2 = 2nd pulse
%%%         * 3 = both pulses (adds two data fields to the struct
%%%
%%% Output:
%%%     * EWithData - same E struct as input, with added fields for data windows
%%%         * E.Trials(x).DataWindow
%%%             * Window
%%%             * ChanIdxs
%%%             * RefIdxs
%%%             * FiltCutoffs
%%%             * Pulse
%%%             * DataFirstPulse - [NxM] N channels, M data points
%%%             * DataSecondPulse - [NxM] N channels, M data points
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EWithData = AddDataWindowsToStruct(E, dataFile, window, chanPrefix, chans, refs, filtCutoffs, pulseSelect, binSize, saveFlag)

% Input defaults and error checking:
if (nargin < 9 || isempty(binSize))
    binSize = 0;
end
if (nargin < 10 || isempty(saveFlag))
    saveFlag = false;
end
if (~ismember(pulseSelect, 1:3))
    error('Please enter a valid pulse selection');
end
if (isempty(chans))
    chans = cell(1,length(chanPrefix));
    refs = num2cell(refs);
else
    chans = num2cell(chans);
    refs = num2cell(refs);
end

[bhi, ahi] = butter(2, filtCutoffs(1)/(E.Info.SampleRate/2), 'high');
[blo, alo] = butter(2, filtCutoffs(2)/(E.Info.SampleRate/2), 'low');

% Find data file for this struct:
if (isempty(dataFile))
    try
        dataFile = [E.Info.EEGFileName(1:end-4), '_Data.mat'];
    catch
        error('Could not find file, either provide a path or move current directory');
    end
end

if (iscell(chanPrefix))
    chanIdxs = cellfun(@(x,y) find(strcmp([y, num2str(x)], E.Info.ChanLabels)), chans, chanPrefix);
    refIdxs = cellfun(@(x,y) find(strcmp([y, num2str(x)], E.Info.ChanLabels)), refs, chanPrefix, 'UniformOutput', false);
else
    chanIdxs = cellfun(@(x) find(strcmp([chanPrefix,num2str(x)], E.Info.ChanLabels)), chans);
    refIdxs = cellfun(@(x) find(strcmp([chanPrefix,num2str(x)], E.Info.ChanLabels)), refs, 'UniformOutput', false);
end

% if (~isnan(str2double(E.Info.ChanLabels{1})))
%     chanIdxs = arrayfun(@(x) find(strcmp(num2str(x), E.Info.ChanLabels)), chans);
%     refIdxs = arrayfun(@(x) find(strcmp(num2str(x), E.Info.ChanLabels)), refs, 'UniformOutput', false);
% else
%     chanIdxs = arrayfun(@(x) find(strcmp(['Ch',num2str(x)], E.Info.ChanLabels)), chans);
%     refIdxs = arrayfun(@(x) find(strcmp(['Ch',num2str(x)], E.Info.ChanLabels)), refs, 'UniformOutput', false);
% %     chanIdxs = arrayfun(@(x) find(strcmp(['ECOG',num2str(x)], E.Info.ChanLabels)), chans);
% %     refIdxs = arrayfun(@(x) find(strcmp(['ECOG',num2str(x)], E.Info.ChanLabels)), refs, 'UniformOutput', false);
% end

% Load and filter data:
Data1 = load(dataFile);

%%% make this faster by using matfile() to only load req. chans

Data = zeros(length(chanIdxs), size(Data1.Data,2));
for i = 1:length(chanIdxs)
    if (~isempty(refIdxs{i}))
        Data(i,:) = filtfilt(blo, alo, filtfilt(bhi, ahi, double(Data1.Data(chanIdxs(i),:))) - filtfilt(bhi, ahi, double(Data1.Data(refIdxs{i}, :))));
    else
        Data(i,:) = filtfilt(blo, alo, filtfilt(bhi, ahi, double(Data1.Data(chanIdxs(i),:))));
    end
end
clear Data1;

% Moving (absolute) average, if requested:
if (binSize > 0)
    
    % Convert bin size in ms to samples (make sure it's odd):
    binsamples = binSize*E.Info.SampleRate/1000;
    binsamples = binsamples - (1 - mod(binsamples,2));
    
    % Rectify and smooth data:
    for i = 1:size(Data,1)
        Data(i,:) = smooth(abs(Data(i,:)'), binsamples)';
        Data(i,:) = Data(i,:) - mean(Data(i,:),2);
    end
end

% Extract windows for each trial:
pulsestrings = {'First', 'Second', 'Both'};
for i = 1:length(E.Trials)
    
    % Set info fields of DataWindow struct:
    E.Trials(i).DataWindow.Window = window;
    E.Trials(i).DataWindow.ChanIdxs = chanIdxs;
    E.Trials(i).DataWindow.RefIdxs = cellfun(@(x) max([length(x), x]), refIdxs);
    E.Trials(i).DataWindow.FiltCutoffs = filtCutoffs;
    E.Trials(i).DataWindow.Pulse = pulsestrings{pulseSelect};
    
    % Extract windows:
    if (pulseSelect == 1 || pulseSelect == 3)
        if (isfield(E.Trials(i), 'FirstStimLoc_sample'))
            t = (E.Trials(i).FirstStimLoc_sample) + ((window(1)*E.Info.SampleRate/1000):(window(2)*E.Info.SampleRate/1000));
        else
            t = (E.Trials(i).StimLoc_sample) + ((window(1)*E.Info.SampleRate/1000):(window(2)*E.Info.SampleRate/1000));
        end
        if (t(1) > 0)
            E.Trials(i).DataWindow.DataFirstPulse = Data(:,t);
        end
    end
    if ((pulseSelect == 2 || pulseSelect == 3) && isfield(E.Trials(i), 'SecondStimLoc_sample'))
        t = (E.Trials(i).SecondStimLoc_sample) + ((window(1)*E.Info.SampleRate/1000):(window(2)*E.Info.SampleRate/1000));
        if (t(1) > 0)
            E.Trials(i).DataWindow.DataSecondPulse = Data(:,t);
        end
    end
        
end

% Save (if asked) and output:
if (saveFlag)
    save([dataFile(1:end-8), '_Trials_withWindows.mat'], 'E');
end
EWithData = E;
