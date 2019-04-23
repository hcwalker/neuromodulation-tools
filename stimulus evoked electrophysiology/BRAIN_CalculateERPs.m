%%%% Calculate ERPs for all conditions (one condition per file):

%% user settings:

subjectID = 'BI004';

electrodeIDs = [1, 2, 3, 4, 5, 6, 13];
chanprefix = 'Ch'; 

chans = [1, 3, 7];
refs =  [2, 4, 6];

amplabels = [1, 2, 3, 4, 5];

filtCutoffs = [50, 1000];
dataWin = [-30, 70];
avgAll = true; % average all channels together after computing ERP (useful for EEG)
phaseRectify = [1,1,1];  % multiplier for phase rectification during averaging
binSize = 0;    % size of bin to smooth data in ms (useful for EMG)

%% load files and calculate ERPs:

[eegfname, eegfpath] = uigetfile('*.vhdr', 'Select an EEG file (*.vhdr)', 'MultiSelect', 'on');
if ~iscell(eegfname)
    eegfname = {eegfname};
end

% convert channel numbers to electrode IDs:
chanIds = electrodeIDs(chans);
refIds = zeros(1,length(refs));
for i = 1:length(refIds)
    if (refs(i) ~= 0)
        refIds(i) = electrodeIDs(refs(i));
    end
end

% collect trials and data from each file:
stimchannames = cell(1,length(eegfname));
sampleRate = 0;
condERPs = zeros(length(eegfname), length(refs), (dataWin(2)-dataWin(1))*100 + 1, 5);
    % stim channel, recording channel, samples, amplitude conditions
for i = 1:length(eegfname)
    
    disp(['-- Loading file ', num2str(i), ' out of ', num2str(length(eegfname))]);
    
    % load data:
    if (~exist(fullfile(eegfpath, [eegfname{i}(1:end-5), '_Trials.mat']), 'file'))
        E = ParseEEGDatav2(fullfile(eegfpath, eegfname{i}), 'single');
    else
        load(fullfile(eegfpath, [eegfname{i}(1:end-5), '_Trials.mat']), 'E');
    end
    datafile = fullfile(eegfpath, [eegfname{i}(1:end-5), '_Data.mat']);
    
    % find stim channel number:
    isnum = isstrprop(eegfname{i}, 'digit');
    first = find([0, isnum], 1)-1; last = find([~isnum(first+1:end), 1] > 0, 1) + first - 1;
    stimchannames{i} = eegfname{i}(first:last);
    
    % extract data windows around both pulses and add to E struct:
    E = AddDataWindowsToStruct(E, datafile, dataWin, chanprefix, chanIds, refIds, filtCutoffs, 1, binSize, false);
    
    % collect together:
    trials = E.Trials;
    
    % identify sampling rate:
    sampleRate = E.Info.SampleRate;
    
    % clear original struct from memory:
    clear E;
    
    % compute ERPs for each amplitude:
    amps = unique(arrayfun(@(x) abs(x.StimAmp(1)), trials));    
    for j = 1:5
        
        % get trials for this amplitude:
        trialsamp = trials(abs(arrayfun(@(x) x.StimAmp(1), trials)) == amps(j));
        
        % gather data windows for these trials:
        erpdata = zeros(length(refs), size(trials(1).DataWindow.DataFirstPulse, 2), length(trialsamp));
        erppos = false(1,length(trialsamp));
        for k = 1:length(trialsamp)
            
            erpdata(:,:,k) = trialsamp(k).DataWindow.DataFirstPulse;
            erppos(k) = trialsamp(k).StimAmp(1) > 0;
            
        end
        
        % remove outliers:
        maxamp = max(abs(erpdata),[],2);
        outliers = maxamp > repmat(3*prctile(maxamp, 90, 3), 1, 1, size(maxamp, 3));
        erpdata(repmat(outliers,1,size(erpdata,2),1)) = NaN;
        numoutliers = sum(sum(sum(outliers)));
        disp(['    -- Stim channel ', stimchannames{i}, ', amplitude ', num2str(j), ...
            ': removed ', num2str(numoutliers), ' total outliers (', num2str(100*numoutliers/size(erpdata,3),2), '%)'])
        
        % compute condition ERPs:
        condERPs(i,:,:,j) = mean(erpdata(:,:,erppos), 3, 'omitnan')+mean(erpdata(:,:,~erppos), 3, 'omitnan');
        
    end
    
end

% average channels together if requested:
if (avgAll)
    condERPs = mean(condERPs.*repmat(phaseRectify, size(condERPs,1), 1, size(condERPs,3), size(condERPs,4)), 2);
    chans = 0; refs = 0;
end

%% create struct for data and save:

info = struct();
info.Subject = subjectID;
info.Files = eegfname;
info.FilePath = eegfpath;
info.FilterCutoffs_Hz = filtCutoffs;
info.DataWindow_ms = dataWin;
info.SamplingRate_Hz = sampleRate;
if (~isempty(chans))
    info.ECoGChans = chans;
else
    info.ECoGChans = 1:length(refs);
end
info.ECoGRefs = refs;
info.ECoGBrainAreas = cell(size(refs));
info.Amplitudes_mA = amplabels;

peaks = repmat(struct('ECoG', ...
        repmat(struct('Baseline', [], 'Latency', [], 'Amplitude', [], 'Group', []), ...
        1, size(condERPs,2))), 1, size(condERPs,1));

Data = struct('Info', info, 'ERPs', condERPs, 'Peaks', peaks);

[file, fpath] = uiputfile('*.mat');

save(fullfile(fpath, file), 'Data');
