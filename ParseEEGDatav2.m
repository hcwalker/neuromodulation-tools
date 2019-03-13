%%% Parse and format data from a *.eeg file
%%%     - uses a modified "pop_loadbv.m" by Andreas Widmann & Arnaud Delorme
%%%     - splits the .eeg file into 2 separate files:
%%%         (1) Trial file: contains formatted stimulus events & experiment info in an E struct
%%%         (2) Data file: contains continuous data in a matlab matrix (chans x samples)
%%%
%%% Input:
%%%     - fname: name (including path) of the BrainVision header file (*.vhdr)
%%%         - if not included, function will pop up file picker
%%%     - paradigm: string indicating if the experiment used single or paired pulses
%%%         - options include 'single' or 'paired'
%%%     - ttlchan: which TTL channel to sync to (can be 1 or 2)
%%%
%%% Output:
%%%     - E: the extracted structure array of trial information
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = ParseEEGDatav2(fname, paradigm, ttlchan)

% Pop up GUI file picker if file name not included:
if (nargin < 1 || isempty(fname))
    [eegfname, eegfpath] = uigetfile('*.vhdr', 'Select an EEG file (*.vhdr)');
else
    [eegfpath, eegfname, ext] = fileparts(fname);
    eegfname = [eegfname, ext];
end
if (nargin < 2 || isempty(paradigm))
    error('Please input an experimental paradigm: either ''single'' or ''paired''');
end
if (nargin < 3 || isempty(ttlchan))
    ttlchan = 1;
end

% Convert BrainVision format to EEGlab format (includes data) and save:
[EEG, hdr] = ImportBV(eegfpath, eegfname);
% save(fullfile(eegfpath, [eegfname(1:strfind(eegfname, '.')-1), '_EEG.mat']), 'EEG');

% Check if data has previously been saved:
datafname = fullfile(eegfpath, [eegfname(1:strfind(eegfname, '.')-1), '_Data.mat']);
saveflag = true;
if (exist(datafname, 'file') == 2)
    S = whos('-file', datafname);
    if (strcmp(S(1).name, 'Data') && all(size(EEG.data) == S.size))
        saveflag = false;
        fprintf('Data file already exists: %s\n', datafname);
    end
end     

% Save data into a separate file if not already and clear out of memory:
if (saveflag)
    fprintf('Saving data file to: %s\n', datafname);
    datafile = matfile(datafname, 'Writable', true);
    datafile.Data = EEG.data;
end
EEG = rmfield(EEG, 'data');

% Reformat EEG struct (make the events more usable):

%% Select and import files:

% User selects corresponding MC_Stimulus file:
[stimfname, stimfpath] = uigetfile(fullfile(eegfpath, '*.txt'), 'Select corresponding Stimulus file'); 

% Load in file:
stim = fileread(fullfile(stimfpath, stimfname));

% Make sure it's a correct MC_Stimulus file:
assert(strcmpi(stim(1:36), 'Multi Channel Systems MC_Stimulus II'), 'File error: Not a valid MC_Stimulus file');

% Extract the name of the original *.eeg file from the EEG struct:
eegfile = EEG.comments(16:end);

%% Convert stimulus file to correct format:

% Extract the output mode (voltage or current) and corresponding units:
idx = strfind(stim, ['output mode:', char(9)]);
outputmode = stim(idx+13:idx+19);
if (strcmpi(outputmode, 'voltage'))
    units = 'mV';
else
    units = 'uA';
end

% Determine maximum number of channels for this stimulator:
maxchans = str2double(stim(strfind(stim, 'output mode:')-3));

% Find TTL section and convert to numbers:

%%%
%%% extract the desired TTL channel
%%% sync to whichever channel(s) match its length/ttl widths

ttlstart = strfind(stim, ['channel:', char(9), num2str(ttlchan+maxchans)]) + 53;
ttlvals = sscanf(stim(ttlstart:end), '%f');
assert(~any(ttlvals(repmat(logical([0 0 1 1 1 0 1]'),length(ttlvals)/7,1))), 'Format error: TTL columns 3, 4, 5, and 7 are not empty');
assert(all(ttlvals(1:7:end) == 1), 'Format error: TTL voltages not all 1');

stimvals = cell(1,maxchans); usedchans = zeros(1,maxchans,'logical');
for i = 1:maxchans
    
    % Find start and end of ith channel stimuli definitions:
    stimstart = strfind(stim, sprintf('channel:\t%d\r\n', i)) + 53;
    stimend = strfind(stim, sprintf('channel:\t%d\r\n', i+1)) - 5;
    
    % Convert this section to numbers:
    stimvals{i} = sscanf(stim(stimstart:stimend), '%f');
    
    % Make sure this channel is used:
    if (isempty(stimvals{i}) || length(stimvals{i}) ~= length(ttlvals))
        stimvals{i} = [];
        continue;
    end
    
    % Make sure the file matches the format we think it's in:
%     assert(length(stimvals{i}) == length(ttlvals), ['Format error channel ', num2str(i), ': Lengths of stimuli and TTL don''t match']);
    assert(mod(length(stimvals{i}),7) == 0, 'Format error channel ', num2str(i), ': Number of columns is wrong');
    assert(~any(stimvals{i}(repmat(logical([0 0 0 0 1 0 1]'),length(stimvals{i})/7,1))), 'Format error channel ', num2str(i), ': Stimulus columns 5 and 7 are not empty');
    
    % Check if this channel is used at all:
    if (any(stimvals{i}(repmat(logical([1 0 1 0 0 0 0]'),length(stimvals{i})/7,1))))
        usedchans(i) = true;
    else
        continue;
    end
    
    % Check if this is a monophasic or biphasic stim file:
    if (~any(stimvals{i}(3:7:end))) % Only 1 stim pulse per line = monophasic stimulation
        isMono = true;
        stimvals{i}(6:7:end) = stimvals{i}(6:7:end) + stimvals{i}(4:7:end); % add any time from the empty 2nd pulse to the ISI time
    else
        isMono = false;
    end
    
    % Get rid of the unused columns in each vector:
    %   * remaining columns: 1 = pulse amplitude, 2 = pulse width, 3 = inter-pulse delay
    stimvals{i} = stimvals{i}(repmat(logical([1 1 ~isMono ~isMono 0 1 0]'),length(stimvals{i})/7,1));
    
    % Check if 
    
end

% Get rid of unused columns in the TTL vector:
ttlvals = ttlvals(repmat(logical([1 1 0 0 0 1 0]'),length(ttlvals)/7,1));

% Determine which channels were actually used:
chanidxs = find(usedchans);
% chanidxs = find(cellfun(@(x) ~isempty(x), stimvals));

%% Extract and measure TTL detection events from the EEG file recording:
       
% Find TTL edge detection events in the EEG file:
events = cellfun(@(x) str2double(x(end-1:end)), {EEG.event(2:end).type}); % all marker numbers
eventsbin = de2bi(events); % convert to binary
trigidxs = find(~all(eventsbin,1)); % which trigger outputs were used for the TTL channel(s)
trig = trigidxs(ttlchan); % which trigger are we syncing to
risingedgeidx = find(diff([0;eventsbin(:,trig)]) > 0) + 1;
fallingedgeidx = find(diff(eventsbin(:,trig)) < 0) + 2;
fprintf('    - Using trigger %d to sync\n', trig);





% risingedgeidx = find(strcmp({EEG.event.type}, EEG.event(2).type));
% fallingedgeidx = find(strcmp({EEG.event.type}, EEG.event(3).type));
% disp(['    - Rising edge event label: ', EEG.event(2).type]);
% disp(['    - Falling edge event label: ', EEG.event(3).type]);
% risingedgeidx = find(strcmp({EEG.event.type}, 'S  7'));
% fallingedgeidx = find(strcmp({EEG.event.type}, 'S  3'));

% Make sure TTL edge detections aren't screwed up:
if (length(fallingedgeidx) == length(risingedgeidx))
    assert(all(fallingedgeidx-risingedgeidx >= 1), 'Event error: Rising and falling edges not adjacent');
elseif (length(fallingedgeidx) == length(risingedgeidx)-1)
    assert(all(fallingedgeidx-risingedgeidx(1:end-1) >= 1), 'Event error: Rising and falling edges not adjacent');
    risingedgeidx = risingedgeidx(1:end-1);
else
    error('Event error: Missing half of a TTL pair');
end


% assert(length(fallingedgeidx) == length(risingedgeidx), 'Event error: Missing a TTL edge');
% if (isMono)
%     assert(length(fallingedgeidx) == length(stimvals{1})/3, 'Event error: Number of detected events doesn''t match stimulus file');
% else
%     assert(length(fallingedgeidx) == length(stimvals{1})/5, 'Event error: Number of detected events doesn''t match stimulus file');
% end
% assert(all(fallingedgeidx-risingedgeidx == 1), 'Event error: Rising and falling edges not adjacent');
% assert(mod(length(fallingedgeidx), 2) == 0, 'Event error: Missing half of a pair');

% Measure width of each TTL pulse:
risingsample = [EEG.event(risingedgeidx).latency];
fallingsample = [EEG.event(fallingedgeidx).latency];
widths = fallingsample - risingsample;

% Get intended TTL widths:
ttlwidths = ttlvals(2:3:end)'*EEG.srate/1e6; % in samples

% Find where intended TTL matches up with our measured values:
for i = 1:length(ttlwidths)
    if (all(abs(ttlwidths(i + (1:length(widths)) - 1) - widths) < 2))
        break;
    end
end

% If we didn't find a match, or there weren't enough captured values, throw an error:
if (i == length(ttlwidths))
    error('Event error: measured TTL widths don''t match stimulus file');
else
    eventstart = i;
end

% % Make sure widths match up with intended pulses:
% ttlvals = ttlvals(1:length(risingedgeidx)*3);
% assert(all(abs(ttlvals(2:3:end)'*EEG.srate/1e6 - widths) < 2), 'Event error: measured TTL widths don''t match stimulus file');

%% Create and fill in new data structure:
%   - NOTE: all time units in milliseconds, except where noted by "..._sample"

% Create new struct to hold data:
E = struct('Info', struct('SampleRate', EEG.srate, 'NumDataPts', EEG.pnts, ...
                          'NumChans', EEG.nbchan, 'ChanLabels', {{EEG.chanlocs.labels}}, ...
                          'StimulatorOutputChannels', chanidxs, ...
                          'StimulatorOutputMode', outputmode, 'StimulatorOutputUnits', units, ...
                          'TimeUnits', 'ms', 'IsMonophasic', isMono, ...
                          'ExperimentParadigm', paradigm, ... 
                          'EEGFileName', eegfile, 'StimFileName', stimfname, ...
                          'SyncTriggerIdx', trig), ...
           'Trials', []);
       
% Number of trials:
n = length(widths);
       
% Split events into trials (paired pulses):
E.Trials = repmat(struct('TrialNumber', [], 'EEGEventIdxs', [], 'StimEventIdx', [], ...
                         'StimLoc_sample', [],  ...
                         'StimWidth', [], ...
                         'StimAmp', [], ...
                         'TTLWidth', [],  ...
                         'NominalPrevISI', [], 'MeasuredPrevISI', [], ...
                         'NominalNextISI', [], 'MeasuredNextISI', []), 1, n);

for i = 1:n
    
    % Record the trial:
    E.Trials(i).TrialNumber = i;
    
    % Which EEG struct/Stim file events does this trial/pair contain:
%     E.Trials(i).EEGEventIdxs = (1:2)+(i-1)*2 + 1;  % [TTL up, TTL down]
    E.Trials(i).EEGEventIdxs = [risingedgeidx(i), fallingedgeidx(i)];  % [TTL up, TTL down]
    E.Trials(i).StimEventIdx = i;     

    % Record width (in ms) of TTL pulse:
    E.Trials(i).TTLWidth = widths(i)*1e3/EEG.srate;
    
    % Record location (in samples) of pulse:
    E.Trials(i).StimLoc_sample = EEG.event(risingedgeidx(i)).latency;  % Pulse TTL up
    
    % Record width (in ms) of stimulus pulse:
    if (isMono)
        E.Trials(i).StimWidth = arrayfun(@(x) stimvals{x}(2 + (i+eventstart-2)*3)/1000, chanidxs);
    else % [channel1phase1 ... channelNphase1; channel1phase2 ... channelNphase2]
        E.Trials(i).StimWidth = cell2mat(arrayfun(@(x) stimvals{x}([2,4] + (i+eventstart-2)*5)/1000, chanidxs, 'UniformOutput', false));
    end
    
    % Record amplitude (in either mV or uA) of stimulus pulse:
    if (isMono)
        E.Trials(i).StimAmp = arrayfun(@(x) stimvals{x}(1 + (i+eventstart-2)*3), chanidxs);
    else % [channel1phase1 ... channelNphase1; channel1phase2 ... channelNphase2]
        E.Trials(i).StimAmp = cell2mat(arrayfun(@(x) stimvals{x}([1,3] + (i+eventstart-2)*5), chanidxs, 'UniformOutput', false));
    end
    
    % Record the inter-stimulus intervals (in ms) for this pulse:
    %   NOTE: PrevISI measured between onset of previous pulse and
    %         onset of current pulse. NextISI measured between onset
    %         of current pulse and onset of next pulse.
    %   NOTE: Nominal ISI = intended ISI (from stim file). 
    %         Measured ISI = actual ISI (from recorded TTL events).
    if (i > 1)
        E.Trials(i).MeasuredPrevISI = (EEG.event(risingedgeidx(i)).latency - EEG.event(risingedgeidx(i-1)).latency)*1e3/EEG.srate;
        if (isMono)
            E.Trials(i).NominalPrevISI = arrayfun(@(x) sum(stimvals{x}([2,3] + (i-1+eventstart-2)*3))/1000, chanidxs);
        else
            E.Trials(i).NominalPrevISI = arrayfun(@(x) sum(stimvals{x}([2,4,5] + (i-1+eventstart-2)*5))/1000, chanidxs);
        end
        if (all(E.Trials(i).NominalPrevISI == E.Trials(i).NominalPrevISI(1)))
            E.Trials(i).NominalPrevISI = E.Trials(i).NominalPrevISI(1);
        else
            warning(['Channel PrevISIs are not equal - check trial ', num2str(i)]);
        end
        assert(abs(E.Trials(i).MeasuredPrevISI - E.Trials(i).NominalPrevISI) < .06, ['Trial ', num2str(i), ' error: Nominal and measured previous ISI mismatch']);
    end
    if (i < n)
        E.Trials(i).MeasuredNextISI = (EEG.event(risingedgeidx(i+1)).latency - EEG.event(risingedgeidx(i)).latency)*1e3/EEG.srate;
        if (isMono)
            E.Trials(i).NominalNextISI = arrayfun(@(x) sum(stimvals{x}([2,3] + (i+eventstart-2)*3))/1000, chanidxs);
        else
            E.Trials(i).NominalNextISI = arrayfun(@(x) sum(stimvals{x}([2,4,5] + (i+eventstart-2)*5))/1000, chanidxs);
        end
        if (all(E.Trials(i).NominalNextISI == E.Trials(i).NominalNextISI(1)))
            E.Trials(i).NominalNextISI = E.Trials(i).NominalNextISI(1);
        else
            warning(['Channel NextISIs are not equal - check trial ', num2str(i)]);
        end
        assert(abs(E.Trials(i).MeasuredNextISI - E.Trials(i).NominalNextISI) < .06, ['Trial ', num2str(i), ' error: Nominal and measured next ISI mismatch']);
    end
    
end

trialfile = matfile(fullfile(eegfpath, [eegfname(1:strfind(eegfname, '.')-1), '_Trials.mat']), 'Writable', true);
trialfile.E = E;
trialfile.HDR = hdr;