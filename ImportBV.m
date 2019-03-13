% pop_loadbv() - load Brain Vision Data Exchange format dataset and
%                return EEGLAB EEG structure
%
% Usage:
%   >> [EEG, com] = pop_loadbv; % pop-up window mode
%   >> [EEG, com] = pop_loadbv(path, hdrfile);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, [], chans);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange, chans);
%
% Optional inputs:
%   path      - path to files
%   hdrfile   - name of Brain Vision vhdr-file (incl. extension)
%   srange    - scalar first sample to read (up to end of file) or
%               vector first and last sample to read (e.g., [7 42];
%               default: all)
%   chans     - vector channels channels to read (e.g., [1:2 4];
%               default: all)
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% Author: Andreas Widmann & Arnaud Delorme, 2004-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: pop_loadbv.m 42 2009-11-12 01:45:54Z arnodelorme $

function [EEG, hdr, com] = ImportBV(path, hdrfile)

com = '';
EEG = [];

if nargin < 2
    [hdrfile path] = uigetfile('*.vhdr', 'Select Brain Vision vhdr-file - pop_loadbv()');
    if hdrfile(1) == 0, return; end
end

% Header file
disp('ImportBV: Reading header file');
hdr = ReadBVConf(path, hdrfile);

% Common Infos
EEG.comments = ['Original file: ' hdr.commoninfos.datafile];
hdr.commoninfos.numberofchannels = str2double(hdr.commoninfos.numberofchannels);
EEG.srate = 1000000 / str2double(hdr.commoninfos.samplinginterval);

% Binary Infos
if strcmpi(hdr.commoninfos.dataformat, 'binary')
    switch lower(hdr.binaryinfos.binaryformat)
        case 'int_16',        binformat = 'int16'; bps = 2;
        case 'uint_16',       binformat = 'uint16'; bps = 2;
        case 'ieee_float_32', binformat = 'float32'; bps = 4;
        otherwise, error('Unsupported binary format');
    end
end

% Channel Infos
chans = 1:hdr.commoninfos.numberofchannels;
EEG.nbchan = hdr.commoninfos.numberofchannels;
    
if isfield(hdr, 'channelinfos')
    for chan = 1:length(chans)
        [EEG.chanlocs(chan).labels, chanlocs(chan).ref, chanlocs(chan).scale, chanlocs(chan).unit] = strread(hdr.channelinfos{chans(chan)}, '%s%s%s%s', 1, 'delimiter', ',');
        EEG.chanlocs(chan).labels = char(EEG.chanlocs(chan).labels);
        chanlocs(chan).scale = str2double(char(chanlocs(chan).scale));
    end
    if isempty([chanlocs.scale])
        chanlocs = rmfield(chanlocs, 'scale');
    end
end

% Coordinates
if isfield(hdr, 'coordinates')
    hdr.coordinates(end+1:length(chans)) = { [] };
    for chan = 1:length(chans)
        if ~isempty(hdr.coordinates{chans(chan)})
            [EEG.chanlocs(chan).sph_radius, theta, phi] = strread(hdr.coordinates{chans(chan)}, '%f%f%f', 'delimiter', ',');
            if EEG.chanlocs(chan).sph_radius == 0 && theta == 0 && phi == 0
                EEG.chanlocs(chan).sph_radius = [];
                EEG.chanlocs(chan).sph_theta = [];
                EEG.chanlocs(chan).sph_phi = [];
            else
                EEG.chanlocs(chan).sph_theta = phi - 90 * sign(theta);
                EEG.chanlocs(chan).sph_phi = -abs(theta) + 90;
            end
        end;
    end
    try,
        [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2topo');
        [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2cart');
    catch, end
end

% Open data file and find the number of data points
% -------------------------------------------------
disp('ImportBV: reading EEG data');
[IN, message] = fopen(fullfile(path, hdr.commoninfos.datafile));
if IN == -1
    [IN, message] = fopen(fullfile(path, lower(hdr.commoninfos.datafile)));
    if IN == -1
        error(message)
    end;
end
if isfield(hdr.commoninfos, 'datapoints')
    hdr.commoninfos.datapoints = str2double(hdr.commoninfos.datapoints);
elseif strcmpi(hdr.commoninfos.dataformat, 'binary')
    fseek(IN, 0, 'eof');
    hdr.commoninfos.datapoints = ftell(IN) / (hdr.commoninfos.numberofchannels * bps);
    fseek(IN, 0, 'bof');
else
    hdr.commoninfos.datapoints = NaN;
end

if ~strcmpi(hdr.commoninfos.dataformat, 'binary') % ASCII
    tmppoint = hdr.commoninfos.datapoints;
    tmpchan = fscanf(IN, '%s', 1);
    tmpdata = fscanf(IN, '%f', inf);
    hdr.commoninfos.datapoints = length(tmpdata);
    chanlabels = 1;
    if str2double(tmpchan) == 0, 
        hdr.commoninfos.datapoints = hdr.commoninfos.datapoints+1; 
        chanlabels = 0;
    end;
    if ~isnan(tmppoint) 
        if tmppoint ~= hdr.commoninfos.datapoints
            error('Error in computing number of data points; try exporting in a different format');
        end;
    end;
end;

% Sample range
srange = [ 1 hdr.commoninfos.datapoints];
EEG.pnts = hdr.commoninfos.datapoints;


% Read data
disp('ImportBV: Reading data file');
if strcmpi(hdr.commoninfos.dataformat, 'binary')
    switch lower(hdr.commoninfos.dataorientation)
        case 'multiplexed'
            if EEG.nbchan == hdr.commoninfos.numberofchannels % Read all channels
                fseek(IN, (srange(1) - 1) * EEG.nbchan * bps, 'bof');
                EEG.data = fread(IN, [EEG.nbchan, EEG.pnts], [binformat '=>float32']);
            else % Read channel subset
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, (srange(1) - 1) * hdr.commoninfos.numberofchannels * bps + (chans(chan) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32'], (hdr.commoninfos.numberofchannels - 1) * bps);
                end
            end
        case 'vectorized'
            if isequal(EEG.pnts, hdr.commoninfos.datapoints) && EEG.nbchan == hdr.commoninfos.numberofchannels % Read entire file
                EEG.data = fread(IN, [EEG.pnts, EEG.nbchan], [binformat '=>float32']).';
            else % Read fraction of file
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, ((chans(chan) - 1) * hdr.commoninfos.datapoints + srange(1) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32']);
                end
            end
        otherwise
            error('Unsupported data orientation')
    end
else % ASCII data
    disp('If this function does not work, export your data in binary format');
    EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]);
    if strcmpi(lower(hdr.commoninfos.dataorientation), 'vectorized')
        count = 1;
        fseek(IN, 0, 'bof');
        len = inf;
        for chan = 1:hdr.commoninfos.numberofchannels
            if chanlabels, tmpchan = fscanf(IN, '%s', 1); end;
            tmpdata = fscanf(IN, '%f', len); len = length(tmpdata);
            if ismember(chan, chans)
                EEG.data(count, :) = tmpdata(srange(1):srange(2))';
                count = count + 1;
            end;
        end;
    elseif strcmpi(lower(hdr.commoninfos.dataorientation), 'multiplexed')
        fclose(IN);
        error('ASCII multiplexed reading not implemeted yet; export as a different format');
    end;
end;

fclose(IN);
EEG.trials = 1;
EEG.xmin   = 0;
EEG.xmax   = (EEG.pnts - 1) / EEG.srate;

% Convert to EEG.data to double for MATLAB < R14
if str2double(version('-release')) < 14
    EEG.data = double(EEG.data);
end

% Scale data
if exist('chanlocs', 'var') && isfield(chanlocs, 'scale')
    disp('ImportBV: Scaling EEG data');
    for chan = 1:EEG.nbchan
        if ~isnan(chanlocs(chan).scale)
            EEG.data(chan, :) = EEG.data(chan, :) * chanlocs(chan).scale;
        end;
    end
end

% Marker file
if isfield(hdr.commoninfos, 'markerfile')
    disp('ImportBV: Reading marker file');
    MRK = ReadBVConf(path, hdr.commoninfos.markerfile);

    if hdr.commoninfos.datafile ~= MRK.commoninfos.datafile
        disp('ImportBV warning: Data files in header and marker files inconsistent.');
    end

    % Marker infos
    if isfield(MRK, 'markerinfos')
        EEG.event = ReadBVMrk(MRK);

        % Correct event latencies by first sample offset
        latency = num2cell([EEG.event(:).latency] - srange(1) + 1);
        [EEG.event(:).latency] = deal(latency{:});

        % Remove unreferenced events
        EEG.event = EEG.event([EEG.event.latency] >= 1 & [EEG.event.latency] <= EEG.pnts);

        % Copy event structure to urevent structure
        EEG.urevent = rmfield(EEG.event, 'urevent');

        % find if boundaries at homogenous intervals
        % ------------------------------------------
        boundaries = strmatch('boundary', {EEG.event.type});
        boundlats = unique([EEG.event(boundaries).latency]);
        if (isfield(hdr.commoninfos, 'segmentationtype') && (strcmpi(hdr.commoninfos.segmentationtype, 'markerbased') || strcmpi(hdr.commoninfos.segmentationtype, 'fixtime'))) && length(boundaries) > 1 && length(unique(diff([boundlats EEG.pnts + 1]))) == 1
            EEG.trials = length(boundlats);
            EEG.pnts   = EEG.pnts / EEG.trials;
            EEG.event(boundaries) = [];

            % adding epoch field
            % ------------------
            for index = 1:length(EEG.event)
                EEG.event(index).epoch = ceil(EEG.event(index).latency / EEG.pnts);
            end

            % finding minimum time
            % --------------------
            tles = strmatch('time 0', lower({EEG.event.code}))';
            if ~isempty(tles)
                [EEG.event(tles).type] = deal('TLE');
                EEG.xmin = -(EEG.event(tles(1)).latency - 1) / EEG.srate;
            end
        end
    end
end

EEG.ref = 'common';

if nargout == 3
    com = sprintf('EEG = pop_loadbv(''%s'', ''%s'', %s, %s);', path, hdrfile, mat2str(srange), mat2str(chans));
end
