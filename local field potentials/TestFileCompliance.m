clc
close all
clearvars

%% Open header file
[fileNameVHDR, pathName] = uigetfile('*.vhdr','Select Header File');
fid = fopen(fullfile(pathName,fileNameVHDR),'r');

%% Find NumberOfChannels value
while true
    str = fgetl(fid);
    if contains(str,'NumberOfChannels')
        break
    end
end
NumberOfChannels = str2double(str(18:end));

%% Find SamplingInterval value
while true
    str = fgetl(fid);
    if contains(str,'SamplingInterval')
        break
    end
end
SamplingInterval = str2double(str(18:end));
Fs = 1000000/SamplingInterval;
if Fs~=25000
    warning('The sampling frequency is not 25KHz')
end

%% find available channel names
while true
    str = fgetl(fid);
    if contains(str,'Ch1')
        break
    end
end

%% get LFP channel numbers
chName = cell(NumberOfChannels,1);
LFPchannels = [];

for chCnt = 1:NumberOfChannels
    tmp = textscan(str,'%s','Delimiter',',');
    tmp = tmp{1,1}{1,1};
    tmp = textscan(tmp,'%s','Delimiter','=');
    chName{chCnt} = tmp{1,1}{2,1};
    if contains(chName{chCnt},'LFP')
        LFPchannels = [LFPchannels;chCnt];
    end
    str = fgetl(fid);
end

%% show impedances out of range
while true
    str = fgetl(fid);
    k = strfind(str,'Impedance');
    if ~isempty(k)
        break
    end
end

cnt = 1;
ImpedanceError = [];
while ~feof(fid)
    str = fgetl(fid);
    k = strfind(str,'Out of Range');
    if ~isempty(k)
        ImpedanceError{cnt} = str(4:end);
        cnt = cnt+1;
    end
end
if ~isempty(ImpedanceError)
    disp(ImpedanceError)
end
fclose(fid); % Close the text file