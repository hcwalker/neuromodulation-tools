function varargout = ERP_Explorer(varargin)
% ERP_EXPLORER MATLAB code for ERP_Explorer.fig
%      ERP_EXPLORER, by itself, creates a new ERP_EXPLORER or raises the existing
%      singleton*.
%
%      H = ERP_EXPLORER returns the handle to a new ERP_EXPLORER or the handle to
%      the existing singleton*.
%
%      ERP_EXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ERP_EXPLORER.M with the given input arguments.
%
%      ERP_EXPLORER('Property','Value',...) creates a new ERP_EXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ERP_Explorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ERP_Explorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ERP_Explorer

% Last Modified by GUIDE v2.5 06-Jun-2018 09:35:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ERP_Explorer_OpeningFcn, ...
                   'gui_OutputFcn',  @ERP_Explorer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ERP_Explorer is made visible.
function ERP_Explorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ERP_Explorer (see VARARGIN)

% Choose default command line output for ERP_Explorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ERP_Explorer wait for user response (see UIRESUME)
% uiwait(handles.GUIFigure);


% --- Outputs from this function are returned to the command line.
function varargout = ERP_Explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load data structure:
[fname, fpath] = uigetfile('*.mat', 'Select a BRAIN data structure file (*.mat)', 'MultiSelect', 'off');
load(fullfile(fpath, fname), 'Data');

% extract some information:
dw = Data.Info.DataWindow_ms;
fc = Data.Info.FilterCutoffs_Hz;

% populate some data fields:
handles.DataPanel.Title = ['Subject: ', Data.Info.Subject, ', File: ', fullfile(fpath, fname)];
handles.FilterCutoffsText.String = sprintf('Filter Cutoffs: [%d, %d]', fc);
handles.DataWindowText.String = sprintf('Data Window: [%d, %d]', dw);
handles.BrainAreaTable.Data = Data.Info.ECoGBrainAreas;
handles.BrainAreaTable.Position(3) = 15*length(Data.Info.ECoGChans)+.4;
if (iscell(Data.Info.ECoGChans))
    handles.BrainAreaTable.ColumnName = Data.Info.ECoGChans;
else
    handles.BrainAreaTable.ColumnName = arrayfun(@(x,y) sprintf('%d-%d', x, y), Data.Info.ECoGChans, Data.Info.ECoGRefs, 'UniformOutput', false);
end
handles.FilePath = fpath;
handles.FileName = fname;

% add data to handles struct to pass it around:
handles.Data = Data;

% create default baseline/peak data if not already there:
for i = 1:length(handles.Data.Peaks)
    for j = 1:length(handles.Data.Peaks(i).ECoG)
        if (isempty(handles.Data.Peaks(i).ECoG(j).Baseline))
            handles.Data.Peaks(i).ECoG(j).Baseline = zeros(2, length(Data.Info.Amplitudes_mA));
        end
        if (isempty(handles.Data.Peaks(i).ECoG(j).Latency))
            handles.Data.Peaks(i).ECoG(j).Latency = NaN(1, length(Data.Info.Amplitudes_mA));
            handles.Data.Peaks(i).ECoG(j).Amplitude = NaN(1, length(Data.Info.Amplitudes_mA));
            handles.Data.Peaks(i).ECoG(j).Group = NaN;
        end
    end
end

% reset plot display check boxes to checked (i.e. display everything):
for i = 1:10
    handles.(['checkbox',num2str(i)]).Value = 1;
end

% plot current view of data:
handles = UpdatePlot(handles);

% update table display:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);

% --- Executes on button press on Line object.
function LineButtonDown_Callback(hObject, eventdata)%, handles)

% get line information:
x = hObject.XData;
y = hObject.YData;
if (isempty(hObject.UserData.BaselineMarker))
    basey = 0;
else
    basey = hObject.UserData.BaselineMarker.YData;
end
% basey = hObject.

% get current handles object:
handles = guidata(ancestor(hObject, 'figure'));

% get parent axis handles:
ax = handles.DataAxes;

% get current point clicked (on the line):
[~, idx] = min(sqrt((x-ax.CurrentPoint(1,1)).^2 + (y-ax.CurrentPoint(1,2)).^2));
point = [x(idx), y(idx)];

% if there's a peak closer than .2ms our point, snap to it:
[dist, idx] = min(abs(hObject.UserData.AllPeaks - point(1)));
if (dist < .2)
    point = [hObject.UserData.AllPeaks(idx), y(x == hObject.UserData.AllPeaks(idx))];
end

% display point coordinates:
handles.XPointText.String = sprintf('X: %.2f', point(1));
handles.YPointText.String = sprintf('Y: %.2f (%.2f)', point(2), point(2) - basey);

% turn off the current point marker by default:
handles.CurrentPointMarker.Visible = 'off';

% determine any additional actions:
switch (handles.CurrentActionButtonGroup.SelectedObject.String)
    case ('No mark')                % display generic mark (just one allowed for all lines)
        handles.CurrentPointMarker.XData = point(1);
        handles.CurrentPointMarker.YData = point(2);
        handles.CurrentPointMarker.Visible = 'on';
        
    case ('Mark baseline')    % mark a baseline point for this line (just one allowed per line)
        if (isempty(hObject.UserData.BaselineMarker))
            handles = UpdateBaseline(hObject, point, handles, 'add');
        end
        
    case ('Delete baseline') % mark the baseline point for this line if there is one
        if (~isempty(hObject.UserData.BaselineMarker))
            handles = UpdateBaseline(hObject, point, handles, 'delete');
        end
        
    case ('Mark peak')        % mark a peak point for this line (multiple peaks allowed per line)
        mrks = hObject.UserData.PeakMarkers;
        if (isempty(mrks) || ~any(abs([mrks.XData] - point(1)) < .2))
            handles = UpdatePeaks(hObject, point, handles, 'add');
        end
        
    case ('Delete peak')           % delete the selected peak point
        % look through peak markers for one clicked:
        mrks = hObject.UserData.PeakMarkers;
        if (~isempty(mrks))
            [dist, idx] = min(abs([mrks.XData] - point(1)));
            if (dist < .1)
                handles = UpdatePeaks(hObject, mrks(idx), handles, 'delete');
            end
        end
        
    case ('Group points')           % group together multiple points
        % set figure callbacks to detect mouse motion and button release:
        handles.GUIFigure.WindowButtonMotionFcn = @GUIFigure_WindowButtonMotionFcn;
        handles.GUIFigure.WindowButtonUpFcn = @GUIFigure_WindowButtonUpFcn;
        
        % initialize a selection rectangle:
        x = handles.DataAxes.CurrentPoint(1,1); y = handles.DataAxes.CurrentPoint(1,2);
        handles.SelectionRect = line(handles.DataAxes, [x, x, x, x, x], [y, y, y, y, y]);
        
end

guidata(handles.GUIFigure, handles);


% --- Executes on button press in SaveDataButton.
function SaveDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname, fpath] = uiputfile(fullfile(handles.FilePath, '*.mat'), ...
    'Select a BRAIN data structure file (*.mat)', handles.FileName);

handles.FilePath = fpath;
handles.FileName = fname;

Data = handles.Data;    %#ok
save(fullfile(fpath, fname), 'Data');

guidata(handles.GUIFigure, handles);


% --- Executes during object creation, after setting all properties.
function DataAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

xlabel(hObject, 'Time (ms)', 'FontSize', 14, 'FontWeight', 'bold'); 
ylabel(hObject, 'Amplitude (\muV)', 'FontSize', 14, 'FontWeight', 'bold');

% --- Plot data on main axes, according to current settings.
function handles = UpdatePlot(handles)

% collect some info:
sr = handles.Data.Info.SamplingRate_Hz;
dw = handles.Data.Info.DataWindow_ms;

% determine time vector:
t = ((dw(1)*sr/1000):(dw(2)*sr/1000))*1000/sr;

% determine which parameter to plot:
switch (handles.PlotByButtonGroup.SelectedObject.String)
    case 'Amplitudes' % single condition, channel:
        cond = str2double(handles.PlotEditText1.String);
        chan = str2double(handles.PlotEditText2.String);
        x = squeeze(handles.Data.ERPs(cond, chan, :, :));
        legstr = arrayfun(@(x) [num2str(x), ' mA'], handles.Data.Info.Amplitudes_mA, 'UniformOutput', false);
        
        % get peak data for markers:
        basedata = handles.Data.Peaks(cond).ECoG(chan).Baseline;
        latdata = cell(1, size(handles.Data.ERPs,4));
        for i = 1:size(latdata,2)
            latdata{i} = handles.Data.Peaks(cond).ECoG(chan).Latency(~isnan(handles.Data.Peaks(cond).ECoG(chan).Latency(:,i)),i);
        end
        
    case 'Channels' % single condition, amplitude:
        cond = str2double(handles.PlotEditText1.String);
        amp = str2double(handles.PlotEditText2.String);
        x = squeeze(handles.Data.ERPs(cond, :, :, amp))';
        if (iscell(handles.Data.Info.ECoGChans))
            legstr = handles.Data.Info.ECoGChans;
        else
            legstr = arrayfun(@(x,y) sprintf('%d-%d', x, y), handles.Data.Info.ECoGChans, handles.Data.Info.ECoGRefs, 'UniformOutput', false);
        end
        
        % gather data:
        basedata = zeros(2, size(handles.Data.ERPs,2));
        latdata = cell(1, size(handles.Data.ERPs,2));
        for i = 1:size(basedata,2)
            basedata(:,i) = handles.Data.Peaks(cond).ECoG(i).Baseline(:,amp);
            latdata{i} = handles.Data.Peaks(cond).ECoG(i).Latency(~isnan(handles.Data.Peaks(cond).ECoG(i).Latency(:,amp)), amp);
        end
        
    case 'Conditions' % single amplitude, channel:
        amp = str2double(handles.PlotEditText1.String);
        chan = str2double(handles.PlotEditText2.String);
        x = squeeze(handles.Data.ERPs(:, chan, :, amp))';
        legstr = arrayfun(@(x) num2str(x), 1:size(handles.Data.ERPs, 1), 'UniformOutput', false);
        
        % gather data:
        basedata = zeros(2, size(handles.Data.ERPs,1));
        latdata = cell(1, size(handles.Data.ERPs,1));
        for i = 1:size(basedata,2)
            basedata(:,i) = handles.Data.Peaks(i).ECoG(chan).Baseline(:,amp);
            latdata{i} = handles.Data.Peaks(i).ECoG(chan).Latency(~isnan(handles.Data.Peaks(i).ECoG(chan).Latency(:,amp)), amp);
        end
end

% make sure rows = samples:
if (size(x,2) > size(x,1))
    x = x';
end

% plot data (delete old data/markers):
hold(handles.DataAxes, 'off');
l = plot(handles.DataAxes, t, x);

% change lines > 8 to dashed to differentiate (colors recycle):
hold(handles.DataAxes, 'on');
for i = 8:size(x,2)
    l(i).LineStyle = '--';
end

% adjust plot display check box visibility and names:
for i = 1:11
    if (i <= size(latdata,2))
        handles.(['checkbox', num2str(i)]).Visible = 'on';
        handles.(['checkbox', num2str(i)]).String = legstr(i);
    else
        handles.(['checkbox', num2str(i)]).Visible = 'off';
    end
end

% set labels/legend:
xlabel(handles.DataAxes, 'Time (ms)', 'FontSize', 14, 'FontWeight', 'bold'); 
ylabel(handles.DataAxes, 'Amplitude (\muV)', 'FontSize', 14, 'FontWeight', 'bold');
handles.DataAxes.LineWidth = 1.5;
legend(legstr, 'AutoUpdate', 'off', 'LineWidth', 1);
handles.DataAxes.XLimMode = 'manual';
handles.DataAxes.YLimMode = 'manual';

% add lines to handles structure for reference:
handles.PlotLines = l;

% add a (blank) current point and selected point(s) markers:
handles.CurrentPointMarker = plot(handles.DataAxes, 0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2, 'Visible', 'off');
handles.SelectedPointMarkers = plot(handles.DataAxes, zeros(length(l),1), zeros(length(l),1), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'Visible', 'off');

% calculate peaks for each line, attach button down callback, display previously marked points:
for i = 1:length(l)
    y = l(i).YData; x = l(i).XData;
    pks = x([false, (y(2:end-1)-y(1:end-2)).*(y(3:end)-y(2:end-1)) < 0]);
    l(i).UserData.AllPeaks = pks; 
    
    l(i).ButtonDownFcn = @LineButtonDown_Callback;
    
    % set line visibility based on check boxes:
    if (handles.(['checkbox', num2str(i)]).Value < 0.5)
        vis = 'off';
    else
        vis = 'on';
    end
    l(i).Visible = vis;
    
    if (any(basedata(:,i))) % only plot marker if one has been marked
        l(i).UserData.BaselineMarker = plot(handles.DataAxes, basedata(1,i), basedata(2,i), 'go', 'Visible', vis);
    else
        l(i).UserData.BaselineMarker = [];
    end
    l(i).UserData.PeakMarkers = matlab.graphics.chart.primitive.Line.empty;
    for j = 1:size(latdata{i},1)
        l(i).UserData.PeakMarkers(j) = plot(handles.DataAxes, latdata{i}(j), y(x==latdata{i}(j)), 'ro', 'Visible', vis);
    end
end

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Update the peaks for a particular line.
function handles = UpdatePeaks(line, point, handles, op)

% determine index of the current line in the plot:
idx = find(handles.PlotLines == line);

% determine which condition/channel/amplitude we're updating:
switch (handles.PlotByButtonGroup.SelectedObject.String)
    case 'Amplitudes'
        cond = str2double(handles.PlotEditText1.String);
        chan = str2double(handles.PlotEditText2.String);
        amp = idx;
    case 'Channels'
        cond = str2double(handles.PlotEditText1.String);
        chan = idx;
        amp = str2double(handles.PlotEditText2.String);
    case 'Conditions'
        cond = idx;
        chan = str2double(handles.PlotEditText2.String);
        amp = str2double(handles.PlotEditText1.String);
end

% update the peak data:
if (strcmp(op, 'add')) % add new row - new peaks are their own group initially
    % update the plot:
    line.UserData.PeakMarkers(end+1) = plot(handles.DataAxes, point(1), point(2), 'ro');
    uistack(line.UserData.PeakMarkers(end), 'bottom');
    
    % re-calculate point amplitude based on current baseline:
    x = line.XData; y = line.YData;
    [~, idx] = min(abs(x - point(1)));
    point(2) = y(idx) - handles.Data.Peaks(cond).ECoG(chan).Baseline(2,amp);
    
    % determine how many groups there are currently:
    numgroups = 0;
    for i = 1:length(handles.Data.Peaks)
        for j = 1:length(handles.Data.Peaks(i).ECoG)
            numgroups = max([numgroups; handles.Data.Peaks(i).ECoG(j).Group], [], 1, 'omitnan');
        end
    end
    
    % update the data:
    handles.Data.Peaks(cond).ECoG(chan).Latency(end+1,:) = [nan(1,amp-1), point(1), nan(1,size(handles.Data.ERPs,4)-amp)];  
    handles.Data.Peaks(cond).ECoG(chan).Amplitude(end+1,:) = [nan(1,amp-1), point(2), nan(1,size(handles.Data.ERPs,4)-amp)];
    handles.Data.Peaks(cond).ECoG(chan).Group(end+1,1) = numgroups + 1;
    
elseif (strcmp(op, 'delete'))
    % update the data:
    idx = find(handles.Data.Peaks(cond).ECoG(chan).Latency(:,amp) == point.XData,1);
    handles.Data.Peaks(cond).ECoG(chan).Latency(idx,amp) = NaN;
    handles.Data.Peaks(cond).ECoG(chan).Amplitude(idx,amp) = NaN;
%     handles.Data.Peaks(cond).ECoG(chan).Group(idx) = NaN;
    
    % update the plot:
    idx = find(line.UserData.PeakMarkers == point);
    delete(line.UserData.PeakMarkers(idx));
    line.UserData.PeakMarkers(idx) = [];
end

% get rid of rows with all nans:
idxs = all(isnan(handles.Data.Peaks(cond).ECoG(chan).Latency) | isnan(handles.Data.Peaks(cond).ECoG(chan).Amplitude), 2);
handles.Data.Peaks(cond).ECoG(chan).Latency(idxs, :) = [];
handles.Data.Peaks(cond).ECoG(chan).Amplitude(idxs, :) = [];
handles.Data.Peaks(cond).ECoG(chan).Group(idxs, :) = [];

% update the peak tables:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);



% --- Update the baseline for a particular line.
function handles = UpdateBaseline(line, point, handles, op)

% determine index of the current line in the plot:
idx = find(handles.PlotLines == line);

% determine which condition/channel/amplitude we're updating:
switch (handles.PlotByButtonGroup.SelectedObject.String)
    case 'Amplitudes'
        cond = str2double(handles.PlotEditText1.String);
        chan = str2double(handles.PlotEditText2.String);
        amp = idx;
    case 'Channels'
        cond = str2double(handles.PlotEditText1.String);
        chan = idx;
        amp = str2double(handles.PlotEditText2.String);
    case 'Conditions'
        cond = idx;
        chan = str2double(handles.PlotEditText2.String);
        amp = str2double(handles.PlotEditText1.String);
end

% update the peak data:
if (strcmp(op, 'add'))
    % update plot:
    line.UserData.BaselineMarker = plot(handles.DataAxes, point(1), point(2), 'go');
    uistack(line.UserData.BaselineMarker, 'bottom');
    
    % update data:
    handles.Data.Peaks(cond).ECoG(chan).Baseline(:,amp) = point';
elseif (strcmp(op, 'delete'))
    % update plot:
    delete(line.UserData.BaselineMarker);
    line.UserData.BaselineMarker = [];
    
    % update data:
    handles.Data.Peaks(cond).ECoG(chan).Baseline(:,amp) = [0;0];
end

% re-calculate peak amplitudes based on new baseline:
x = line.XData; y = line.YData;
for i = 1:size(handles.Data.Peaks(cond).ECoG(chan).Amplitude,1)
    pklat = handles.Data.Peaks(cond).ECoG(chan).Latency(i,amp);
    if (~isnan(pklat))
        [~, idx] = min(abs(x - pklat));
        handles.Data.Peaks(cond).ECoG(chan).Amplitude(i,amp) = y(idx) - handles.Data.Peaks(cond).ECoG(chan).Baseline(2,amp);
    end
end

% update the baseline tables:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);


% -- update the data in the tables,
function handles = UpdateTables(handles)

% update peak tables, depending on which plot we're on:
switch (handles.PlotByButtonGroup.SelectedObject.String)
    case 'Amplitudes'
        cond = str2double(handles.PlotEditText1.String);
        chan = str2double(handles.PlotEditText2.String);
        
        % gather data:
        basedata = handles.Data.Peaks(cond).ECoG(chan).Baseline;
        ampdata = handles.Data.Peaks(cond).ECoG(chan).Amplitude;
        latdata = handles.Data.Peaks(cond).ECoG(chan).Latency;
        grpdata = repmat(handles.Data.Peaks(cond).ECoG(chan).Group, 1, size(ampdata,2));
    case 'Channels'
        cond = str2double(handles.PlotEditText1.String);
        amp = str2double(handles.PlotEditText2.String);
        
        % gather data:
        basedata = zeros(2, size(handles.Data.ERPs,2));
        ampdata = cell(1, size(handles.Data.ERPs,2));
        latdata = cell(1, size(handles.Data.ERPs,2));
        grpdata = cell(1, size(handles.Data.ERPs,2));
        for i = 1:size(basedata,2)
            basedata(:,i) = handles.Data.Peaks(cond).ECoG(i).Baseline(:,amp);
            ampdata{i} = handles.Data.Peaks(cond).ECoG(i).Amplitude(:, amp);
            latdata{i} = handles.Data.Peaks(cond).ECoG(i).Latency(:, amp);
            grpdata{i} = handles.Data.Peaks(cond).ECoG(i).Group;
        end
        ampdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,ampdata))-length(x),1)], ampdata, 'UniformOutput', false)); % pad with NaNs
        latdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,latdata))-length(x),1)], latdata, 'UniformOutput', false)); % pad with NaNs
        grpdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,grpdata))-length(x),1)], grpdata, 'UniformOutput', false)); % pad with NaNs
    case 'Conditions'
        chan = str2double(handles.PlotEditText2.String);
        amp = str2double(handles.PlotEditText1.String);
        
        % gather data:
        basedata = zeros(2, size(handles.Data.ERPs,1));
        ampdata = cell(1, size(handles.Data.ERPs,1));
        latdata = cell(1, size(handles.Data.ERPs,1));
        grpdata = cell(1, size(handles.Data.ERPs,1));
        for i = 1:size(basedata,2)
            basedata(:,i) = handles.Data.Peaks(i).ECoG(chan).Baseline(:,amp);
            ampdata{i} = handles.Data.Peaks(i).ECoG(chan).Amplitude(:, amp);
            latdata{i} = handles.Data.Peaks(i).ECoG(chan).Latency(:, amp);
            grpdata{i} = handles.Data.Peaks(i).ECoG(chan).Group;
        end
        ampdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,ampdata))-length(x),1)], ampdata, 'UniformOutput', false)); % pad with NaNs
        latdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,latdata))-length(x),1)], latdata, 'UniformOutput', false)); % pad with NaNs 
        grpdata = cell2mat(cellfun(@(x) [x;nan(max(cellfun(@length,grpdata))-length(x),1)], grpdata, 'UniformOutput', false)); % pad with NaNs
end

% remove data rows with all NaNs:
% idx = all(isnan(ampdata) | isnan(latdata), 2);
% ampdata(idx, :) = []; latdata(idx, :) = []; grpdata(idx, :) = [];

% re-organize rows by group:
% allgroups = unique(grpdata(~isnan(grpdata)));
allgroups = unique(grpdata(~isnan(grpdata) & ~isnan(ampdata) & ~isnan(latdata)));
newampdata = NaN(length(allgroups), size(ampdata,2));
newlatdata = NaN(length(allgroups), size(latdata,2));
for i = 1:length(allgroups)
%     newampdata(i,:) = max(ampdata(grpdata == allgroups(i),:),[],1,'omitnan');
%     newlatdata(i,:) = max(latdata(grpdata == allgroups(i),:),[],1,'omitnan');
    idx = find(grpdata == allgroups(i)); [~, c] = ind2sub(size(ampdata), idx);
    newampdata(i,c) = ampdata(idx);
    newlatdata(i,c) = latdata(idx);
end

% update baseline tables:
handles.BaseLatTable.Data = sprintfc('%.2f', basedata(1,:));
handles.BaseLatTable.ColumnName = sprintfc('%d', 1:size(ampdata,2));
handles.BaseLatTable.ColumnWidth = num2cell(40*ones(1,size(ampdata,2)));
handles.BaseAmpTable.Data = sprintfc('%.2f', basedata(2,:));
handles.BaseAmpTable.ColumnName = sprintfc('%d', 1:size(ampdata,2));
handles.BaseAmpTable.ColumnWidth = num2cell(40*ones(1,size(ampdata,2)));

% update peak tables:
if (isrow(allgroups))
    group = sprintfc('%d', allgroups');
else
    group = sprintfc('%d', allgroups);
end
handles.PeakLatTable.Data = [group, strrep(sprintfc('%.2f', newlatdata), 'NaN', ' ')];
handles.PeakLatTable.ColumnName = [{'Grp'}, sprintfc('%d', 1:size(ampdata,2))];
handles.PeakLatTable.ColumnWidth = [{25}, num2cell(40*ones(1,size(ampdata,2)))];
handles.PeakAmpTable.Data = [group, strrep(sprintfc('%.2f', newampdata), 'NaN', ' ')];
handles.PeakAmpTable.ColumnName = [{'Grp'}, sprintfc('%d', 1:size(ampdata,2))];
handles.PeakAmpTable.ColumnWidth = [{25}, num2cell(40*ones(1,size(ampdata,2)))];

% update handles structure:
guidata(handles.GUIFigure, handles);



function PlotEditText_Callback(hObject, eventdata, handles)
% hObject    handle to PlotEditText1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update plot:
handles = UpdatePlot(handles);

% update table display:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes during object creation, after setting all properties.
function PlotEditText1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotEditText1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function PlotEditText2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotEditText2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in PlotByButtonGroup.
function PlotByButtonGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in PlotByButtonGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% determine which condition/channel/amplitude we're on now:
switch (hObject.String)
    case 'Amplitudes'
        handles.PlotText1.String = 'Condition:';
        handles.PlotText2.String = 'Channel:';
    case 'Channels'
        handles.PlotText1.String = 'Condition:';
        handles.PlotText2.String = 'Amplitude:';
    case 'Conditions'
        handles.PlotText1.String = 'Amplitude:';
        handles.PlotText2.String = 'Channel:';
end

% reset plot display check boxes to checked (i.e. display everything):
for i = 1:10
    handles.(['checkbox',num2str(i)]).Value = 1;
end

% update plot:
handles = UpdatePlot(handles);

% update table display:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes when entered data in editable cell(s) in BrainAreaTable.
function BrainAreaTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to BrainAreaTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% save string to data structure:
handles.Data.Info.ECoGBrainAreas{eventdata.Indices(2)} = eventdata.EditData;

% get rid of selected cell background color (ugh):
hObject.Data = cell(size(handles.Data.Info.ECoGBrainAreas));
hObject.Data = handles.Data.Info.ECoGBrainAreas;

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes when selected object is changed in CurrentActionButtonGroup.
function CurrentActionButtonGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in CurrentActionButtonGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set defaults:
zoom(handles.DataAxes, 'off');
handles.DataAxes.ButtonDownFcn = [];

% determine current action:
switch (hObject.String)
    case 'Zoom'
        zoom(handles.DataAxes, 'on');
        
    case {'Group points - merge', 'Group points - new'}
        handles.DataAxes.ButtonDownFcn = @DataAxes_ButtonDownFcn;
        
end

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes on mouse press over axes background.
function DataAxes_ButtonDownFcn(hObject, eventdata)
% hObject    handle to DataAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current handles object:
handles = guidata(ancestor(hObject, 'figure'));

% set figure callbacks to detect mouse motion and button release:
handles.GUIFigure.WindowButtonMotionFcn = @GUIFigure_WindowButtonMotionFcn;
handles.GUIFigure.WindowButtonUpFcn = @GUIFigure_WindowButtonUpFcn;

% initialize a selection rectangle:
x = hObject.CurrentPoint(1,1); y = hObject.CurrentPoint(1,2);
handles.SelectionRect = line(hObject, [x, x, x, x, x], [y, y, y, y, y]);

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes on mouse motion over figure - except title and menu.
function GUIFigure_WindowButtonMotionFcn(hObject, eventdata)
% hObject    handle to GUIFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current handles object:
handles = guidata(ancestor(hObject, 'figure'));

% update selection rectangle size:
x = handles.DataAxes.CurrentPoint(1,1); y = handles.DataAxes.CurrentPoint(1,2);
handles.SelectionRect.XData(3:4) = x;
handles.SelectionRect.YData(2:3) = y;

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function GUIFigure_WindowButtonUpFcn(hObject, eventdata)
% hObject    handle to GUIFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current handles object:
handles = guidata(ancestor(hObject, 'figure'));

% clear figure callbacks:
handles.GUIFigure.WindowButtonMotionFcn = [];
handles.GUIFigure.WindowButtonUpFcn = [];

% get currently displayed data:
lats = str2double(strrep(handles.PeakLatTable.Data(:,2:end), ' ', 'NaN'));
amps = str2double(strrep(handles.PeakAmpTable.Data(:,2:end), ' ', 'NaN'));
baseamp = repmat(str2double(handles.BaseAmpTable.Data), size(lats,1), 1);
grps = str2double(handles.PeakLatTable.Data(:,1));

% determine which markers are inside the selection rectangle:
idxs = inpolygon(lats, amps+baseamp, handles.SelectionRect.XData, handles.SelectionRect.YData);
[rows, cols] = find(idxs);

% make sure all selected points are visible:
notvislines = find(arrayfun(@(x) strcmp(x.Visible, 'off'), handles.PlotLines));
idxs = find(ismember(cols, notvislines));
rows(idxs) = []; cols(idxs) = [];

% clear selection rectangle:
delete(handles.SelectionRect);
handles.SelectionRect = [];

% stop if no points selected:
if (isempty(cols))
    return;
end

% don't allow two peaks on the same amplitude/channel/condition group:
if (length(unique(cols)) ~= length(cols))
    msgbox('Two peaks on the same ERP can''t be in the same group.', 'Error', 'error', 'modal');
    return;
end

% which grouping function are we on (merge vs new group):
if (strcmp(handles.CurrentActionButtonGroup.SelectedObject.String, 'Group points - merge')) % merge rows/groups together
    
    % make sure there isn't overlap between the groups to merge:
    for i = 1:length(handles.Data.Peaks)
        for j = 1:length(handles.Data.Peaks(i).ECoG)
            grpmem = ismember(handles.Data.Peaks(i).ECoG(j).Group, grps(rows));
            if (sum(grpmem) > 1 && any(~isnan(sum(handles.Data.Peaks(i).ECoG(j).Latency(grpmem,:),1))))
                badidx = find(~isnan(sum(handles.Data.Peaks(i).ECoG(j).Latency(grpmem,:),1)), 1);
                msgbox(sprintf('Group conflict at Condition %d, Channel %d, Amplitude %d - merge not possible', i, j, badidx), 'Error', 'error', 'modal');
                return;
            end
        end
    end
    
    % merge groups for each condition/channel:
    grouplabels = grps(rows);
    for i = 1:length(handles.Data.Peaks)
        for j = 1:length(handles.Data.Peaks(i).ECoG)
            
            grprows = find(ismember(handles.Data.Peaks(i).ECoG(j).Group, grouplabels));
            if (length(grprows) > 1)    % merge rows
                handles.Data.Peaks(i).ECoG(j).Amplitude(grprows(1),:) = max(handles.Data.Peaks(i).ECoG(j).Amplitude(grprows,:), [], 1, 'omitnan');
                handles.Data.Peaks(i).ECoG(j).Amplitude(grprows(2:end),:) = [];
                
                handles.Data.Peaks(i).ECoG(j).Latency(grprows(1),:) = max(handles.Data.Peaks(i).ECoG(j).Latency(grprows,:), [], 1, 'omitnan');
                handles.Data.Peaks(i).ECoG(j).Latency(grprows(2:end),:) = [];
                
                handles.Data.Peaks(i).ECoG(j).Group(grprows(1)) = grouplabels(1);
                handles.Data.Peaks(i).ECoG(j).Group(grprows(2:end)) = [];
            elseif (length(grprows) == 1)   % just change group number
                handles.Data.Peaks(i).ECoG(j).Group(grprows(1)) = grouplabels(1);
            end
            
        end
    end
    
else    % make selected points into a new group
    
    % determine how many groups there are currently:
    numgroups = 0;
    for i = 1:length(handles.Data.Peaks)
        for j = 1:length(handles.Data.Peaks(i).ECoG)
            numgroups = max([numgroups; handles.Data.Peaks(i).ECoG(j).Group], [], 1, 'omitnan');
        end
    end
    
    % which peaks have we selected:
%     groups = grps(rows);
    switch (handles.PlotByButtonGroup.SelectedObject.String)
        case 'Amplitudes'
            conds = str2double(handles.PlotEditText1.String);
            chans = str2double(handles.PlotEditText2.String);
            amps = cols;
        case 'Channels'
            conds = str2double(handles.PlotEditText1.String);
            chans = cols;
            amps = str2double(handles.PlotEditText2.String);
        case 'Conditions'
            conds = cols;
            chans = str2double(handles.PlotEditText1.String);
            amps = str2double(handles.PlotEditText2.String);
    end
    groups = zeros(size(handles.Data.ERPs,1),size(handles.Data.ERPs,2),size(handles.Data.ERPs,4));
    groups(conds, chans, amps) = grps(rows);
    
    for i = 1:length(conds)
        for j = 1:length(chans)
            % which points are contained in this condition/channel:
%             groups = grps(rows(ismember(grps(rows), handles.Data.Peaks(conds(i)).ECoG(chans(j)).Group)));
            
            % which rows/groups are we changing in this condition/channel:
            rowidxs = arrayfun(@(x) find(handles.Data.Peaks(conds(i)).ECoG(chans(j)).Group == x,1), squeeze(groups(conds(i),chans(2),amps)));
            
            % get the indices of the selected points before and after adding the new group/row:
            pointidxs = sub2ind(size(handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude), rowidxs, amps);
            newpointidxs = pointidxs + (amps-1);
            
            % add points to new row/group:
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude(end+1,amps) = handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude(pointidxs);
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude(end,setdiff(1:size(handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude,2), amps)) = nan;
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency(end+1,amps) = handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency(pointidxs);
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency(end,setdiff(1:size(handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency,2), amps)) = nan;
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Group(end+1,1) = numgroups + 1;
            
            % delete points from old rows/groups:
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude(newpointidxs) = nan;
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency(newpointidxs) = nan;
            
            % delete empty groups:
            allnan = all(isnan(handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude), 2);
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Amplitude(allnan,:) = [];
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Latency(allnan,:) = [];
            handles.Data.Peaks(conds(i)).ECoG(chans(j)).Group(allnan) = [];
        end
    end
    
end

% % determine which condition/channel/amplitude those peaks belong to:
% switch (handles.PlotByButtonGroup.SelectedObject.String)
%     case 'Amplitudes'
%         
%         % which grouping function are we on (merge vs new group):
%         if (strcmp(handles.CurrentActionButtonGroup.SelectedObject.String, 'Group points - merge')) % merge rows/groups together
%             
%             % make sure there isn't overlap between the groups to merge:
%             grouplabels = grps(rows);
%             for i = 1:length(handles.Data.Peaks)
%                 for j = 1:length(handles.Data.Peaks(i).ECoG)
%                     grpmem = ismember(handles.Data.Peaks(i).ECoG(j).Group, grouplabels);
%                     if (length(grpmem) > 1 && any(~isnan(sum(handles.Data.Peaks(i).ECoG(j).Latency(grpmem,:),1))))
%                         msgbox(sprintf('Group conflict at Condition %d, Channel %d - merge not possible', i, j), 'Error', 'error', 'modal');
%                         return;
%                     end
%                 end
%             end
% 
%             % merge groups for each condition/channel:
%             grouplabels = grps(rows);
%             for i = 1:length(handles.Data.Peaks)
%                 for j = 1:length(handles.Data.Peaks(i).ECoG)
%                     
%                     grprows = find(ismember(handles.Data.Peaks(i).ECoG(j).Group, grouplabels));
%                     if (length(grprows) > 1)    % merge rows
%                         handles.Data.Peaks(i).ECoG(j).Amplitude(grprows(1),:) = max(handles.Data.Peaks(i).ECoG(j).Amplitude(grprows,:), [], 1, 'omitnan');
%                         handles.Data.Peaks(i).ECoG(j).Amplitude(grprows(2:end),:) = [];
%                         
%                         handles.Data.Peaks(i).ECoG(j).Latency(grprows(1),:) = max(handles.Data.Peaks(i).ECoG(j).Latency(grprows,:), [], 1, 'omitnan');
%                         handles.Data.Peaks(i).ECoG(j).Latency(grprows(2:end),:) = [];
%                         
%                         handles.Data.Peaks(i).ECoG(j).Group(grprows(1)) = grouplabels(1);
%                         handles.Data.Peaks(i).ECoG(j).Group(grprows(2:end)) = [];
%                     elseif (length(grprows) == 1)   % just change group number
%                         handles.Data.Peaks(i).ECoG(j).Group(grprows(1)) = grouplabels(1);
%                     end
%                     
%                 end
%             end
%         
%         else    % make selected points into a new group
%             
%             % determine how many groups there are currently:
%             numgroups = 0;
%             for i = 1:length(handles.Data.Peaks)
%                 for j = 1:length(handles.Data.Peaks(i).ECoG)
%                     numgroups = max([numgroups; handles.Data.Peaks(i).ECoG(j).Group], [], 1, 'omitnan');
%                 end
%             end
%             
%             % which condition/channel are we on:
%             cond = str2double(handles.PlotEditText1.String);
%             chan = str2double(handles.PlotEditText2.String);
%             
%             % add points to new row/group:
%             handles.Data.Peaks(cond).ECoG(chan).Amplitude(end+1,cols) = handles.Data.Peaks(cond).ECoG(chan).Amplitude(idxs);
%             handles.Data.Peaks(cond).ECoG(chan).Amplitude(end,setdiff(1:size(handles.Data.Peaks(cond).ECoG(chan).Amplitude,2), cols)) = nan;
%             handles.Data.Peaks(cond).ECoG(chan).Latency(end+1,cols) = handles.Data.Peaks(cond).ECoG(chan).Latency(idxs);
%             handles.Data.Peaks(cond).ECoG(chan).Latency(end,setdiff(1:size(handles.Data.Peaks(cond).ECoG(chan).Latency,2), cols)) = nan;
%             handles.Data.Peaks(cond).ECoG(chan).Group(end+1,1) = numgroups + 1;
%             
%             % delete points from old rows/groups:
%             handles.Data.Peaks(cond).ECoG(chan).Amplitude(rows,cols) = nan;
%             handles.Data.Peaks(cond).ECoG(chan).Latency(rows,cols) = nan;
%             
%             % delete empty groups:
%             allnan = all(isnan(handles.Data.Peaks(cond).ECoG(chan).Amplitude), 2);
%             handles.Data.Peaks(cond).ECoG(chan).Amplitude(allnan,:) = [];
%             handles.Data.Peaks(cond).ECoG(chan).Latency(allnan,:) = [];
%             handles.Data.Peaks(cond).ECoG(chan).Group(allnan) = [];
%             
%         end
%         
%     case 'Channels'
%         cond = repmat(str2double(handles.PlotEditText1.String), length(cols), 1);
%         chan = cols;
% %         amp = repmat(str2double(handles.PlotEditText2.String), length(cols), 1);
%         
%         %%% these should be easier, i just need to change the group number,
%         %%% right????????
%         
%         for i = 2:length(rows)
%             handles.Data.Peaks(cond(i)).ECoG(chan(i)).Group(handles.Data.Peaks(cond(i)).ECoG(chan(i)).Group == grps(rows(i))) = grps(rows(1));
%         end
%         
%     case 'Conditions'
%         cond = cols;
%         chan = repmat(str2double(handles.PlotEditText2.String), length(cols), 1);
%         amp = repmat(str2double(handles.PlotEditText1.String), length(cols), 1);
% end

% figure out which groups are present and accounted for:
groups = cell(length(handles.Data.Peaks), length(handles.Data.Peaks(1).ECoG));
for i = 1:length(handles.Data.Peaks)
    for j = 1:length(handles.Data.Peaks(i).ECoG)
        groups{i,j} = handles.Data.Peaks(i).ECoG(j).Group(~isnan(handles.Data.Peaks(i).ECoG(j).Group));
    end
end
groupnums = unique(cell2mat(reshape(groups, numel(groups), 1)));

% re-name groups to keep a sequential order:
newgroupnums = (1:length(groupnums))';
for i = 1:length(handles.Data.Peaks)
    for j = 1:length(handles.Data.Peaks(i).ECoG)
        handles.Data.Peaks(i).ECoG(j).Group = changem(handles.Data.Peaks(i).ECoG(j).Group, newgroupnums, groupnums);
    end
end

% update table display:
handles = UpdateTables(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes on button press in any checkbox.
function checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to calling checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear any table selection so its not confusing:
olddata = handles.PeakAmpTable.Data;
handles.PeakAmpTable.Data = cell(size(olddata));
handles.PeakAmpTable.Data = olddata;
olddata = handles.PeakLatTable.Data;
handles.PeakLatTable.Data = cell(size(olddata));
handles.PeakLatTable.Data = olddata;

% update plot:
handles = UpdatePlot(handles);

% update handles structure:
guidata(handles.GUIFigure, handles);


% --- Executes when selected cell(s) is changed in PeakLatTable.
function Table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to PeakLatTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

% protect against inadvertent callbacks:
if (isempty(eventdata.Indices))
    return;
end

% clear the other table selection so its not confusing:
if (hObject == handles.PeakLatTable)
    olddata = handles.PeakAmpTable.Data;
    handles.PeakAmpTable.Data = cell(size(olddata));
    handles.PeakAmpTable.Data = olddata;
else
    olddata = handles.PeakLatTable.Data;
    handles.PeakLatTable.Data = cell(size(olddata));
    handles.PeakLatTable.Data = olddata;
end

% get all the latencies for this group (whole row):
lats = str2double(handles.PeakLatTable.Data(eventdata.Indices(1), 2:end));
lineidxs = find(~isnan(lats));
lats = lats(~isnan(lats));

% gather the x/y data for the selected markers:
xdata = zeros(1,length(lats));
ydata = zeros(1,length(lats));
for i = 1:length(xdata)
    
    if (strcmp(handles.PlotLines(lineidxs(i)).Visible, 'off'))
        xdata(i) = NaN;
        ydata(i) = NaN;
    else
        [~, markidx] = min(abs(lats(i) - [handles.PlotLines(lineidxs(i)).UserData.PeakMarkers.XData]));
        xdata(i) = handles.PlotLines(lineidxs(i)).UserData.PeakMarkers(markidx).XData(1);
        ydata(i) = handles.PlotLines(lineidxs(i)).UserData.PeakMarkers(markidx).YData(1);
    end
    
end

% set the selection x/y data and make it visible:
handles.SelectedPointMarkers.XData = xdata;
handles.SelectedPointMarkers.YData = ydata;
handles.SelectedPointMarkers.Visible = 'on';
