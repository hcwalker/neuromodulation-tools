%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HPF 2Hz to remove voltage drift from EEG/EMG/ECoG signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = filterVoltageDrift(x,Fs)
x = double(x);

%% HPF 2Hz
filterOrder = 3;
cutOffFreq = 2;
[Bcoef, Acoef] = butter(filterOrder, cutOffFreq/(Fs/2)); % Generate filter coefficients
y = x-filtfilt(Bcoef, Acoef, x); % Apply filter to data using zero-phase filtering
