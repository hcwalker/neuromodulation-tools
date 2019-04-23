%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter high frequency noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = filterNoise(x,Fs)
%y = sgolayfilt(x,3,11);
y = wiener2(x,[1 5]);