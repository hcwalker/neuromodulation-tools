%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove 60Hz and its harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = filter60Hz(x,Fs)
q = 65; 
bw = (120/Fs)/q;
[Bcoef,Acoef] = iircomb(round(Fs/60),bw,'notch'); 
y = filtfilt(Bcoef,Acoef,x);
