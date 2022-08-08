function [freq,y] = my_dft(x,Fs)

dt=1/Fs;
l_x=length(x);
tnew=0:dt:l_x*dt-dt;
%%% FFT Parameters
ntemps = length(tnew);
puis2 = nextpow2(ntemps);
newn = 2^puis2;
nfourier = newn;
npt = nfourier/2;
freq = 1/dt*(0:npt-1)/nfourier; % Frequency vector
%%% Compute FFT of xnew
tff = fft(x,nfourier);
tff(npt+1:end) = [];
y = (tff)/(length(x)/nfourier/2);

end