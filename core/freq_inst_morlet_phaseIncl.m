function [tnew,interval_freq,module,resut3] = freq_inst_morlet_phaseIncl(xnew,FS,fi,ff,nf,Fo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet Transform using Morlet Wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs
% ------
% xnew = signal with dimensions N x 1
% Fs = sampling frequency
% fi = low frequency bound
% ff = high frequency bound
% nf = number of frequency points (higher is results in finer plots, 100 
% or 200 are good starting points)
% Fo = Mother wavelet frequency (2 or 4 is best)
%
% Outputs
% -------
% tnew = Time vector based on length of xnew and FS
% interval_freq = Frequency vector
% module = Absolute value of wavelet transform of xnew

%% Transform Parameters
dt=1/FS;
l_x=length(xnew);
tnew=0:dt:l_x*dt-dt;
df = (ff-fi)/nf;
interval_freq=fi:df:ff;
a = Fo./interval_freq;

%% FFT Parameters
ntemps = length(tnew);
puis2 = nextpow2(ntemps);
newn = 2^puis2;
nfourier = newn;
npt = nfourier/2;
freq = 1/dt*(0:npt-1)/nfourier; % Frequency vector

%% Compute FFT of xnew
tff = fft(xnew,nfourier);
tff(npt+1:end) = [];

%% Vectorized Computation of Wavelet Transform of xnew
noyau2 = bsxfun(@times,conj(bsxfun(@times,(2^0.5)*exp(-0.5*(2*pi*(bsxfun(@times,freq',a)-Fo)).^2),sqrt(a))),tff);
noyau2(freq == 0,:) = 0;
resut2 = ifft(noyau2,nfourier);
resut3 = resut2(1:ntemps,:);
module  = abs(resut3);
