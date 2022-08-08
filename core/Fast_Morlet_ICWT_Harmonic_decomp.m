function [t,x_new,interval_freq,modes] = Fast_Morlet_ICWT_Harmonic_decomp(x,representative_signal,Fs,dt,Fo,fi,ff,nF,N_modes,WT_pow,WT_plot,N_contourf)
% x is a collection of column vectors
% representative_signal: Corresponds to the signal with the richest frequency content (in case x is a collection of column vectors)
% Fs: Sampling frequency
% dt: time increment duration
% Fo: Mother (Morlet) wavelet center frequency
% fi: Lower boundary of the wavelet transform's frequency
% ff: Upper boundary of the wavelet transform's frequency
% nF: Number of frequencies in the [fi, ff] interval
% N_modes: Number modes to be extracted from the original signal
% WT_pow: Power that the wavelet transform is rased to, so that more (WT_pow < 1) or less (WT_pow > 1) details are observed.
% WT_plot = 'contourf' or 'imagesc'


if size(x,2)==1
    representative_signal = 1;
end

NDOF = size(x,2);
l_new = 2^nextpow2(size(x,1));
x_new = zeros(l_new,NDOF);
x_new(1:size(x,1),:) = x;
t = 0:dt:l_new*dt-dt;
nfourier = l_new;
npt = nfourier/2;
freq = Fs*([0:nfourier-1])/nfourier;
FREQ = Fs*([0:npt-1 npt:-1:1])/nfourier;
interval_freq = FREQ;
a = Fo./interval_freq;
fft_MW = conj(bsxfun(@times,pi^(1/4)*(2^0.5)*(exp(-0.5*(2*pi*(bsxfun(@times,FREQ',a)-Fo)).^2)-1*exp(-0.5*(2*pi^2*(bsxfun(@times,FREQ',a).^2+Fo.^2)))),sqrt(a)));

[~,frequency_WT,MODS] = freq_inst_morlet(x_new(:,representative_signal),Fs,fi,ff,nF,Fo);
figure
if WT_plot == "imagcesc"
imagesc((MODS').^WT_pow)
set(gca,'ydir','nor')
colormap(1-gray.^(1/2))
elseif WT_plot == "contourf"
contourf(MODS'.^WT_pow, N_contourf,'linestyle','none')
colormap(1-gray.^(1/2))
end

modes = zeros(length(t),N_modes,NDOF);
for counter = 1:N_modes
    h = impoly;
    poss = getPosition(h);
    %     Mask = poly2mask(poss(:,1),poss(:,2),size(MODS',1),size(MODS',2));
    freq_repmat = repmat(freq,length(poss(:,2)),1);
    patch_freq = frequency_WT(1,ceil(poss(:,2)));
    indices = sum(freq_repmat<=patch_freq',2);
    Mask2 = poly2mask(poss(:,1),indices,length(freq)/2,length(t));
    for counter_NDOF = 1:NDOF
        tff = fft(x_new(:,counter_NDOF),nfourier);
        noyau2 = bsxfun(@times,fft_MW,tff);
        resut2 = ifft(noyau2,nfourier);
        Patch = [Mask2' fliplr(Mask2')].*resut2;
        fft_Patch = fft(Patch,nfourier);
        fft_mode_2D = fft_Patch./fft_MW;
        fft_mode = diag(fft_mode_2D);
        fft_mode(isnan(fft_mode)) = 0;
        modes(:,counter,counter_NDOF) = ifft(fft_mode,nfourier,'symmetric');
    end
    delete(h)
end
end