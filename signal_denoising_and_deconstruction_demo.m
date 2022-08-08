clear
clc
opengl hardware
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)
%% Constructing signal
t = linspace(0,70,2^13)';
dt=(t(2)-t(1));
Fs = 2*pi/dt;
xx = zeros(length(t),3);
ww1 = 1 + 4*exp(-.025*t);
td = 50;
ww2 = (5+2*(cos(pi/td*t).*(t<=td)-(t>td)));
xx(:,1) = 1/4*exp(-.025*t).*sin(cumtrapz(t,ww1));
xx(:,2) = exp(-.025*t).*cos(cumtrapz(t,ww2));
x = sum(xx,2);
figure
subplot(2,1,1)
plot(t,x)
xlabel('Time')
ylabel('x(t)')
subplot(2,1,2)
[~,freq_recon,MODS_orig] = freq_inst_morlet(x,Fs,0,15,200,2);
contourf(t,freq_recon,MODS_orig', 30, 'linestyle', 'none')
colormap(1-gray.^.5)
xlabel('Time, t')
ylabel('Frequency, $\omega$')
%% ICWT parameters
dt=(t(2)-t(1));
Fs = 2*pi/dt;
representative_signal = 1;
Fo = 2;
fi = 0;
ff = 10;
nF = 200;
N_modes = 2;
WT_pow = .75;
WT_plot = "contourf";
N_contourf = 30;
%% ICWT functions
[t,x_new,interval_freq,modes] = Fast_Morlet_ICWT_Harmonic_decomp(x,representative_signal,Fs,dt,Fo,fi,ff,nF,N_modes,WT_pow,WT_plot,N_contourf);
%% Extracted modes comparison
NNN = N_modes;
OMEGA_SIG = [ww2 ww1];
figure
for counterf = 1:size(modes,2)
    subplot(NNN,2,(counterf-1)*2+1)
    plot(t,xx(:,counterf),'k',t,modes(:,counterf),'-.r','linewidth',1.2)
    xlabel('Time, t')
%     ylabel(sprintf('mode %i',counterf))
%     title('timeseries')
    legend('Exact','Extracted')
    subplot(NNN,2,(counterf-1)*2+2)
    [~,FF,MODS] = freq_inst_morlet(modes(:,counterf),Fs,0,15,200,Fo);
    contourf(t,FF,MODS')    
    Fnyq = Fs/2;
    Fc = 2/Fnyq;
    [b,a] = butter(3,Fc,'low');
    
    Sig = hilbert(modes(:,counterf));
    phase = unwrap(angle(Sig));
    Frequency = gradient(phase)/(t(2)-t(1));
    Frequency = filtfilt(b,a,double(Frequency));
    hold on
    plot(t,Frequency,'r',t,OMEGA_SIG(:,3-counterf),'-.g','linewidth',2)
    legend('','$\tilde{\Omega}$','${\Omega}$')
    xlabel('Time, t')
    ylabel('Frequency, $\omega$')
    title('Extracted')
    colormap(1-gray.^.5)
end
x_recon = sum(modes,2);
S1 = sum((x-mean(x)).^2);
S2 = sum((x-x_recon).^2);
R2 = 1-S2/S1;
figure
plot(t,x,'k',t,x_recon,'-.r','linewidth',2)
xlabel('t')
ylabel('signal')
legend('Original','Reconstructed')
title(sprintf('R-squared: %f',R2))