clear
clc
opengl hardware
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)
%% Constructing signal
t = linspace(0,70,2^10)';
dt = t(2)-t(1);
Fs = 2*pi/dt;
xx = zeros(length(t),3);
ww1 = 1 + 4*exp(-.025*t);
td = 50;
ww2 = (5+2*(cos(pi/td*t).*(t<=td)-(t>td)));
xx(:,1) = 1/4*exp(-.05*t).*sin(cumtrapz(t,ww1));
xx(:,2) = exp(-.05*t).*cos(cumtrapz(t,ww2));
xnew = sum(xx,2);
figure
subplot(2,1,1)
plot(t,xnew, 'k')
xlabel('Time')
ylabel('x(t)')
subplot(2,1,2)
[~,freq_recon,MODS_orig] = freq_inst_morlet(xnew,Fs,0,15,200,2);
contourf(t,freq_recon,MODS_orig', 30, 'linestyle', 'none')
colormap(1-gray.^.5)
xlabel('Time, t')
ylabel('Frequency, $\omega$')
%% Padding signal based on its frequency, phase and amplitude
sample_time_interval = 2;
t_sample = t(1:sum(t<=sample_time_interval));
x_sample = xnew(1:sum(t<=sample_time_interval));

[xData, yData] = prepareCurveData( t_sample, x_sample );

% Set up fittype and options.
ft = fittype( 'sin2' ); % sink where k is the number of harmonics in xnew
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf];
opts.StartPoint = [0.7650630511096 6.33270794999479 2.23058102767865 0.257345195058598 9.49906192499219 -0.739534210172849];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


t_new = [transpose(-512:-1)*dt;t];
x_new = [fitresult(transpose(-512:-1)*dt);xnew];
%% Padded signal
dt=(t(2)-t(1));
FS = 2*pi/dt;
Fo = 2;
figure
subplot(2,1,1)
plot(t_new,x_new,'k','linewidth',1.5)
xlabel('Time, t')
ylabel('y(t)')
% xlim([0 50])
subplot(2,1,2)
[~,frequency_WT,MODS] = freq_inst_morlet(x_new,FS,0,15,200,Fo);
contourf(t_new,frequency_WT,sqrt((MODS')),'linestyle','none')
xlabel('Time, t')
ylabel('Frequency, $\omega$')
title('CWT of $y$(t)')
colormap(1-gray.^(0.5))
