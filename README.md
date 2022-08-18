# Continuous Wavelet Transform based Signal Processing tool for time-frequency denoising and deconstruction of time series
Time-frequency analysis of any given time series or signals in general, can provide valuable information about their nature. Information such as frequency (mode) content, stationarity and non-stationarity, noise pollution, noise color, etc. One of the most powerful techniques for time-frequency analysis and manipulation of time series is Empirical Mode Decomposition ([EMD](https://www.ripublication.com/irph/ijeee_spl/ijeeev7n8_14.pdf)) and its derivatives, however, they fail in many scenarios where signal-to-noise ratios are relatively small or where the signal has a high-velocity non-stationary nature or the signal has weak but important underlying frequency content (modes). Because of this, I developed a [Wavelet-Based time-frequency analysis techinque](https://www.sciencedirect.com/science/article/abs/pii/S0888327021000868) that overcomes all of the aforementioned issues.

In this project I take advantage of Parseval's theorem, fast Fourier transfor and its inverse for fast and accurate computation of the continuous wavelet transform of 1-dimensional signals and the inverse continuous wavelet transforms, see `freq_inst_morlet.m` and `Fast_Morlet_ICWT_Harmonic_decomp.m` source codes. I have used these codes in the following applications:

- Time-frequency decomposition of time series (an improved and more accurate replacement for EMD-based mothods);

- Signal denoising;

- Custom time-frequency filtration.

In the process of developing this technique (please refer to `signal_denoising_and_deconstruction_demo.m` to see a demo of it implementation) I noticed that [boundary effects](https://www.mathworks.com/help/wavelet/ug/boundary-effects-and-the-cone-of-influence.html) of wavelet transform can sometimes make it inconvenient to properly interpret time-frequency contents of the signals and more importantly, interfere with accurate computation of the inverse wavelet transforms. As I failed to find any articles, books or codes that addressd this problem outside the scope of traditional signal padding, I developed a time-frequency extrapolation algorithm that eliminates the aforementioned boundary effects (please see `Boundary_effect_treatment.m` for a simple implementaion demo).

## Getting Started
The codes were written, run and tested by MATLAB 2018a.

The codes are self-contained, so everything should run smoothly.

## Tutorial
For the sake of brevity I'll provide a tutorial on the `signal_denoising_and_deconstruction_demo.m` code. The boundary effect treatment code (`Boundary_effect_treatment.m`) is concisely written and should be relatively easy to undestand on its own.

In the demo code included in this repository I show the implementation of the time-frequency analysis technique I developed to separate all the nonstationary time-scales of a time series. The figure below is an artificially generated time series with two known nonstationary time-scales, so that we can quantify the accuracy of my time-frequency analysis method.

```Matlab
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
```

<img src="https://user-images.githubusercontent.com/110791799/185438351-6dc453ba-5682-43b9-8c11-7a787396c192.png" alt="obj_fcn" width="500"/>

1- Once we have the time series as a column vector (Note: the code can handle an array of column vectors), the following code creates a wavelet transform surface that can be used to manually choose and separate the time-scales.

```Matlab
dt=(t(2)-t(1)); % Time increment duration (needs to be constant!)
Fs = 2*pi/dt; % Sampling frequency
representative_signal = 1; % The signal with the richest time-frequency content (in case the input is an array of column vectors)
Fo = 2; % Mother Wavelet frequency (Small/large values lead to high temporal/frequency resolution wavelet trasnform surface)
fi = 0; % Lower bound of the frequency range of interest
ff = 10; % Upper bound of the frequency range of interest
nF = 200; % Number of frequencies in the [fi, ff] interval
N_modes = 2; % Number of time-scales (modes) to be extracted from the original signal
WT_pow = .75; % Intensity of low amplitude noise, time-scales and modes (values less/larger than 1 intensify/lessen the low-amplitude details)
WT_plot = "contourf"; % Plotting function for the wavelet transform surface (contourf: accurate but costly, imagesc: less accurate but cheap) 
N_contourf = 30; % resolution of the wavelet transform surface contour plot if WT_plot = "contourf"

% Runs the ICWT code
[t,x_new,interval_freq,modes] = Fast_Morlet_ICWT_Harmonic_decomp(x,representative_signal,Fs,dt,Fo,fi,ff,nF,N_modes,WT_pow,WT_plot,N_contourf);
```

Running the code above gives us the following figure:

<img src="https://user-images.githubusercontent.com/110791799/185442226-bb7e7c69-e81d-4819-95a4-9f41e8a1971b.png" alt="obj_fcn" width="500"/>

After the window with the figure above shows up, we can manually select the time-scales of interest as below:

Time-scale1:

<img src="https://user-images.githubusercontent.com/110791799/185444074-fe91ae51-5493-45f4-be19-87bbd4703194.gif" alt="obj_fcn" width="500"/>

Time-scale2:

<img src="https://user-images.githubusercontent.com/110791799/185444091-b9170493-18cf-4c68-a820-db755dff8551.gif" alt="obj_fcn" width="500"/>

Note: This step can easily be automated if needed.

After this step, we get the time-scales we just selected in time domain. For comparison reason, I am comparing the extracted time-scales (in time and time-frequency domains) with the exact time-scales from the artificial signal:

<img src="https://user-images.githubusercontent.com/110791799/185445045-14d596de-5ca9-424a-944c-c384087183d5.png" alt="obj_fcn" width="750"/>

More importantly, the extracted time-scales are able to reconstruct the original signal with 99.98% accuracy:

<img src="https://user-images.githubusercontent.com/110791799/185445346-673e8c59-8534-4fc7-a507-3ad3f5512f81.png" alt="obj_fcn" width="500"/>

## Further Reading
The details of this method and some of the scenarios and applications it can be implemented in are [here](https://www.sciencedirect.com/science/article/abs/pii/S0888327021000868).

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
