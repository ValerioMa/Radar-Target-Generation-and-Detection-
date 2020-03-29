%% Compute the FFT of a generated signal
% Generate a signal containing 2 sinusoidal wave of different frequencies
% and compute the FFT

%% Init console
close all force;
clc;

%% Params
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

%% Signal generation
% Form a signal containing a 77 Hz sinusoid of amplitude 0.7 and a 43Hz sinusoid of amplitude 2.
S = 0.7*sin(2*pi*77*t) + 2*sin(2*pi*43*t);

% Corrupt the signal with noise 
X = S + 2*randn(size(t));

% Plot the noisy signal in the time domain. It is difficult to identify the frequency components by looking at the signal X(t). 
f_id = figure(1);
subplot(2,1,1); hold on; grid on; box on;
plot(1000*t(1:50) ,X(1:50))
plot(1000*t(1:50) ,S(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
legend('with noise','without noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

%% FFT
% Compute the Fourier transform of the signal. 
Y = fft(X);    % noisy fft
Ys = fft(S);   % no noise fft

% Take the ampliture of the normalized signal
P2 = abs(Y/L);
P2s = abs(Ys/L);

% TODO : Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P1  = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1s  = P2s(1:L/2+1);
P1s(2:end-1) = 2*P1s(2:end-1);

% Plotting
figure(f_id)
subplot(2,1,2); hold on; grid on; box on;
f = Fs*(0:(L/2))/L;
plot(f,P1) 
plot(f,P1s) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('with noise','without noise')