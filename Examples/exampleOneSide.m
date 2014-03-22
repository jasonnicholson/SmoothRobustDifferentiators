% clean up
clc;close all; clear;
%% Signal Parameters
sampleRate = 1000; % Hz
dt = 1/sampleRate;
t = 0:dt:1;
signalFrequency = 40; % Hz
noiseMeanFrequency = 350; % Hz
noiseMeanAmplitude = 0.2;

%% Calculate Signal
pureSignal = sin(2*pi*signalFrequency*t);

% Noise has an amplitude range of 75% to 125% of mean amplitude
% Noise varies in frequency from 90% to 110% of mean frequency
% Noise has a phase delay of 0.3rad/s
noisySignal = sin(2*pi*signalFrequency*t) + noiseMeanAmplitude*(0.75+0.5*rand(size(t))).*cos(2*pi*noiseMeanFrequency*(0.9+0.1*rand(size(t))).*t-0.3);

%% Calculate Derivative
signal_dot = 2*pi*signalFrequency*cos(2*pi*signalFrequency*t);

signal_dot_estimate5 = robustDiffOneSide(noisySignal,dt,5);
signal_dot_estimate6 = robustDiffOneSide(noisySignal,dt,6);
signal_dot_estimate9 = robustDiffOneSide(noisySignal,dt,9);
signal_dot_estimate10 = robustDiffOneSide(noisySignal,dt,10);
signal_dot_estimate15 = robustDiffOneSide(noisySignal,dt,15);
signal_dot_estimate21 = robustDiffOneSide(noisySignal,dt,21);

%% Plots Signals
subplot(211)
plot(t,pureSignal,t,noisySignal);
xlim([0 4/signalFrequency])
xlabel('time [sec]')
ylabel('Amplitude')
title('y')
legend('Pure Signal','Noisy Signal');

%% Plots derivative of Signals
% Please note the one sided formulas have a significant phase delay.  This
% means you should take this into account for your signal.
subplot(212)
plot(t, signal_dot, ...
     t, signal_dot_estimate5, ...
     t, signal_dot_estimate6, ...
     t, signal_dot_estimate9, ...
     t, signal_dot_estimate10, ...
     t, signal_dot_estimate15, ...
     t, signal_dot_estimate21);
xlim([0 4/signalFrequency])
xlabel('time [sec]')
ylabel('Amplitude')
title('y''')
legend('Pure Signal Derivative', ...
       'N = 5', ...
       'N = 6', ...
       'N = 9', ...
       'N = 10', ... 
       'N = 15', ...
       'N = 21');
disp(['Please note the one sided formulas have a significant phase delay.  This ' ...
      'means you should take this into account for your signal.']);