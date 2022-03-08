clear all;

% load the signal as variable 'm'
load 'mod_resp_sig.mat'
% sampling frequency of the downloaded signal
f_s = 10000; % Hz
N = length(m);
% time array of the plot
t = linspace(0,N/f_s,N);
% frequency array of the plot
f = linspace(0,f_s/2,N/2);
% fourier transform of the signal
fourier_m = fft(m);
fourier_m = abs(fourier_m(1:length(fourier_m)/2));
figure
subplot(3,1,1)
plot(t,m)
title('Time domain')
xlabel('time (s)')
ylabel('amplitude')

subplot(3,1,2)
plot(f,fourier_m)
title('Frequency domain')
xlabel('time (Hz)')
ylabel('amplitude')

subplot(3,1,3)
plot(f,fourier_m)
title('Zoomed-in frequency domain')
xlabel('time (Hz)')
ylabel('amplitude')
xlim([1800,2200])
ylim([0,100])