clear all;

% load the signal as variable 'm'
load 'mod_resp_sig.mat'
% sampling frequency of the downloaded signal
f_s = 10000; % Hz
N =length(m);
T = 1/f_s;
% time array of the plot
t =  (0:N-1)*T;
% frequency array of the plot
f =  (f_s*(0:(N/2))/N);
% fourier transform of the signal
fourier_m = fft(m) ;
P2 = abs(fourier_m/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
samplerange =linspace(1800,2200);
figure
subplot(3,1,1)
plot(t,m )
title('Time domain')
xlabel('time (s)')
ylabel('amplitude')

subplot(3,1,2)
plot(f,P1)
title('Frequency domain')
xlabel('time (Hz)')
ylabel('amplitude')

subplot(3,1,3)
plot(f,P1); 
title('Zoomed-in frequency domain')
xlabel('time (Hz)')
ylabel('amplitude')
xlim([1800,2200])
ylim([0,100])