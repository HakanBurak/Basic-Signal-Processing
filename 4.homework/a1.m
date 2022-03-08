clear all
load BSP_P300.mat
fs = 2000; % Hz
% cropping of the data
important_start = 178*1000;
data = data(important_start:end,:);
Ns = size(data,1); % number of samples


BNC = data(:,1);
EEG = data(:,2);
Button = data(:,3);

% time array
t = linspace(0,(Ns-1)/fs,Ns) ;
% frequency array
f = linspace(0,fs/2,round(Ns/2)) ;
% Spectrum of the EEG signal
fEEG = fft(EEG);
fEEG = abs(fEEG(1:round(Ns/2)));
figure
plot(f,fEEG )
xlabel('frequency (Hz)')
ylabel('spectrum')


fc = 50;
[b,a] =  butter(9,fc/(fs/2),'low');
% Inspect the Bode diagram of the filter
figure
freqz(b,a,fs,fs);

% use filtfilt to avoid signal delay
EEG_prefilt = filtfilt(b,a,EEG);
fEEG_prefilt = fft(EEG_prefilt);
fEEG_prefilt = abs(fEEG_prefilt(1:round(Ns/2)));


figure(5)
alpha = mean(fEEG_prefilt(8:14));
beta = mean(fEEG_prefilt(14:30)); 
teta = mean(fEEG_prefilt(4:8));
delta = mean(fEEG_prefilt(1:4));
gama = mean(fEEG_prefilt(30:50));
%y = [alpha beta teta delta gama];
%x = 1:10:50;
%bar(x,y);

plot(f,fEEG_prefilt(1:round(Ns/2)))
xlabel('frequency (Hz)')
ylabel('spectrum')
hold off 
dampened_magnitude = 0.1*filter(b,a,fEEG_prefilt(100));

figure
plot(t,BNC)
xlabel('ttime (s)')
ylabel('V')

%finding all the peaks
[pks, locs] = findpeaks(BNC,'MinPeakHeight',0.2,"MinPeakDistance",0.5*fs);
hold on
plot(t(locs),pks,'k*')

% picking only the higher amplitude, deviant peaks
dev_locs = locs(pks> 0.5);
dev_pks = pks(pks>0.5);
plot(t(dev_locs),dev_pks,'r*')

% picking only the lower amplitude, standard peaks
std_locs = locs(pks< 0.39);
std_pks = pks(pks<0.39);
plot(t(std_locs),std_pks,'g*')


% number of samples in 800 ms-long recording
Np300 = length(EEG_prefilt(1:800)) ;

% averaging the standard EEG responses
EEG_std = zeros(Np300,1);
EEG_std = EEG_prefilt(1:Np300);
sum1 =[];
ind =[];
sum2 =[];
ind2 =[];

% this might actually be multiple lines...
%for i = 1:length(std_locs)
 %   result = find(EEG_prefilt==std_locs(i));
%    sum = [sum,EEG_prefilt(result)];
 %   ind = [ind,result];
%end
for i = 1:length(std_locs)
    result = EEG_prefilt(std_locs(i)) ;
    sum1  = [sum1,result];
end
S = sum(sum1,'all');
avg = S/length(std_locs);

for i = 1:length(std_locs)
    EEG_std(i) = avg;
end
% averaging the deviant EEG responses
EEG_dev = zeros(Np300,1);
EEG_dev = EEG_prefilt(1:Np300);
% this might actually be multiple lines...
 for j = 1:length(dev_locs)
    result1 =EEG_prefilt(std_locs(j));
    sum2 = [sum2,(result1)];
 end
S2 = sum(sum2,'all');
avg2 = S2/length(dev_locs);
for i = 1:length(std_locs)
    EEG_dev(i) = avg2;
end


t_short = linspace(0,Np300/fs,Np300);

figure
plot(t_short*1000,EEG_std)
hold on
plot(t_short*1000,EEG_dev)
legend('normal','deviant')
xlabel('time (ms)')
ylabel('uV')
title('filtered')
