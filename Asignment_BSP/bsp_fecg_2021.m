function     [fetal_QRSAnn_est,QT_Interval]=bsp_fecg_2021(tm,ecgs,Fs)

% number of samples
%disp(tm);
%initilazing first values
fs = Fs;
N = size(ecgs(:, 1), 1);
t=(0:1/Fs:(length(ecgs)-1)/Fs)';
%to see what we are deailing with
% plot(t,ecgs(:, 1));
% xlim([0 100]);

%using prewhiten data
mixdata = prewhiten(ecgs);
% figure
% histogram(mixdata(:,1))
q = 4;
%Applying ica on the data
Mdl = rica(mixdata,q,'NonGaussianityIndicator',ones(4,1));
%to use Mdl we need to transform it
unmixed = transform(Mdl,ecgs);
variances = [];
%To select correct mother peaks we thought using minumum varianced channel
%would give more accurate results
for i=1:4
    sig=unmixed(:,i);
    [pks, locs2] = findpeaks(abs(sig),'MinPeakProminence',mean(abs(sig))*1.5,"MinPeakDistance",0.4*fs);
    difference = diff(locs2);
    variances = [variances var(difference)];
end
[~, minvar] = min(variances>0);
%motherECG = ecgs(:,minvar);


[pks, locs2] =findpeaks(unmixed(:,minvar),'MinPeakProminence',mean(abs(sig))*1.5,"MinPeakDistance",0.4*fs);
border =80;
fetal_ecg = ecgs;
%Also we used for loop to get rid of mother peaks in ecgs(mixed signals)
%with regards to location points on every channel
for j=1:size(ecgs,2)
    for i=1:length(locs2)
        lower_boundry =locs2(i)-border;
        upper_boundry =locs2(i)+border;
        if lower_boundry < 1
            lower_boundry =1;
        end
        if upper_boundry>length(fetal_ecg)
            upper_boundry = length(fetal_ecg);
        end
        fetal_ecg(lower_boundry:upper_boundry,j) = 0;%we tried changing it with mean however sometimes mean can go over fetal peaks so we changed it with 0 %mean(envelopeABP(locs1(i)-border:locs1(i)+border))
    end
end
%to see if we can get fetal_ecg accurate and we can
% figure
% plot(t,ecgs(:,1));
% hold on;
% plot(t, fetal_ecg(:,1));
% legend('Mixed ECG','Fetal ECG')

%we used Pan Tompkins algorithm which is showed below in practice and in the other assignments  %taking bandpass variable
fetal_variances = [];
%as we did in the ignoring mother ECG part we are going through every
%channel on fetal ECG channel to get fetal variances. So if the fetal ECG
%variance and mother ECG varience is similiar that means we can wrong
%locations while finding fetal peaks.
for i = 1:size(fetal_ecg,2)
    current_fetal = fetal_ecg(:,i);
    [b,a] = butter(3,[6, 42]/(fs/2),'bandpass');
    %filtering current_fetal channel (bandpass filtering)
    filtfetal = filtfilt(b, a, current_fetal);
%     just to see if we can get it correctly
%     plot(t/60, filtfetal);
%     xlim([7.65, 7.85])
%     applying diff oparator on filtered signal
    aDiff= 1;
    bDiff= (1/8)*[2, 1, 0, -1, -2];
    difffetal = filtfilt(bDiff, aDiff, filtfetal);
%     to check
%     plot(t/60, difffetal)
%     getting the envolope to select peaks more easily.
    bEnvelope = ones(1, 30)/30;
    aEnvelope = 1;
    envelopefetal = filtfilt(bEnvelope, aEnvelope, difffetal.^2);
    %above we used Pan Tompkins method for getting peaks more accuracetly. We selected
    %Minpeak distance smaller then mother because fetal's BPM is much more from mothers
    [pks1, locs_current] = findpeaks(envelopefetal,'MinPeakProminence',mean(abs(envelopefetal))*2,"MinPeakDistance",0.2*fs);
    difference = diff(locs_current);
    fetal_variances = [fetal_variances var(difference)];
end
%getting the least variance and not same with mothers variances
min_fetal_ecg = fetal_variances(1);
for i=1:size(fetal_variances)
    if fetal_variances(i) < min_fetal_ecg && abs(min_fetal_ecg-variances(minvar)) > 10^-4
        min_fetal_ecg = fetal_variances(i);
    end
    
end
%we select those points to get going
min_fetal_loc = find(min_fetal_ecg ==fetal_variances);
%there was some error with min_fetal_loc and we couldn't figure it out so
%instead we said if you directly use the second one it works better than any other.
current_fetal = fetal_ecg(:,2);%min_fetal_loc

%applying filter
filtfetal = filtfilt(b, a, current_fetal);
% to check if can separate those 2 signal
 plot(t/60, filtfetal);
 xlim([7.65, 7.85])
legend('Fetal signal after filtering')
% arranging diff values
aDiff= 1;
bDiff= (1/8)*[2, 1, 0, -1, -2];
% applying dif on our fetal signal
difffetal = filtfilt(bDiff, aDiff, filtfetal);
% plotting diff we can get it correct

 plot(t/60, difffetal)
legend('Fetal signal after diff')
 xlim([7.65, 7.85])
% envolope variables
bEnvelope = ones(1, 30)/30;
aEnvelope = 1;
%getting envelope of the fetal to get peaks
envelopefetal = filtfilt(bEnvelope, aEnvelope, difffetal.^2);
%  plot(t/60, envelopefetal)
% legend('Envolope of Fetal signal')
% xlim([7.65, 7.85])




[pks1, fetal_peaks] = findpeaks(envelopefetal,'MinPeakProminence',mean(abs(envelopefetal))*2,"MinPeakDistance",0.2*fs);
temp =fetal_ecg(:,2);
figure(1)
subplot(4, 1, 1)
plot(t/60, filtfetal)
xlim([7.65, 7.85])
title("fetal after bandpass filtering")
subplot(4, 1, 2)
plot(t/60, difffetal)
xlim([7.65, 7.85])
title("fetal after differentiation")
subplot(4, 1, 3)
plot(t/60, envelopefetal)
xlim([7.65, 7.85])
title("fetal- envelope")
subplot(4, 1, 4)
plot(t/60, temp)
hold on
plot(t(fetal_peaks)/60, temp(fetal_peaks), 'g*')
xlim([7.65, 7.85])
title("fetal with systolic peaks")



%this was to understand if it works correctly and it does
figure(2)
plot(t/60,ecgs(:,2),'b');
hold on
plot(t/60, fetal_ecg(:,2),'r');
plot(t(fetal_peaks)/60,fetal_ecg(fetal_peaks,2),'g*')
xlim([7.66 7.8])
title('Fetal peaks')
legend('Mixed ECG','Fetal ECG', 'Fetal Peaks')
%assigning relevant variables to output variables
fetal_QRSAnn_est=fetal_peaks;
QT_Interval=(t(fetal_peaks(2:end))-t(fetal_peaks(1:end-1)));
end

%the codes below unsued but tried multiple times(days) 


% first = unmixed(:,2);
% max_first = max(first);
% second = unmixed(:,1);
% max_second = max(second);
% if max_first < max_second
%     fetal = first;
%     mother = second;
% else
%     fetal = second;
%     mother = first;
% end
% disp(max(fetal));
% % disp(max(mother));
% %1/tm;
% fetal = ecgs;
% [b,a] = butter(3,[1, 10]/(fs/2),'bandpass');
%
% filtfetal = filtfilt(b, a, fetal);
% plot(t/60, filtfetal);
% xlim([7.65, 7.85])
% aDiff= 1;
% bDiff= (1/8)*[2, 1, 0, -1, -2];
% difffetal = filtfilt(bDiff, aDiff, filtfetal);
% plot(t/60, difffetal)
% bEnvelope = ones(1, 30)/30;
% aEnvelope = 1;
% envelopefetal = filtfilt(bEnvelope, aEnvelope, difffetal.^2);
%
% plot(t/60, envelopefetal)
% xlim([7.65, 7.85])
%
% [pks1, locs1] = findpeaks(envelopefetal,'MinPeakHeight',5*10^(-8),"MinPeakDistance",0.2*fs);
%
% std_locs = locs(pks1 < 10^(-7));
% std_pks = pks(pks1 < 10^(-7));
% plot(t(std_locs)/60,std_pks,'r*')
% border =25;
% for i=1:length(locs2)
%     lower_boundry =locs2(i)-border;
%     upper_boundry =locs2(i)+border;
%     if lower_boundry < 1
%         lower_boundry =1;
%     end
%     if upper_boundry>length(envelopefetal)
%     upper_boundry = length(envelopefetal);
%     end
%     envelopefetal(lower_boundry:upper_boundry) = 0;%mean(envelopeABP(locs1(i)-border:locs1(i)+border))
% end
% plot(t/60, envelopefetal)
% xlim([7.670, 7.685])
%
% subplot(5,1,1)
% plot(t/60,fetal);
% hold on
%
% [pks, locs2] = findpeaks(envelopefetal,'MinPeakHeight',1*10^(-8),"MinPeakDistance",0.2*fs);
%
% subplot(5,1,5)
% plot(t/60,fetal);
% hold on
% plot(t(locs2)/60,fetal(locs2),'g*')
% xlim([7.66 7.8])
% title('peaks')
% plot(t/60, fetal)
% hold on
% plot(t(locs2)/60, fetal(locs2), 'g*')
% xlim([7.65, 7.85])
% title("peaks")
% %plotting part
%
% subplot(5, 1, 2)
% plot(t/60, filtfetal)
% xlim([7.65, 7.85])
% title("fetal after bandpass filtering")
% subplot(5, 1, 3)
% plot(t/60, difffetal)
% xlim([7.65, 7.85])
% title("fetal after differentiation")
% subplot(5, 1, 4)
% plot(t/60, envelopefetal)
% xlim([7.65, 7.85])
% title("fetal- envelope")
% subplot(5, 1, 5)
% plot(t/60, fetal)
% hold on
% plot(t(locs1)/60, fetal(locs1), 'g*')
% xlim([7.65, 7.85])
% title("fetal with systolic peaks")

% disp(size(unmixed));
% figure
% subplot(1,2,1);
% plot(time, unmixed(:,1))
% xlim([0 100]);
%
% mecg=ecgs(:,1);
%
%
% xlim([0 100]);
% subplot(1,2,2);
% plot(time, unmixed(:,2))
%[min_vol,min_loc]=findpeaks(mecg,'MinPeakDistance',Fs*0.6);


%
% 