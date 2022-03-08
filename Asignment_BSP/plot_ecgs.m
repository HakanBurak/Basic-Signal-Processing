function plot_ecgs(Fs,ecgs,fetal_QRSAnn_est,QT_Interval)

figure()
len=size(ecgs,2);
[ha] = tight_subplot(len,1,[.03 .03],[.08 .01],[.03 .03]);

for i=1:len
    ecg=ecgs(:,i);
    axes(ha(i));
    t=(0:1/Fs:(length(ecg)-1)/Fs)';
    plot(t,ecg,t(fetal_QRSAnn_est),ecg(fetal_QRSAnn_est),'r*');
    %xlim([7.65, 7.85])

    
    sQT_Interval=sort(QT_Interval);
    avg_dp=mean(sQT_Interval(1:round(end*70/100)));
    mean_FHR=round(1/avg_dp*60);
    legend(['FHR ',num2str(mean_FHR),' bpm']);
end
xlabel('Time (s)')