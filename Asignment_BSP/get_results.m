close all;
clear all;

set_name='set-ar'; % (change this to use other datasets)
rec_ext='dat'; % (using the WFDB binary dataset)
records=dir([set_name '\*.' rec_ext]);
I=length(records);
display(['Processing ' num2str(I) ' records ...'])

for i=1:I
    record_id=records(i).name(1:end-4);
    disp(records(i));
    fname=records(i).name;  
    %[tm,ecgs]=rdsamp([set_name '\' record_id]);
    [ecgs, Fs, tm]=rdsamp([set_name '\' record_id]);
    fid = fopen([set_name '\' record_id '.hea']);
    text = fgetl(fid);
    wsp=strfind(text,' ');
    Fs=str2num(text(wsp(2):wsp(3)));
    [fetal_QRSAnn_est,QT_Interval]=bsp_fecg_2021(tm,ecgs,Fs);
    disp(QT_Interval);
    
    plot_ecgs(Fs,ecgs,fetal_QRSAnn_est,QT_Interval);
end    