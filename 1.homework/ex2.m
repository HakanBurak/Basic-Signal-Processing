signal = load("tekli.mat");
    signal = signal.val;
    % number of samples
    N = size(signal);
    % number of channels
    Ch = size(signal,1);
    % sampling  - read from the .info file
    fs = 500;
    % base(s) and gain(s) - copy-paste data from .info file, for all channels if
    % necessary
    base = 759.926000705;
    gain = -3475;
    units = 'uV';
    channel_names = ['EEG Fp1'];
    
    %base = [759.926000705,665.572325678, 573.024187942	,631.69163231	,789.269357091	,390.107862279, 222.987346527,	547.763543346, 527.555401568, 656.48690402, 560.753949022, 469.788751694,  496.93170171, 577.685346555, 197.322543698, 437.220424275, 680.442260993, 624.598992865, 202.349833174, 1180.34654117, 47502.3484796, 32767.5];
        %gain = [-3475, -3016, -513, -2019, -2530, -6292, -1361, -456, -1115, 2539, 2884, -1545, -3431, 990, -2476,	-3611,	943,	1209,	-2092,	4558,	-10770,-1];
        % unit(s) - copy-paste data from .info file, for all channels if necessary
        %units = {'uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','uV','mV','nd'};
        % channel names - copy-paste data from .info file, for all channels if
        % necessary
        %channel_names = {'EEG Fp1',	'EEG Fp2', 'EEG F3', 'EEG F4', 'EEG F7', 'EEG F8', 'EEG T3', 'EEG T4', 'EEG C3', 'EEG C4', 'EEG T5', 'EEG T6', 'EEG P3', 'EEG P4', 'EEG O1', 'EEG O2', 'EEG Fz', 'EEG Cz', 'EEG Pz'	'EEG A2-A1', 'ECG ECG','EDF Annotations'};
        
        % create the appropriate time array, starting from 0s
        time = 0:0.002:59.998;
        
        figure(2)
        for I = 1: Ch
        %correct for base and gain
        signal = (signal- base)/gain;
            
        subplot(Ch,1,I)
        % plotting the upcoming channel data
        plot(time,signal(I,:))
        % title and axes
        title(channel_names)
        xlabel('time (s)')
        ylabel(units)
        end
% create the frequency array you will plot onto. It should start at 0 Hz
    % and end at half of the sampling frequency.
    frequency = 0:0.0083333333:249.992666666666;
    
    % initialize the spectrum variable
    spectra = [];
    
    figure(3)
    for I = 1: Ch
        % create the spectrum of your signal channel, and append it to the
        % previous spectra
        spectrum = fft(signal);
        spectra = [spectra; spectrum];
        
        subplot(Ch,1,I)
        % plotting the upcoming channel data. Be aware: you should have the
        % same number of samples in the frequency and spectrum arrays.
        plot(frequency,spectra(I,:))
        % title and axes
        title(['Spectrum of ',channel_names(I)])
        xlabel('frequency (Hz)')
        ylabel('amplitude')
    end
   