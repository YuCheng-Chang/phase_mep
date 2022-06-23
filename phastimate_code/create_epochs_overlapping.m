function [epochs, truephase] = create_epochs_overlapping(signal, fs,nr_seconds,passband)
addpath('G:\我的雲端硬碟\Documents\110上學期\研究\phase_mep');
% disp('create_epochs_overlapping');
% disp('signal');
% disp(size(signal));
assert(nr_seconds>0,'buffer length should be positive');
if isempty(nr_seconds)
    nr_seconds = 2;
end
nr_samples = nr_seconds * fs;
epoch_overlap = 0.75*nr_samples; 
[epochs, ~] = buffer(signal, nr_samples, epoch_overlap,'nodelay');
[epochs2,~] = buffer(signal,nr_samples*2,epoch_overlap+nr_samples,'nodelay');
epochs=epochs(:,1:size(epochs2,2));
 % set-up equivalent filter objects for given peak frequency
filter_objects = {};
mid_frequency=mean(passband,'all');
for ord = [2 3 4 5] % FIR - windowed sinc
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/mid_frequency)), 'CutoffFrequency1', passband(1)-1, 'CutoffFrequency2', passband(2)+1, 'SampleRate', fs, 'DesignMethod', 'window')};
end
for ord = [3 4 5] % FIR - least squares (equiripple is similar)
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/mid_frequency)), 'StopbandFrequency1', passband(1)-4, 'PassbandFrequency1', passband(1)-1, 'PassbandFrequency2', passband(2)+1, 'StopbandFrequency2', passband(2)+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
end
for ord = [4 8 12] % IIR - butterworth
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', passband(1), 'HalfPowerFrequency2', passband(2), 'SampleRate', fs, 'DesignMethod', 'butter')};
end
for ord = [4 6 8] % IIR - chebychev I
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', passband(1)-1.5, 'PassbandFrequency2', passband(2)+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
end
for attenuation = [10 20] % IIR - elliptic
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', passband(1)-2, 'PassbandFrequency1', passband(1)-1, 'PassbandFrequency2', passband(2)+1, 'StopbandFrequency2', passband(2)+2, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
end    
[truephase, ~, ~, ~] = phastimate_truephase(epochs2, filter_objects);
% y=buffer(x,n); The number of columns in y is ceil(L/(n-p)).L is the
% length of x
% tmp=zeros([1,size(epochs2)]);
% tmp(1,:,:)=epochs2;
% nr_epochs=size(tmp,3);
% latency=num2cell(size(epochs2,1)*(1:nr_epochs));
% EEG = pop_importdata('dataformat','array','nbchan',1,'data','epochs2','srate',fs,'pnts',size(tmp,2),'xmin',-nr_seconds);
% EEG.event=struct('latency',latency,'type','TMS','chan',[],'epoch',num2cell(1:nr_epochs));
% EEG.epoch=struct('event',num2cell(1:nr_epochs),'eventlatency',0,'eventtype','TMS','eventchan',[]);
% pop_eegplot(EEG,1,1,1);
% % remove some starting samples
% epochs = epochs(:, ceil(nr_samples/(nr_samples-epoch_overlap)):end);
% 
% % create timeaxis in milliseconds
% diff_t = 1000/fs; 
% time = 0:diff_t:(nr_samples-1)*diff_t;
% time = time-time(floor(numel(time)/2));

% linear detrending
epochs = detrend(epochs, 'constant');
disp('epochs');
disp(size(epochs));

    
end