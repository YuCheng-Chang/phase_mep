% use a family of equivalent filters to calculate the true phase
function [true_phases]=get_true_phases2(EEG,passband,varargin)
% input
% EEG: single channel(C3) continuous EEG
% passband: 2x1 vector containing lower and higher cutoff frequency
%     p = inputParser;
%     addParameter(p,'doPlot',false,...
%         @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
%     parse(p,varargin{:});
%     s=p.Results;
    EEG = pop_epoch( EEG, {  'TMS'  }, [-1 1], 'epochinfo', 'yes');
    epochs=squeeze(EEG.data);
    % set-up equivalent filter objects for given peak frequency
    filter_objects = {};
    fs=EEG.srate;
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
    [true_phases, ~, ~, ~] = phastimate_truephase(epochs, filter_objects);
    
%     EEG = pop_eegfiltnew(EEG, 'locutoff',passband(1),'hicutoff',passband(2),'plotfreqz',0);
%     if s.doPlot
%         pop_eegplot(EEG,1,1,1);
%     end
% %     title('true');
%     tmp=angle(hilbert(squeeze(EEG.data)));
%     true_phases=nan(1,length(EEG.event)-1);%exclude boundary event
%     for iEvent=2:length(EEG.event)
%         true_phases(iEvent-1)=tmp(EEG.event(iEvent).latency );
%     end
end
