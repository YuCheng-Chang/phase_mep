function [true_phases]=get_true_phases(EEG,passband,varargin)
% input
% EEG: single channel(C3) continuous EEG
% passband: 2x1 vector containing lower and higher cutoff frequency
    p = inputParser;
    addParameter(p,'doPlot',false,...
        @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
    parse(p,varargin{:});
    s=p.Results;
    EEG = pop_eegfiltnew(EEG, 'locutoff',passband(1),'hicutoff',passband(2),'plotfreqz',0);
    if s.doPlot
        pop_eegplot(EEG,1,1,1);
    end
%     title('true');
    tmp=angle(hilbert(squeeze(EEG.data)));
    true_phases=nan(1,length(EEG.event)-1);%exclude boundary event
    for iEvent=2:length(EEG.event)
        true_phases(iEvent-1)=tmp(EEG.event(iEvent).latency );
    end
end
