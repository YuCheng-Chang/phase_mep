function [phase]=get_phase2(EEG,TMS_loc,passband,artifactTimespan,varargin)
% input:
% EEG: single channel epoched EEG
% TMS_loc: location of TMS from start of the epoch in sec
% passband: 2x1 vector containing lower and higher cutoff frequency
    p = inputParser;
    addParameter(p,'doPlot',false,...
        @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
    parse(p,varargin{:});
    s=p.Results;
    
    EEG = c_TMSEEG_applyModifiedBandpassFilter(EEG,...
        'lowCutoff', passband(1),...
        'highCutoff',passband(2),...
        'artifactTimespan', artifactTimespan*3,...
        'timeToExtend',1,...
        'doDebug',true,...
        'iCh',1,...
        'iTr',1);
    if s.doPlot
        pop_eegplot(EEG,1,1,1);
    end
%     title('m_BPF');
    EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
        'timespanToReplace', artifactTimespan,...
        'method', 'ARExtrapolation',...
        'prePostFitDurations', [20 20]*1e-3);
    if s.doPlot
        pop_eegplot(EEG,1,1,1);
    end
    phases=angle(hilbert(squeeze(EEG.data)));
    phase=phases(round(1+TMS_loc.*EEG.srate),:);
end