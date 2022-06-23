function [phase]=get_phase1(EEG,TMS_loc,passband,varargin)
% input:
% EEG: single channel epoched EEG
% window: number of points to fit the AR model. if window==-1, use all the
% prestimulus points to fit the AR model.
% TMS_loc: location of TMS from start of the epoch in sec
% passband: 2x1 vector containing lower and higher cutoff frequency
% arOrd is the order of the autoregressive yule-walker model
% edge: edge to remove due to filter edge artifact
% hilbertwindow is the length of the window used for the hilberttransform
% filterOrd is filter order of the digital filter. default: 10

    p = inputParser;
    addParameter(p,'window',-1,...
        @(x) isPositiveInt(x) || all(x==-1.0));
    
%     addParameter(p,'window',-1);
    addParameter(p,'doPlot',false,...
        @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
    addParameter(p,'fvPlot',true,...
        @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
    addParameter(p,'edge',37,...
        @isPositiveInt);
    addParameter(p,'arOrd',13,...
        @isPositiveInt);
    addParameter(p,'hilbertwindow',[],...
        @isPositiveInt);
    addParameter(p,'filterOrd',[],...
        @isPositiveInt);
%     disp(varargin);
    parse(p,varargin{:});
    
    s=p.Results;
    if s.doPlot
        tmpEEG=EEG;
    end
    if EEG.xmax<0
        EEG = pop_select(EEG,'time',[-TMS_loc EEG.xmax]);%prestimulus eeg
    else
        EEG = pop_select(EEG,'time',[-TMS_loc 0]);%prestimulus eeg
    end
    chunk=squeeze(EEG.data);%the 1st dimension=1
    window=s.window;
    if window==-1
        window=length(chunk);
    else
        assert(window<=size(chunk,1),'window=%d is larger then chunk size=%d',window,size(chunk,1))
        chunk=chunk(end-window+1:end,:);
    end

    %chunk : time window x trial
    fs=EEG.srate;
    filterOrd=s.filterOrd;
    if filterOrd>=size(chunk,1)/3
        filterOrd=floor(size(chunk,1)/3)-1;
        warning('filterOrd is too large. Reset edge to %d',filterOrd);
    end
    D = design_phastimate_filter(filterOrd, passband, fs);
    if s.fvPlot
        fvtool(D)
    end
    chunk=filtfilt(D,double(chunk));
%     chunk = bandpass(chunk,passband,fs);
    
    if s.doPlot% plot the BPF signal in the window before TMS pulse
%         disp(size(tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-window+1:round(1+TMS_loc*tmpEEG.srate),:)))

        len=length(tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-window+1,:));
        offset=reshape(tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-window+1,:),1,len)-reshape(chunk(1,:),1,len);
        tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-window+1:round(1+TMS_loc*tmpEEG.srate),:)=chunk+offset;
        pop_eegplot(tmpEEG,1,1,1);
    end
    arOrd = s.arOrd;
    edge = s.edge;
    if 2*edge>window-1
        edge=floor((window-1)/2);
        warning('edge is too long. Reset edge to %d',edge);
    end
    hilbertwindow=s.hilbertwindow;
    if isempty(hilbertwindow)
        forwardsample = edge+63;
    else
        forwardsample = edge + ceil(hilbertwindow/2);
    end
    coeffs = aryule(chunk(edge+1:end-edge,:), arOrd); 
    coeffs = -coeffs(:,end:-1:2);
    nextvalues = zeros(arOrd+forwardsample,size(chunk,2));
    nextvalues(1:arOrd,:) = chunk(end-arOrd-edge+1:end-edge,:);

    for i = 1:forwardsample
        nextvalues(arOrd+i,:) = sum(coeffs'.*nextvalues(i:arOrd+i-1,:),1);
    end
%     figure;
%     y=cat(1,chunk(1:end-arOrd-edge,:),nextvalues);
%     plot(y);
%     xline(arOrd+size(chunk(1:end-arOrd-edge,:),1));
    tmp=nextvalues(arOrd+1:end,:);
    tmp=cat(1,chunk(1:end-edge,:), tmp);%concatenate along time axis
    if isempty(hilbertwindow)
        phase = angle(hilbert(tmp));
        phase = phase(end-forwardsample+edge,:);
    elseif size(tmp,1)<=hilbertwindow
%             fprintf('tmp length=%d\n',size(tmp,1));
        org_hilbertwindow=hilbertwindow;
        hilbertwindow=size(tmp,1);
        warning('Edge=%d and forwardsample=%d .hilbertwindow is too large. Reset hilbertwindow from %d to tmp length(%d)',...
            edge,forwardsample,org_hilbertwindow,hilbertwindow);
        clear org_hilbertwindow;
        phase = angle(hilbert(tmp));
        phase = phase(end-forwardsample+edge,:);
    else
        tmp=tmp(end-hilbertwindow+1:end,:);
        phase = angle(hilbert(tmp));
        phase=phase(end-forwardsample+edge,:);
    end
%     if isempty(hilbertwindow)
%         phase = angle(hilbert(tmp));
%         phase = phase(edge,:);
%     else
%         
%         
%     end
        
    if s.doPlot%plot the AR predicted values
        disp('tmpEEG.data');
        disp(size(tmpEEG.data));
        disp('replaced interval');
        disp(size(nextvalues(arOrd+1:end,:)));
        len=length(tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-edge+1,:));
        offset=reshape(tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-edge+1,:),1,len)-reshape(nextvalues(arOrd+1,:),1,len);
        tmpEEG.data(:,round(1+TMS_loc*tmpEEG.srate)-...
            edge+1:round(1+TMS_loc*tmpEEG.srate)-edge+forwardsample,:)=nextvalues(arOrd+1:end,:)+offset;
        pop_eegplot(tmpEEG,1,1,1);
    end
end