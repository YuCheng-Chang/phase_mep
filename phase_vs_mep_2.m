addpath('G:\我的雲端硬碟\Documents\110上學期\研究\EMG-Feature-Extraction-Toolbox-master');
addpath(genpath('G:\我的雲端硬碟\Documents\110上學期\研究\AARATEPPipeline\Common'))
% root_folders={'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0227\誠',...
%     'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0312\誠',...
%     'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0409\誠',...
%     'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0415\誠'};
root_folders={'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0409\誠',...
    'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0415\誠'};
phases=[];
zscores=[];% z score of mep
Muscles={'APB','FDI'};
s.artifactTimespan= [-0.002, 0.8];
s.doDecayRemovalPerTrial=true;
s.bandpassFreqSpan=[8 13];
s.epochTimespan=[-0.5 1];
% num_trials=ntrials(root_folders);
for i =1:length(root_folders)
    root_folder=root_folders{i};
%     cd(root_folder);
    myFiles = dir(fullfile(root_folder,'*.set'));
    for k =1:length(myFiles)
        file_name=myFiles(k).name;
        EEG=pop_loadset(file_name,root_folder);
%         EEG = pop_epoch( EEG, {  'TMS'  }, [-0.5           1], 'epochinfo', 'yes');
        APB = pop_select(EEG,'channel',{'T8','TP10'});
        APB = pop_eegfiltnew(APB, 'locutoff',10,'hicutoff',[]);
        APB = pop_eegfiltnew(APB, 'locutoff',58,'hicutoff',62,'revfilt',1,'plotfreqz',0);
        APB = pop_reref(APB,'TP10');
        APB = pop_epoch( APB, {  'TMS'  }, s.epochTimespan, 'epochinfo', 'yes');
        APB = pop_select(APB,'time',[0 0.05]);

        C3 = pop_select(EEG,'channel',{'C3'});
        C3 = pop_epoch( C3, {  'TMS'  }, s.epochTimespan, 'epochinfo', 'yes');
        C3 = c_TMSEEG_applyModifiedBandpassFilter(C3,...
            'lowCutoff', s.bandpassFreqSpan(1),...
            'highCutoff',s.bandpassFreqSpan(2),...
            'artifactTimespan', s.artifactTimespan*3);
        pop_eegplot(C3,1,1,1)
%         C3 = pop_eegfiltnew(C3, 'locutoff',58,'hicutoff',62,'revfilt',1,'plotfreqz',0);
%         C3 = c_EEG_ReplaceEpochTimeSegment(C3,...
%             'timespanToReplace', artifactTimespan,...
%             'eventType','TMS',...
%             'method', 'zero');
%         C3 = pop_eegfiltnew(C3, 'locutoff',58,'hicutoff',62,'revfilt',1,'plotfreqz',0);
        C3 = c_EEG_ReplaceEpochTimeSegment(C3,...
            'timespanToReplace', s.artifactTimespan,...
            'method', 'ARExtrapolation',...
            'prePostFitDurations', [20 20]*1e-3);

%         
%         phase=now_phase(squeeze(C3.data),C3.srate);
        pop_eegplot(C3,1,1,1)
        phases = [phases get_phase2(C3,-s.epochTimespan(1))];
        emgdata=APB.data;
        disp('emgdata');
        disp(size(emgdata));
        p2p = squeeze(abs(max(emgdata)-min(emgdata)));
        disp(size(zscore(p2p)));
        zscores = [zscores zscore(p2p)'];
    end
end


%assinging 15 intervals [0 24]、 [24 48]、...[336 360]
n_interval=15;
ang = deg2rad(linspace(0,360,n_interval+1));% angles
vals = cell(1,n_interval);
i_interval=floor(wrapTo360(rad2deg(phases))/(360/n_interval))+1;%from 1 to n_interval
for i = 1:length(zscores)
    vals{1,i_interval(i)}=[vals{1,i_interval(i)} zscores(i)];
end
disp('vals');
disp(size(vals));
for i = 1:length(vals)
    if isempty(vals{i})
        vals{i}=[inf];
    end
end
vals = cellfun(@double,vals,'UniformOutput',false);
vals = cellfun(@mean, vals);

vals_norm = (vals-min(zscores))/(max(zscores)-min(zscores)); 
vals_norm(vals_norm == inf) = 0;


figure;
scatter(phases, zscores);
hold on;
fo = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[1 0 1 0]);
ft = fittype('a*cos(b*x+c)+d');
%https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
[curve,gof] = fit(phases',zscores',ft,fo)
plot(curve);
bar_interval =  360/n_interval;
bar_x = linspace(0,360,n_interval+1);% angles
bar_x = bar_x(1:end-1);
bar_x = bar_x+bar_interval/2;
bar_x = deg2rad(bar_x);
%wrap btw -pi to pi
for i=1:length(bar_x)
    if bar_x(i)>=pi
        bar_x(i)=bar_x(i)-2*pi;
    end
end
bar_y = vals;
bar_y(bar_y==inf)=0;
bar(bar_x,bar_y);
xlabel('phase');
ylabel('z_{mep}');
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2', '0', '\pi/2', '\pi'});
hold off;

zscores_norm = (zscores-min(zscores))/(max(zscores)-min(zscores));
figure;
subplot(2,1,1);
polarscatter(phases, zscores_norm);
subplot(2,1,2);
polarhistogram('BinEdges',ang,'BinCounts',vals_norm);



function [phase]=get_phase(EEG)
%input:
%EEG: single channel continuous EEG
%     EEG=pop_eegfiltnew(EEG,'locutoff',8,'hicutoff',13);
    nevent=length(EEG.event);
    phases=angle(hilbert(squeeze(EEG.data)));
    phase=[];
    for i=1:nevent
        if strcmp(EEG.event(i).type,'TMS')
            phase=[phase phases(EEG.event(i).latency)];
        end
    end
end
function [phase]=get_phase2(EEG,TMS_loc)
% input:
% EEG: single channel epoched EEG
% TMS_loc: location of TMS from start of the epoch in sec
%     EEG=pop_eegfiltnew(EEG,'locutoff',8,'hicutoff',13);
    
    phases=angle(hilbert(squeeze(EEG.data)));
    phase=phases(TMS_loc.*EEG.srate,:);
end
function [ntrials]=ntrials(root_folders)
    ntrials=0;
    for i =1:length(root_folders)
        root_folder=root_folders{i};
        myFiles = dir(fullfile(root_folder,'*.set'));
        for k =1:length(myFiles)
            file_name=myFiles(k).name;
            EEG=pop_loadset(file_name,root_folder);
            disp(EEG.trials)
            ntrials = ntrials + EEG.trials;
        end
    end
end

% TODO: test the phase estimation algorithm to see the phase distribution
% https://github.com/mastaneht/SPIS-Resting-State-Dataset