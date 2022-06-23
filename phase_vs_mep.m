addpath('G:\我的雲端硬碟\Documents\110上學期\研究\EMG-Feature-Extraction-Toolbox-master');
% root_folders={'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0227\誠',...
%     'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0312\誠'};
root_folders={'G:\我的雲端硬碟\Documents\110上學期\研究\TMS data\0227\誠'};
phases=[];
zscores=[];% z score of mep
Muscles={'APB','FDI'};
% num_trials=ntrials(root_folders);
for i =1:length(root_folders)
    root_folder=root_folders{i};
%     cd(root_folder);
    myFiles = dir(fullfile(root_folder,'*.set'));
    for k =1:length(myFiles)
        file_name=myFiles(k).name;
        EEG=pop_loadset(file_name,root_folder);
        EEG = pop_epoch( EEG, {  'TMS'  }, [-0.5           1], 'epochinfo', 'yes');
        APB = pop_select(EEG,'channel',{'T8','TP10'});
        APB = pop_reref(APB,'TP10');
        APB = pop_select(APB,'time',[0 0.02]);
        FDI = pop_select(EEG,'channel',{'T7','TP9'});
        FDI = pop_reref(FDI,'TP9');
        FDI = pop_select(FDI,'time',[0 0.02]);
        EEG = pop_select(EEG,'nochannel',{'T8','TP10','T7','TP9'});
        EEG = eeg_checkset( EEG );
        C3 = pop_select(EEG,'channel',{'C3'});
        C3 = pop_select(C3,'time',[-0.5 0]);%prestimulus eeg
%         C3 = pop_eegfilt(C3,8,13);
        eegdata=squeeze(C3.data);
        phase=now_phase(eegdata,C3.srate);
        phases = [phases phase];
        emgdata=APB.data;
        disp('emgdata');
        disp(size(emgdata));
        p2p = squeeze(abs(max(emgdata)-min(emgdata)));
        disp(size(zscore(p2p)));
        zscores = [zscores zscore(p2p)'];
    end
end
figure;
scatter(phases, zscores);
zscores_norm = (zscores-min(zscores))/(max(zscores)-min(zscores));
figure;
subplot(2,1,1);
polarscatter(phases, zscores_norm);
ang = deg2rad(linspace(0,360,361));% angles
vals = cell(1,360);
discrete_phases=mod(floor(rad2deg(phases)),360);
for i = 1:length(zscores)
    vals{1,discrete_phases(i)+1}=[vals{1,discrete_phases(i)+1} zscores(i)];
end
for i = 1:length(vals)
    if isempty(vals{i})
        vals{i}=[inf];
    end
end
disp('vals');
disp(size(vals));
vals = cellfun(@double,vals,'UniformOutput',false);
vals = cellfun(@mean, vals);

vals_norm = (vals-min(zscores))/(max(zscores)-min(zscores)); 
vals_norm(vals_norm == inf) = 0;
subplot(2,1,2);
polarhistogram('BinEdges',ang,'BinCounts',vals_norm);
function [phase]=now_phase(chunk,fs)
%chunk : time points x trial
    chunk = bandpass(chunk,[8 13],fs);
    p = 13;
    edge = 37;
    forwardsample = 100;
    coeffs = aryule(chunk(edge+1:end-edge,:), p); 
    coeffs = -coeffs(:,end:-1:2);
    nextvalues = zeros(p+forwardsample,size(chunk,2));
    nextvalues(1:p,:) = chunk(end-p-edge+1:end-edge,:);

    for i = 1:forwardsample
        nextvalues(p+i,:) = sum(coeffs'.*nextvalues(i:p+i-1,:),1);
    end
%     figure;
%     y=cat(1,chunk(1:end-p-edge,:),nextvalues);
%     plot(y);
%     xline(p+size(chunk(1:end-p-edge,:),1));
    phase = angle(hilbert(nextvalues(p+1:end,:)));
    phase = phase(edge,:);
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

% https://stackoverflow.com/questions/41556971/creating-a-plot-in-polar-coordinates-with-magnitude-vectors
% ang = deg2rad(linspace(0,360,25));% angles
% vals = 1:24; % values
% polarhistogram('BinEdges',ang,'BinCounts',vals)