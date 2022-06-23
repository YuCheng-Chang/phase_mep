%% variable declaration
close all;
clear all;
clc;
addpath(genpath('G:\我的雲端硬碟\Documents\110上學期\研究\AARATEPPipeline\Common'));
addpath('G:\我的雲端硬碟\Documents\110上學期\研究\phase_mep\phastimate_code');
root_folder="G:\我的雲端硬碟\Documents\110上學期\研究\phase_mep\SPIS-Resting-State-Dataset-master\Pre-SART EEG";
addpath(root_folder);
s.artifactTimespan= [-0.002, 0.2];
s.doDecayRemovalPerTrial=true;
s.bandpassFreqSpan=[12 13];
s.epochTimespan=[-2 1];
myFiles = dir(fullfile(root_folder,'*.mat'));
true_phases=[];% true phases for algo1 and algo2.
true_phases_algo3=[];% true phases for algo3 ('EO' only).
phases1=[];% estimated phases from algo1.
phases2=[];% estimated phases from algo2.
phases3=[];% estimated phases from algo2.
opt_params=containers.Map;
mat_path=fullfile('G:\我的雲端硬碟\Documents\110上學期\研究\phase_mep\opt_params',...
    sprintf('%d_%d',s.bandpassFreqSpan(1),s.bandpassFreqSpan(2)),...
    'opt_params.mat');
rng(1);

%% optimized parameters
for k =1:length(myFiles)
    % from mat file to EEG structure
    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',fullfile(myFiles(k).folder,myFiles(k).name),'srate',256,'pnts',0,'xmin',0);
    EEG = add_pseudo_TMS(EEG,3,0.2);
    channelList = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';...
        'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz'; 'Fp2';'AF8';'AF4';'Afz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';...
        'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2';'EOG1';'EOG2';'EOG3';'tirgger'};
    chanlocs = struct('labels', channelList);
    EEG.chanlocs=chanlocs;
    C3 = pop_select(EEG,'channel',{'C3'});
    if contains(myFiles(k).name,'_EC') 
        [epochs,truephase] = create_epochs_overlapping(squeeze(C3.data), C3.srate,-s.epochTimespan(1),s.bandpassFreqSpan); % from continuous data
        tmp=zeros(1, size(epochs,1), size(epochs,2));
        tmp(1,:,:)=epochs;
        epochs=tmp;%(1,pnts,epochs)
        clear tmp;
        
        EEG = pop_importdata('dataformat','array','nbchan',size(epochs,1),...
            'data','epochs','srate',C3.srate,'pnts',size(epochs,2),'xmin',s.epochTimespan(1));
        
        
        bounds_window = [10 floor(-s.epochTimespan(1)*C3.srate)-1];%[10 511]
        filter_order_range = 100:floor(-s.epochTimespan(1)*C3.srate/3);%[100 170]
        filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
        for ord = filter_order_range
            filter_objects_by_order{ord} = design_phastimate_filter(ord, s.bandpassFreqSpan, EEG.srate);
        end
        bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
        
        bounds_edge = [5 120];
        bounds_ar_order = [5 60];
        bounds_hilbertwindow=[50 110];
        
        % phastimate_optimize(EEG,passband, truephase, filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, hilbertwindow)
        [optimal_parameters, ga_output] = phastimate_optimize(EEG,s.bandpassFreqSpan, truephase, filter_objects_by_order,...
            bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, bounds_hilbertwindow);
       opt_params(file2subject(myFiles(k).name))=optimal_parameters;
    end
end
clear optimal_parameters;
save(mat_path,"opt_params");
%% phase distribution

load(mat_path);
for k =1:length(myFiles)
    % from mat file to EEG structure
    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',fullfile(myFiles(k).folder,myFiles(k).name),'srate',256,'pnts',0,'xmin',0);
    EEG = add_pseudo_TMS(EEG,3,0.2);
    channelList = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';...
        'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz'; 'Fp2';'AF8';'AF4';'Afz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';...
        'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2';'EOG1';'EOG2';'EOG3';'tirgger'};
    chanlocs = struct('labels', channelList);
    EEG.chanlocs=chanlocs;
    C3 = pop_select(EEG,'channel',{'C3'});
    % true phase
%     true_phases=[true_phases get_true_phases(C3,s.bandpassFreqSpan,'doPlot',true)];
    phases_to_append=get_true_phases2(C3,s.bandpassFreqSpan);
    true_phases=[true_phases phases_to_append];
    
    % epoching
    C3 = pop_epoch( C3, {  'TMS'  }, s.epochTimespan, 'epochinfo', 'yes');
    % phase est algo1.
    phases1 = [phases1 get_phase1(C3,-s.epochTimespan(1),s.bandpassFreqSpan,'window',-1,'doPlot',true,'hilbertwindow',100)];
    % phase est algo2.
    phases2 = [phases2 get_phase2(C3,-s.epochTimespan(1),s.bandpassFreqSpan,s.artifactTimespan,'doPlot',true)];
    if contains(myFiles(k).name,'_EO')
        subject=file2subject(myFiles(k).name);
        optimal_parameters=opt_params(subject);
        true_phases_algo3 = [true_phases_algo3 phases_to_append];
        phases3 = [phases3 get_phase1(C3,-s.epochTimespan(1),...
            s.bandpassFreqSpan,...
            'window',optimal_parameters.window_length,...
            'filterOrd',optimal_parameters.filter_order,...
            'edge', optimal_parameters.edge,...
            'arOrd', optimal_parameters.ar_order,...
            'hilbertwindow',optimal_parameters.hilbertwindow)];


    end
end
%% plot histogram
figure;
subplot(3,1,1);
title('phase distribution from algo1');
polarhistogram(phases1,36);
subplot(3,1,2);
title('phase distribution from algo2');
polarhistogram(phases2,36);
subplot(3,1,3);
title('phase distribution from algo3');
polarhistogram(phases3,36);
figure;
title('true phase distribution');
polarhistogram(true_phases,36);
%% calculate average phase error of 2 different algo.
disp('algo1');
disp(avg_error(true_phases,phases1))
% disp(mean(abs(true_phases-phases1)));
disp('algo2');
disp(avg_error(true_phases,phases2))
% disp(mean(abs(true_phases-phases2)));
disp('algo3');
disp(avg_error(true_phases_algo3,phases3))
% disp(mean(abs(true_phases_algo3-phases3)));

%% test
%     p = inputParser;
% %     paramName = 'myParam';
% %     defaultVal = 1;
% %     errorMsg = 'Value must be positive, scalar, and numeric.'; 
% %     validationFcn = @(x) isnumeric(x) && isscalar(x) ...
% %         && (x > 0);
%     addParameter(p,'myParam',false,...
%         @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
%     addRequired(p,'p1',@(x) isnumeric(x))
%     parse(p,100);
%     s=p.Results;
% % [phase, amplitude] = phastimate(data, D, edge, ord, hilbertwindow, varargin)
%% function definition

function [EEG]=add_pseudo_TMS(EEG,ISI,jitter)
    event_num=1;
    latency=0.0;
    EEG.event=struct('latency',1,'type','boundary','chan',[]);
%     EEG = pop_editeventvals(EEG,'append',{1,[],[],[]},'changefield',{event_num,'latency',latency},'changefield',{event_num,'type','boundary'});
    event_num=event_num+1;
    latency=latency+ISI-jitter+2*jitter.*rand(1);
    while latency<EEG.xmax-5
        EEG = pop_editeventvals(EEG,'append',{event_num-1,[],[],[]},'changefield',{event_num,'latency',latency},'changefield',{event_num,'type','TMS'});
        event_num=event_num+1;
        latency=latency+ISI-jitter+2*jitter.*rand(1);
    end
    t=num2cell(round([EEG.event.latency]));
    [EEG.event.latency]=t{:};
    EEG = eeg_checkset(EEG, 'eventconsistency');
end