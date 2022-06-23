% EEGLAB history file generated on the 02-May-2022
% ------------------------------------------------

EEG.etc.eeglabvers = '2022.0'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data','G:\\我的雲端硬碟\\Documents\\110上學期\\研究\\phase_mep\\SPIS-Resting-State-Dataset-master\\Pre-SART EEG\\S02_restingPre_EC.mat','srate',256,'pnts',0,'xmin',0);
EEG.setname='test';
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = pop_editeventvals(EEG,'append',{1,[],[],[]},'changefield',{2,'latency',0.5},'changefield',{2,'type','Test'});
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = pop_editeventvals(EEG,'append',{2,[],[],[]},'changefield',{3,'latency',1},'changefield',{3,'type','test'});
EEG = eeg_checkset( EEG );
