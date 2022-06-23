% EEGLAB history file generated on the 26-May-2022
% ------------------------------------------------
test=rand(1,256*2,4);
eeglab;
EEG.etc.eeglabvers = '2022.0'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_importdata('dataformat','array','nbchan',1,'data','test','srate',256,'pnts',256*2,'xmin',-2);
EEG.setname='test';
EEG.epoch=struct('event',num2cell(1:4),'eventlatency',0,'eventtype','TMS','eventchan',[]);
EEG.event=struct('latency',{256*2*1,256*2*2,256*2*3,256*2*4},'type','boundary','chan',[],'epoch',num2cell(1:4));
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);

