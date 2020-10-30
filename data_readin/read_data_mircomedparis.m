%% readin in micromed data from paris (p_sub08)
addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults

path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
mkdir(path_datainfo)
mkdir(path_out)

sel_sub='p_sub08';

load(strcat(path_trialinfo,sel_sub,'_trlinfo.mat'))
sel_folder=strcat(path_data,sel_sub,'\ieeg\');
cd (sel_folder)
filename=dir('*.TRC');
cfg=[];
cfg.dataset      = filename.name;
cfg.continuous = 'yes';
cfg.channel    = {'all','-sti1+','-sti2+'};
data=ft_preprocessing(cfg);



cfg=[];
cfg.dataset      = filename.name;
cfg.continuous = 'yes';
cfg.channel={'sti1+','sti2+'};
trigger=ft_preprocessing(cfg)


% mark beginning and end of data using figure mark
f=figure
plot(trigger.trial{1}(1,:))
%waitfor(f)
%%%% here for p_sub08
learn_end_sp=6709116;
learn_start_sp=5167953;

test_end_sp=7236124;
test_start_sp=6783301;


% match learn log and trig, in ms samplerate
learn_chan=trigger.trial{1}(1,learn_start_sp:learn_end_sp);
learn_chan_sp_org=find(diff(learn_chan)>100);
learn_chan_sp_org(find(diff(learn_chan_sp_org)<500)+1)=[];
learn_chan_sp=learn_chan_sp_org./(data.fsample/1000);
learn_log_sp=trlinfo(trlinfo(:,2)<3,11)./10;
% create two vectors
learn_chan_vec=zeros(round(learn_chan_sp(end))+4000,1);

learn_chan_vec(round(learn_chan_sp))=1;
learn_log_vec=zeros(round(learn_log_sp(end))+4000,1);
learn_log_vec(round(learn_log_sp))=1;

numel(find(learn_chan_vec))


% convolve log vec for some uncertainty
learn_log_vec=conv(learn_log_vec,[zeros(400,1);ones(400,1);zeros(400,1)]);
[cross_corr,lag_corr]=xcorr(learn_log_vec,learn_chan_vec);
ind=find(cross_corr==max(cross_corr));
vec_to_chan=lag_corr(ind(4));


% there are inconsistencies in the data, add sample points to chan_vec to fix
log_sp1=812123;
chan_sp1=810471;
missing_sp=log_sp1-chan_sp1;
learn_chan_vec=[learn_chan_vec(1:chan_sp1-1);zeros(missing_sp,1);learn_chan_vec(chan_sp1:end)];

log_sp2=851293;
chan_sp2=849635;
missing_sp=log_sp2-chan_sp2;
learn_chan_vec=[learn_chan_vec(1:chan_sp2-1);zeros(missing_sp,1);learn_chan_vec(chan_sp2:end)];

log_sp3=1083553;
chan_sp3=1080236;
missing_sp=log_sp3-chan_sp3;
learn_chan_vec=[learn_chan_vec(1:chan_sp3-1);zeros(missing_sp,1);learn_chan_vec(chan_sp3:end)];


log_sp4=1442656;
chan_sp4=1418785;
missing_sp=log_sp4-chan_sp4;
learn_chan_vec=[learn_chan_vec(1:chan_sp4-1);zeros(missing_sp,1);learn_chan_vec(chan_sp4:end)];

figure
plot(learn_log_vec(vec_to_chan:end))
hold on
plot(learn_chan_vec)
ylim([-1,3])

for i=1:numel(learn_log_sp)
text(learn_log_sp(i)-vec_to_chan,1.5,num2str(i))
end

aligned_learn_log_vec=learn_log_vec(vec_to_chan+10:end);
end_sp=numel(learn_chan_vec)-numel(aligned_learn_log_vec);
aligned_learn_log_vec=[aligned_learn_log_vec;zeros(end_sp,1)];

matched_trig_vec=aligned_learn_log_vec.*learn_chan_vec;

trigger_learn=learn_chan_sp_org(find(aligned_learn_log_vec(find(learn_chan_vec))));
trigger_learn=trigger_learn+learn_start_sp;

learn_trl=ones(144,1);
learn_trial(133:134)=0;


figure
plot(trigger.trial{1}(1,:))
hold on
scatter(trigger_learn,ones(size(trigger_learn)))

clear aligned_learn_log_vec chan_sp1 chan_sp2 chan_sp3 chan_sp4
clear cross_corr end_sp ind lag_corr log_sp1 log_sp2 log_sp3 log_sp4
clear matched_trig_vec missing_sp vec_to_chan
close all
%%%%%test

% match learn log and trig, in ms samplerate
test_chan=trigger.trial{1}(1,test_start_sp:test_end_sp);
test_chan_sp_org=find(diff(test_chan)>100);
test_chan_sp_org(find(diff(test_chan_sp_org)<500)+1)=[];
test_chan_sp=test_chan_sp_org./(data.fsample/1000);
test_log_sp=trlinfo(trlinfo(:,2)==3,11)./10;
% create two vectors
test_chan_vec=zeros(round(test_chan_sp(end))+4000,1);

test_chan_vec(round(test_chan_sp))=1;
test_log_vec=zeros(round(test_log_sp(end))+4000,1);
test_log_vec(round(test_log_sp))=1;


% convolve log vec for some uncertainty
test_log_vec=conv(test_log_vec,[zeros(100,1);ones(100,1);zeros(100,1)]);
[cross_corr,lag_corr]=xcorr(test_log_vec,test_chan_vec);
ind=find(cross_corr==max(cross_corr));
vec_to_chan=lag_corr(ind(21));


% there are inconsistencies in the data, add sample points to chan_vec to fix
log_sp1=418362;
chan_sp1=417549;
missing_sp=log_sp1-chan_sp1;
test_chan_vec=[test_chan_vec(1:chan_sp1-1);zeros(missing_sp,1);test_chan_vec(chan_sp1:end)];


figure
plot(test_log_vec(vec_to_chan:end))
hold on
plot(test_chan_vec)
ylim([-1,3])

for i=1:numel(test_log_sp)
text(test_log_sp(i)-vec_to_chan,1.5,num2str(i))
end

aligned_test_log_vec=test_log_vec(vec_to_chan+10:end);
end_sp=numel(test_chan_vec)-numel(aligned_test_log_vec);
aligned_test_log_vec=[aligned_test_log_vec;zeros(end_sp,1)];

matched_trig_vec=aligned_test_log_vec.*test_chan_vec;

trigger_test=test_chan_sp_org(find(aligned_test_log_vec(find(test_chan_vec))));
trigger_test=trigger_test+test_start_sp;

test_trl=ones(48,1);

figure
plot(trigger.trial{1}(1,:))
hold on
scatter(trigger_test,ones(size(trigger_test)))

clear aligned_test_log_vec chan_sp1 chan_sp2 chan_sp3 chan_sp4
clear cross_corr end_sp ind lag_corr log_sp1 log_sp2 log_sp3 log_sp4
clear matched_trig_vec missing_sp vec_to_chan
close all



% save data and write triggerinfo in datainfo file    
% check trialdef format in fieldtrip
% save data and write triggerinfo in datainfo file    
% check trialdef format in fieldtrip
datainfo.trigger.trigger_type=(trlinfo(:,1)+100)';
datainfo.trigger.trigger_type=datainfo.trigger.trigger_type([learn_trl;test_trl])
trigger_sp=[trigger_learn,trigger_test];
datainfo.trigger.trigger_sp=trigger_sp;
datainfo.trigger.originalfile_trigger.file=filename.name;
datainfo.trigger.originalfile_trigger.trigger_channel='sti1';

datainfo.trialinfo=trlinfo([learn_trl;test_trl],:);
% 
% 
% there is a lot of data around the paradigm in the dataset, save diskspace
% and read/write time by recutting data
task_sp=[trigger_sp(1)-data.fsample.*60,trigger_sp(end)+data.fsample.*60];
num_sp=numel(task_sp(1):task_sp(end));

data.sampleinfo=[1,num_sp];
data.trial{1,1}=data.trial{1,1}(:,task_sp(1):task_sp(2));
data.time{1,1}=0:(1/data.fsample):((num_sp-1)*(1/data.fsample));
data.cfg.data_recut_to_task_sp_in_orgfile=task_sp;


datainfo.trigger.originalfile_trigger=datainfo.trigger;
datainfo.trigger.trigger_sp=trigger_sp-task_sp(1);


save(strcat(path_datainfo,sel_sub,'_datainfo.mat'),'datainfo')
save(strcat(path_out,sel_sub,'_data.mat'),'data','-v7.3')
