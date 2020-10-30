%% readin neuralynx macro data
addpath('D:\matlab_tools\fieldtrip-20190108')
ft_defaults
addpath('D:\matlab_tools\MatlabImportExport_v6.0.0')
%%
% path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
% path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
% path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
% path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% mkdir(path_datainfo)
% mkdir(path_out)
% 
% allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};
% %allsubs = {'p_sub01','p_sub02','p_sub03','p_sub06','p_sub07'};
% 
% % check p_sub04 (two sess) and p_sub05 extra (discontunuous) 
% for sub=7:length(allsubs)
% sel_sub=allsubs{sub};
% 
% load(strcat(path_trialinfo,sel_sub,'_trlinfo.mat'))
% 
% if strcmp(sel_sub,'p_sub04')
% dataset1=strcat(path_data,sel_sub,'\ieeg\sess1');
% dataset2=strcat(path_data,sel_sub,'\ieeg\sess2');
% 
% hdr1=ft_read_header(dataset1,'headerformat','neuralynx_ds')
% data_raw1=ft_read_data(dataset1);
% event1=ft_read_event(dataset1)
% 
% hdr2=ft_read_header(dataset2,'headerformat','neuralynx_ds')
% data_raw2=ft_read_data(dataset2);
% event2=ft_read_event(dataset2)
% 
% % there is no sample point missing between the two files
% %(double(hdr1.LastTimeStamp)-double(hdr2.FirstTimeStamp)): 244 equals
% %exactly the interval between two sample points
% 
% % rawdata can be concatenated:
% clear data_raw data
% data_raw=[data_raw1,data_raw2];
% clear data_raw1 data_raw2
% % check transition
% % figure
% % plot(data_raw(10,hdr1.nSamples-100:hdr1.nSamples+100))
% 
% % event files need to be merged:
% % event1 has all information and correct samplepoint info
% event=event1;
% hdr=hdr1;
% hdr.originalhdr.hdr1=hdr1;
% hdr.originalhdr.hdr2=hdr2;
% elseif strcmp(sel_sub,'p_sub05')
% % there is a "data is discontonuous warning"
% dataset=strcat(path_data,sel_sub,'\ieeg\');
% hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')
% event=ft_read_event(dataset)
% data_raw=ft_read_data(dataset);
% 
% data.cfg.missingchannels={'F1a_5','F1a_6'};
% else
% dataset=strcat(path_data,sel_sub,'\ieeg\');
% hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')
% data_raw=ft_read_data(dataset);
% event=ft_read_event(dataset)
% end
% 
% 
% 
% % combine this info to fieldtrip like data structure
% data.label=hdr.label;
% data.fsample=hdr.Fs;
% data.hdr=hdr;
% data.sampleinfo=[1,hdr.nSamples];
% data.trial{1,1}=data_raw;
% data.time{1,1}=0:(1/hdr.Fs):((hdr.nSamples-1)*(1/hdr.Fs));
% data.cfg.dataset=dataset;
% data.cfg.hdr=hdr;
% data.cfg.event=event;
% 
% % check triggers and save in datainfo
% % vidA 101, vidB 102, vidC103
% ind_vidonset=find([event(:).value]==101|[event(:).value]==102|[event(:).value]==103);
% trigger_value=[event(ind_vidonset).value];
% trigger_sp=[event(ind_vidonset).sample];
% 
% if strcmp(sel_sub,'p_sub05')
% ind_vidonset(1:8)=[];
% trigger_value(1:8)=[];
% trigger_sp(1:8)=[];
% end
% 
% % compare trigger with info from logfile
%     % check number of triggers
%     if numel(ind_vidonset)~=size(trlinfo,1)
%         error('check why number of trigger do not match number of trials in trlinfo')
%     end
%     
%     % check diff times between triggers
%     trigger_diff=diff(trigger_sp)'.*(1000/data.fsample); %convert ms
%     trlinfo_diff=diff(trlinfo(:,11))./10; % convert ms
%     f=figure
%     plot(1:numel(trigger_diff),trigger_diff,1:numel( trlinfo_diff),trlinfo_diff)
%  %   waitfor(f); % waits for you to close the figure
%     %check_logfilematch=input('Do trialinfo and trigger match? Type 1 for yes 0 for no');
% %     if check_logfilematch==0
% %     error('check why logfile and trialinfo are not matching')
% %     end   
% 
%     
% % save data and write triggerinfo in datainfo file    
% % check trialdef format in fieldtrip
% % save data and write triggerinfo in datainfo file    
% % check trialdef format in fieldtrip
% datainfo.trigger.trigger_type=trigger_value;
% datainfo.trigger.trigger_sp=trigger_sp;
% datainfo.trigger.trigger_timestamp=[event(ind_vidonset).timestamp];
% datainfo.trialinfo=trlinfo;
% 
% 
% % there is a lot of data around the paradigm in the dataset, save diskspace
% % and read/write time by recutting data
% task_sp=[trigger_sp(1)-data.fsample.*60,trigger_sp(end)+data.fsample.*60];
% num_sp=numel(task_sp(1):task_sp(end));
% task_timestamp=[((task_sp(1)-1).*hdr.TimeStampPerSample)+event(1).timestamp,((task_sp(end)-1).*hdr.TimeStampPerSample)+event(1).timestamp];
% 
% data.sampleinfo=[1,num_sp];
% data.trial{1,1}=data.trial{1,1}(:,task_sp(1):task_sp(2));
% data.time{1,1}=0:(1/hdr.Fs):((num_sp-1)*(1/hdr.Fs));
% data.cfg.data_recut_to_task_sp_in_orgfile=task_sp;
% data.cfg.data_recut_to_task_timestamp_in_orgfile=task_timestamp;
% 
% 
% datainfo.trigger.originalfile_trigger=datainfo.trigger;
% datainfo.trigger.trigger_sp=trigger_sp-task_sp(1);
% datainfo.trigger.trigger_timestamp=[event(ind_vidonset).timestamp]-task_timestamp(1);
% 
% 
% %save(strcat(path_datainfo,sel_sub,'_datainfo.mat'),'datainfo')
% save(strcat(path_out,sel_sub,'_data.mat'),'data','-v7.3')
% 
% 
% end
% 
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
%% readin china edf data
path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
mkdir(path_datainfo)
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

      
for sub=5%:length(allsubs)
sel_sub=allsubs{sub};

load(strcat(path_trialinfo,sel_sub,'_trlinfo.mat'))
sel_folder=strcat(path_data,sel_sub,'\ieeg\');
cd (sel_folder)
filename=dir('*.edf');

% define trials
cfg            = [];
cfg.dataset    = filename.name;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data_all           = ft_preprocessing(cfg);

% cfg            = [];
% cfg.viewmode   = 'vertical';
% ft_databrowser(cfg, data_all);

% here paradigm has been restarted? or something like that, cut first part
% of the data 
if strcmp(sel_sub, 'c_sub06')
    start_sp=995811;
    data_all.trial{1}=data_all.trial{1}(:,start_sp:end);
    data_all.time{1}=data_all.time{1}(start_sp:end)-data_all.time{1}(start_sp);
    data_all.sampleinfo(2)=data_all.sampleinfo(2)-start_sp;
end
% get rid of 'POL' in label names
for i=1:numel(data_all.label)
   data_all.label{i}=data_all.label{i}(5:end); 
   if isempty(data_all.label{i})
       data_all.label{i}='POL';
   end
end

if strcmp(sel_sub, 'c_sub03')|| strcmp(sel_sub, 'c_sub04') || strcmp(sel_sub, 'c_sub05') || strcmp(sel_sub, 'c_sub06')  || strcmp(sel_sub, 'c_sub08')...
        || strcmp(sel_sub, 'c_sub11') || strcmp(sel_sub, 'c_sub12')|| strcmp(sel_sub, 'c_sub14') || strcmp(sel_sub, 'c_sub15') || strcmp(sel_sub, 'c_sub18')
 for i=1:numel(data_all.label)
    refmatch=strfind(data_all.label{i},'-Ref');
   if    ~isempty(refmatch)
       data_all.label{i}=data_all.label{i}(1:refmatch-1);
   end
end
end

if strcmp(sel_sub,  'c_sub15')
data_all.label{108}='S9';
end

if strcmp(sel_sub, 'c_sub10')  
    data_all.label{84}='J10';
end

if strcmp(sel_sub, 'c_sub20')  
    clean_chan=data_all.label([1:19,21:37,42:45,48:145]);
else
    
file_chan=dir('*_chaninfo.mat');
load(file_chan.name)
end

% select only oked channels (& EKG)
cfg            = [];
if strcmp(sel_sub,'c_sub20')
cfg.channel    = [clean_chan;{'SEKG-1';'SEKG-2'};{'-DC09';'-DC10';'-DC11';'-DC12'}];
data           = ft_preprocessing(cfg,data_all);
data.label{end-1}='EKG1';
data.label{end}='EKG2';
elseif strcmp(sel_sub,'c_sub08')
cfg.channel    = [clean_chan;{'EKG'};{'-DC09';'-DC10';'-DC11';'-DC12'}];
data           = ft_preprocessing(cfg,data_all);
data.label{end}='EKG1';
else
cfg.channel    = [clean_chan;{'EKG1';'EKG2'};{'-DC09';'-DC10';'-DC11';'-DC12'}];
data           = ft_preprocessing(cfg,data_all);
end
% 

% trigger part commented out: already saved and checked
% % select trigger channels
% % trigger channels: POL DC09-POL DC12
% cfg            = [];
% cfg.channel    = {'DC09';'DC10';'DC11';'DC12'};
% trigger_chan=ft_preprocessing(cfg,data_all);
% 
%     % recalculate trigger numbers (binary coding)
%     % select onset of trigger as event start (sanity check this)
%         % Phase AB    
%         % int vid_trigA=11;
%         % int vid_trigB=13;
%         % 
%         % int us_trig=3;
%         % int nous_trig=5;
%         % phase C
%         % int vid_trigC=15;
%         % int nous_trig=5;
% %     
% trig_tresh=100000;
% trigger_chan.trial{1}=diff(trigger_chan.trial{1},1,2)>trig_tresh;
% trig_onset_all=(trigger_chan.trial{1}(1,:))+(trigger_chan.trial{1}(2,:).*2)+(trigger_chan.trial{1}(3,:).*4)++(trigger_chan.trial{1}(4,:).*8);
% if strcmp(sel_sub,'c_sub02')||strcmp(sel_sub,'c_sub10')
% trig_onset=trig_onset_all.*(trig_onset_all==11|trig_onset_all==13|trig_onset_all==15|trig_onset_all==9);
% elseif strcmp(sel_sub,'c_sub09')
% trig_onset=trig_onset_all.*(trig_onset_all>=1);
%     
% else
% trig_onset=trig_onset_all.*(trig_onset_all==11|trig_onset_all==13|trig_onset_all==15);
% end
% [~,ind,trig]=(find(trig_onset));
% 
% if strcmp(sel_sub,'c_sub09')
%     ind(1:40)=[];
%     trig(1:40)=[];
%     ind(1499:1575)=[];
%     trig(1499:1575)=[];
% end
% % somehow some triggers for cue onset are recoded to 11, remove them to
% % only leave video onset triggers
% if strcmp(sel_sub,'c_sub09')
% diff_tresh=4.5.*data.fsample; % 2.5 sec min diff
% else
% diff_tresh=2.5.*data.fsample; % 2.5 sec min diff
% end
% diff_trig=[ind(1)+diff_tresh,diff(ind)];
% sel_trig=(find(diff_trig>diff_tresh));
% trigger_value=trig(sel_trig);
% trigger_sp=ind(sel_trig);
% 
% % 
% figure
% plot(trig_onset_all)
% hold on 
% scatter(trigger_sp,ones(size(trigger_sp)))
% 
% % compare trigger with info from logfile
%     % check number of triggers
%     if numel(sel_trig)~=size(trlinfo,1)
%         error('check why number of trigger do not match number of trials in trlinfo')
%     end
%     % check diff times between triggers
%     trigger_diff=diff(trigger_sp)'.*(1000/data.fsample); %convert ms
%     trlinfo_diff=diff(trlinfo(:,11))./10; % convert ms
%     f=figure
%     plot(1:numel(trigger_diff),trigger_diff,1:numel( trlinfo_diff),trlinfo_diff)
%     waitfor(f); % waits for you to close the figure
%     check_logfilematch=input('Do trialinfo and trigger match? Type 1 for yes 0 for no');
%     if check_logfilematch==0
%     error('check why logfile and trialinfo are not matching')
%     end   
%     
%     
% % save data and write triggerinfo in datainfo file    
% % check trialdef format in fieldtrip
% datainfo.trigger.trigger_type=trigger_value;
% datainfo.trigger.trigger_sp=trigger_sp;
% datainfo.trialinfo=trlinfo;
% save(strcat(path_datainfo,sel_sub,'_datainfo.mat'),'datainfo')


save(strcat(path_out,sel_sub,'_data.mat'),'data','-v7.3')
end

%% check trigger onsets with ERPs
addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults

path_data='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08',...
           'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

% trial_definition
t_pre=1; % in sec
t_post=8; % in sec
      
for sub=8%:numel(allsubs)
sel_sub=allsubs{sub};

load(strcat(path_data,sel_sub,'_data.mat'))
load(strcat(path_info,sel_sub,'_datainfo.mat'))

% cut data in trials using trigger_sp in datainfo
    %trial definition "trl" is an Nx3 matrix, N is the number of trials.
    % The first column contains the sample-indices of the begin of each trial
    % relative to the begin of the raw data, the second column contains the
    % sample-indices of the end of each trial, and the third column contains
    % the offset of the trigger with respect to the trial. An offset of 0
    % means that the first sample of the trial corresponds to the trigger. A
    % positive offset indicates that the first sample is later than the trigger,
    % a negative offset indicates that the trial begins before the trigger.
        sp_pre=t_pre*data.fsample;
        sp_post=t_post*data.fsample;

        trl_def(:,1)=double(datainfo.trigger.trigger_sp-sp_pre);
        trl_def(:,2)=double(datainfo.trigger.trigger_sp+sp_post);
        trl_def(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(0-sp_pre);

        cfg=[];
        cfg.trl=trl_def;
        data=ft_redefinetrial(cfg,data);
                
        cfg=[];
        cfg.keeptrials='no';
        erp=ft_timelockanalysis(cfg,data);
        cfg=[];
        cfg.baseline=[-0.2 0];
        erp=ft_timelockbaseline(cfg,erp);
        
        cfg=[];
        cfg.layout = 'ordered'
        layout = ft_prepare_layout(cfg, data)
        
        cfg=[];
        cfg.layout = layout;
        ft_multiplotER(cfg, erp)
        check=input('trigger timing ok? ')
end