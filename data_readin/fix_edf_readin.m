%% readin china edf data
path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
mkdir(path_datainfo)
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

for sub=1:length(allsubs)
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
