%% readin neuralynx macro data

% this script is a test to read in data with discontuous data 

addpath('D:\matlab_tools\fieldtrip-20190108')
ft_defaults
addpath('D:\matlab_tools\MatlabImportExport_v6.0.0')
addpath ('D:\Extinction\iEEG\scripts\additional_functions');

path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
mkdir(path_datainfo)
mkdir(path_out)

allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};
%allsubs = {'p_sub01','p_sub02','p_sub03','p_sub06','p_sub07'};

% check p_sub04 (two sess) and p_sub05 extra (discontunuous) 
for sub=5%:length(allsubs)
sel_sub=allsubs{sub};


load(strcat(path_trialinfo,sel_sub,'_trlinfo.mat'))

dataset=strcat(path_data,sel_sub,'\ieeg');
hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')
event=ft_read_event(dataset)
%data_raw=ft_read_data(dataset);

data=ft_read_neuralynx_interp(hdr.filename);

tmp_time=data.time;
data=rmfield(data,'time');
data.time{1}=tmp_time;

data.fsample=hdr.Fs;
ind_vidonset=find([event(:).value]==101|[event(:).value]==102|[event(:).value]==103);
trigger_value=[event(ind_vidonset).value];
trigger_sp=[event(ind_vidonset).sample];

ind_vidonset(1:8)=[];
trigger_value(1:8)=[];
trigger_sp(1:8)=[];

% compare trigger with info from logfile
    % check number of triggers
    if numel(ind_vidonset)~=size(trlinfo,1)
        error('check why number of trigger do not match number of trials in trlinfo')
    end
    
    % check diff times between triggers
    trigger_diff=diff(trigger_sp)'.*(1000/data.fsample); %convert ms
    trlinfo_diff=diff(trlinfo(:,11))./10; % convert ms
    f=figure
    plot(1:numel(trigger_diff),trigger_diff,1:numel( trlinfo_diff),trlinfo_diff)
    waitfor(f); % waits for you to close the figure
    %check_logfilematch=input('Do trialinfo and trigger match? Type 1 for yes 0 for no');
%     if check_logfilematch==0
%     error('check why logfile and trialinfo are not matching')
%     end   

    
% save data and write triggerinfo in datainfo file    
% check trialdef format in fieldtrip
% save data and write triggerinfo in datainfo file    
% check trialdef format in fieldtrip
datainfo.trigger.trigger_type=trigger_value;
datainfo.trigger.trigger_sp=trigger_sp;
datainfo.trigger.trigger_timestamp=[event(ind_vidonset).timestamp];
datainfo.trialinfo=trlinfo;


% there is a lot of data around the paradigm in the dataset, save diskspace
% and read/write time by recutting data
task_sp=[trigger_sp(1)-data.fsample.*60,trigger_sp(end)+data.fsample.*60];
num_sp=numel(task_sp(1):task_sp(end));
task_timestamp=[((task_sp(1)-1).*hdr.TimeStampPerSample)+event(1).timestamp,((task_sp(end)-1).*hdr.TimeStampPerSample)+event(1).timestamp];


data.sampleinfo=[1,num_sp];
data.trial{1,1}=data.trial{1,1}(:,task_sp(1):task_sp(2));
data.time=0:(1/hdr.Fs):((num_sp-1)*(1/hdr.Fs));
data.cfg.data_recut_to_task_sp_in_orgfile=task_sp;
data.cfg.data_recut_to_task_timestamp_in_orgfile=task_timestamp;

data=rmfield(data,'hdr')
data.cfg=cfg_tmp;
cfg            = [];
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);




end