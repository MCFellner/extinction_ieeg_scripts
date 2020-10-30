%% data checker
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
mkdir(path_artifactinfo)
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
    
for sub=27:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)  
    load(strcat(path_preproc,sel_sub,'_data.mat'))

    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);

    % add arti_vec as additional channel (for easy artifact removal)
    data.trial{1}(end+1,:)=datainfo.artifact_info.browsercheck_bip.artifact_vec;
    data.label(end+1)={'artifact'};

    trlinfo=datainfo.trialinfo;
% segment data in different trial parts
% item window: -1 to 4
pre_item=1;
post_item=4;
trl_item(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_item);
trl_item(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
trl_item(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_item.*data.fsample);
datainfo.trigger.trl_item=trl_item;

pre_us=1;
post_us=4;
us_onset=datainfo.trigger.trigger_sp'+round(((trlinfo(:,13)-trlinfo(:,11))./10000).*data.fsample);
trl_us(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_us);
trl_us(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_us);
trl_us(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_us.*data.fsample);
datainfo.trigger.trl_us=trl_us;
cfg=[];
cfg.trl=trl_item;
data_item=ft_redefinetrial(cfg,data);
cfg=[];
cfg.keeptrials='yes';
erp_item=ft_timelockanalysis(cfg,data_item)
artifree_ind=squeeze(sum(erp_item.trial(:,end,:),3))==0;
clear erp_item
cfg=[];
cfg.trials=find(artifree_ind);
data_item=ft_preprocessing(cfg,data_item)


cfg=[];
cfg.viewmode='vertical';
ft_databrowser(cfg,data_item)
close all
delete(gcf)
clear
end