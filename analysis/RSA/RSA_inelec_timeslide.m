addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% rsa for iEEG

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';

mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

    
 contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
        'cs_specific','cs_specific_block1','cs_specific_block2',...
        'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
toi=[2 4];
win_pow=0.05; % in sec, power estimated for every win
feature='powlogscale';
norm='z_crosstrials';

% downsample to smallest sr
sr=1000;
% time window definitions
% define sliding window
sr_pow=1/(win_pow);
win=0.200; % duration of item
slide=0.05;

step=win/(1/sr_pow)+1;

cfg_rsa.tois=toi(1):1/sr_pow:toi(2);
cfg_rsa.t1=toi(1):slide:(toi(2)-win);
cfg_rsa.t2=cfg_rsa.t1+win;
cfg_rsa.ind_t1=1:slide/(1/sr_pow):((numel(cfg_rsa.tois)-win/(1/sr_pow)));
cfg_rsa.ind_t2=cfg_rsa.ind_t1+win/(1/sr_pow);
cfg_rsa.n_bins=numel(cfg_rsa.t1);

    switch feature
        case 'powlogscale'
            cfg_freq=[];
            cfg_freq.output='pow';
            %cfg.toi =pre_item:0.1:post_item;
            cfg_freq.keeptrials  = 'yes';
            cfg_freq.pad='nextpow2';
            cfg_freq.foi     = logspace(log10(2),log10(200), 50);
            cfg_freq.toi=cfg_rsa.tois;
            cfg_freq.method='wavelet';
           cfg_freq.width = 5;
    end

for c=1:numel(contrasts)
        contrast=contrasts{c};
        folder_out=fullfile(path_out,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
        mkdir(folder_out)

for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
        sel_folder=fullfile(path_out,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);        
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
        cfg_rsa.contrast_mat=getfield(contrast_def,contrast);
        cfg_rsa.sortind=contrast_def.sortind_org2usedtrlinfo;
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);
    clear montage
    
    % cut trials
    trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
    trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
    cfg=[];
    cfg.trl=trl;
    data=ft_redefinetrial(cfg,data);
    clear trl
    
    % downsample to common sampling rate
    cfg=[];
    cfg.resamplefs      = sr;
    cfg.detrend='yes';
    data=ft_resampledata(cfg,data);
    
    % only select artifree trials in trialinfo
    trlinfo=datainfo.trialinfo;
    trlinfo=trlinfo(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree,:);
    
    cfg=[]
    cfg.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
    data=ft_selectdata(cfg,data);
    data.trlinfo=trlinfo;
  all_labels=data.label;
   parfor chan=1:numel(data.label)
        sel_chan=all_labels{chan};
        all_stat{chan}=mcf_chanwise_rsa_stat(sel_chan,data,cfg_freq,cfg_rsa,norm)
    end        
        save(fullfile(folder_out,strcat(sel_sub,'_rsastat')),'all_stat')
    clear all_stat
end
end 

