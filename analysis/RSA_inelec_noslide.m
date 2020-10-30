addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% rsa for iEEG

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\rsa\';

mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
toi=[3 4]; 
win=0.1; % in sec, power estimated for every win
feature='powlogscale';
norm='z_crosstrials';

% downsample to smallest sr
sr=1000;


for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
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
    
    
    % define features
    switch feature
        case 'powlogscale'
        cfg=[];
        cfg.output='pow';
        %cfg.toi =pre_item:0.1:post_item;
        cfg.keeptrials  = 'yes';
        cfg.pad='nextpow2';
        cfg.foi     = logspace(log10(2),2, 40);
        cfg.toi=toi(1):win:toi(2);
        cfg.method='wavelet';
        cfg.width = 5;       
        freq=ft_freqanalysis(cfg,data);
        clear data
    end
    
    num_trials=size(trlinfo,1);
    num_freq=numel(freq.freq);
    num_time=numel(freq.time);
    num_chan=numel(freq.label);
    
    switch norm
        case 'z_crosstrials'
         m=nanmean(freq.powspctrm,1);
         s=nanstd(freq.powspctrm,1);
         freq.powspctrm=(freq.powspctrm-repmat(m,num_trials,1,1,1))./repmat(s,num_trials,1,1,1);
    clear m s
    end
 
    
    % reshape in feature vec    
    tmp_vec=permute(freq.powspctrm,[3,4,1,2]);
    tmp_vec=reshape(tmp_vec,[],num_trials,num_chan);
    tmp_vec=permute(tmp_vec,[3,2,1]);
    rsa.dim='chan_trial_trial';
    
    % correlate
    % loop over chan
    rsa_mat=zeros(num_chan,num_trials,num_trials);
    for e=1:num_chan
    rsa_mat(e,:,:)=corr(squeeze(tmp_vec(e,:,:))','Type','Spearman');
    end
    
    % fisher z
    rsa_mat=0.5.*log(((ones(size(rsa_mat))+rsa_mat)./(ones(size(rsa_mat))-rsa_mat)));
       
    rsa.rsa_mat=rsa_mat;
    rsa.cfg.corr='Spearman';
    rsa.cfg.norm=norm;
    rsa.cfg.pow=cfg;
    rsa.trlinfo=trlinfo;
    rsa.label=freq.label;
    sel_folder=fullfile(path_out,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
    mkdir(sel_folder)
    save(fullfile(sel_folder,strcat(sel_sub,'_rsa_')),'rsa')
    clear rsa rsa_mat rmp_vec
end
    
    
%% stats for elecwise rsa mat
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_rsa='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

toi=[3 4]; 
win=0.1; % in sec, power estimated for every win
feature='powlogscale';
norm='z_crosstrials';
sel_folder=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
nrand=1000;

contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
            'cs_specific','cs_specific_block1','cs_specific_block2',...
            'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};

for c=1:numel(contrasts)
contrast=contrasts{c};

for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
% load rsa matrix
load(fullfile(sel_folder,strcat(sel_sub,'_rsa_')))
load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))

contrast_mat=getfield(contrast_def,contrast);

% sort rsa to match contrast_mat
sortind=contrast_def.sortind_org2usedtrlinfo;
rsa.trlinfo=rsa.trlinfo(sortind,:);
rsa.rsa_mat=rsa.rsa_mat(:,sortind,:);
rsa.rsa_mat=rsa.rsa_mat(:,:,sortind);


num_trial=size(rsa.trlinfo,1);
num_chan=numel(rsa.label);

% sel_col=rsa.trlinfo(:,col);
% all_cat=unique(sel_col);
% contrast_mat=zeros(num_trial);
% for i=1:numel(all_cat)
% tmp=sel_col==all_cat(i);
% tmp_mat=tmp*tmp';
% contrast_mat=contrast_mat+tmp_mat;
% end
% clear tmp tmp_mat
% % remove upper triangle (no same trial, no doubeling)
% upper=triu(ones(size(contrast_mat)));
% contrast_mat(logical(upper))=NaN;
% 
% % % sort for check
% % [~,ind]=sort(sel_col);
% % reshape design and rsa mat to vector
contrast_vec=reshape(contrast_mat,[],1);


rsa_tmp=permute(rsa.rsa_mat,[2,3,1]);
rsa_tmp=reshape(rsa_tmp,[],num_chan);

% select only non nan trials
sel_ind=~isnan(contrast_vec);
contrast_vec=contrast_vec(sel_ind);
rsa_tmp=rsa_tmp(sel_ind,:);
clear sel_ind

num_ind=numel(contrast_vec);

% data stat
[~,~,~,dat_stat]=ttest2(rsa_tmp(contrast_vec==1,:),rsa_tmp(contrast_vec==0,:));
data_t=dat_stat.tstat;
% permutation stats
all_rand_t=zeros(nrand,num_chan);
for r=1:nrand
rand_vec=contrast_vec(randperm(num_ind));

[~,~,~,rand_stat]=ttest2(rsa_tmp(rand_vec==1,:),rsa_tmp(rand_vec==0,:));
all_rand_t(r,:)=rand_stat.tstat;
end
all_rand_t=sort(all_rand_t);

for c=1:num_chan
    rank(c)=nearest(all_rand_t(:,c),data_t(c));
end

stat.rank=rank;
stat.label=rsa.label;
stat.nrand=nrand;
stat.contrast_mat=contrast_mat;

clear rank all_rand_t  constrast_mat rand_vec

folder_out=fullfile(sel_folder,contrast);
mkdir(folder_out)
save(fullfile(folder_out,strcat(sel_sub,'_rsastat')),'stat')
clear stat 
end
end
%%
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_rsa='D:\Extinction\iEEG\analysis\rsa\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

toi=[3 4]; 
feature='powlogscale';
norm='z_crosstrials';
contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
            'cs_specific','cs_specific_block1','cs_specific_block2',...
            'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};
                
sel_folder=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));

roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};
roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};

roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};

roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};



alpha_def=0.05;
for con=1:numel(contrasts)
    contrast=contrasts{con}

folder_out=fullfile(sel_folder,contrast);


sig_ind=[];
all_label=[];
all_pos=[];
for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    load(fullfile(folder_out,strcat(sel_sub,'_rsastat')))
    
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
       
        sig_tmp=zeros(size(stat.rank));
    sig_tmp=sig_tmp+(stat.rank>((1-alpha_def)*stat.nrand));
    sig_tmp=sig_tmp-(stat.rank<((alpha_def)*stat.nrand));

    sig_def{sub}=sig_tmp;
    all_elec{sub}=stat.label;
        
% get positions        
[~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');     
all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,1);

sig_ind=[sig_ind;sig_tmp'];
all_label=[all_label;all_elec{sub}'];
all_pos=[all_pos;all_elec_pos{sub}];
end

[region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def,roi)
save(fullfile(folder_out,'results_table'),'region_count','subject_count','roi_count')

%
num_all=sum([region_count.absnumelecinregion{:}])
num_sig=sum([region_count.absnumsigelec_pos{:}])
rel_sig=sum([region_count.absnumsigelec_pos{:}])/sum([region_count.absnumelecinregion{:}])

clear num_all num_sig rel_sig region_count 
% plot electrodes and count electrodes per region/pat
% sort elecs for subfunction in sig and no sig group
sort_ind=sig_ind;
all_ind=unique(sort_ind);
for i=1:numel(all_ind)
 elec_pos{i}=all_pos(sig_ind==all_ind(i),:);
 elec_label{i}=all_label(sig_ind==all_ind(i));
end

cfg.elec_pos =elec_pos;
cfg.elec_label=elec_label;
cfg.col_def=[0,0,1;0.5 0.5 0.5;1,0,0];
cfg.trans_def=[1 0.2 1];
views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
cfg.views=views;
cfg.hemispheres={'left','right'};
cfg.hemisphere_surf={'D:\matlab_tools\fieldtrip-20200130\template\anatomy\surface_pial_left.mat',...
    'D:\matlab_tools\fieldtrip-20200130\template\anatomy\surface_pial_right.mat'};

cfg.path_fig=folder_out;
mcf_ieegelectrodeplotter(cfg)
clear cfg elec_pos elec_label all_ind sort_ind sig_def all_elec_pos all_elec_label
end
