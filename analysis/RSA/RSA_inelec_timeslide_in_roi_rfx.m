addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\extinction_ieeg_scripts\additional_functions')

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
win=0.200; % duration of item
slide=0.05;

pow_feature='powlogscale';
norm='z_crosstrials';

% downsample to smallest sr
sr=1000;



nrand=1000;


roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};
roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
%roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
%roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};

rois=fieldnames(roi);
sel_rois=getfield(roi,rois{1});
sel_rois={sel_rois{:}};

contrast=contrasts{1};

%%%%%%%%%%%%%%
load('D:\matlab_tools\jet_grey.mat')

folder_out=fullfile(path_out,[pow_feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)],'stats',contrast);
mkdir(folder_out)

% for every sub check for elecs in the roi (then run rsa/stats)
for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    cfg=[];
    cfg.atlas='freesurferDK';
    cfg.ana_labels='nearestGMlabelfreesurferDK';
    cfg.region= sel_rois;
    cfg.distance_region= 0;
    cfg.ref=  'bipolar'
    cfg.clean='clean' % currently only for bipolar
    elec_selection=mcf_elec_region_selector(cfg,datainfo)
    all_roi_labels{sub}=[elec_selection.labels_per_roi{:}];
    clear cfg
end




% subs with electrodes in this roi
sel_subs=allsubs(cellfun(@isempty, all_roi_labels(:))==0);
%preallocate for speed
n_bins=numel(toi(1):slide:(toi(2)-win));
rand_rsa=zeros(numel(sel_subs),2,nrand,n_bins,n_bins);
cond_rsa=zeros(numel(sel_subs),2,n_bins,n_bins);

for sub=1:length(sel_subs)
    sel_sub=sel_subs{sub};
    sub_ind=find(strcmp(sel_sub,allsubs));
    load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat_sym')))
    contrast_mat=getfield(contrast_def,contrast);
    
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    load(strcat(path_preproc,sel_sub,'_data.mat'))
        
    trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
    trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
    
    cfg_preproc=[];
    cfg_preproc.montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    cfg_preproc.resamplefs      = sr;
    cfg_preproc.trl_def=trl;
    cfg_preproc.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
    cfg_preproc.channel     = all_roi_labels{sub_ind};
    cfg_preproc.trialinfo=datainfo.trialinfo;
    data=mcf_preproc(cfg_preproc, data);
    clear cfg_preproc
            
   % defintion for features
    cfg_rsa.freqdef=pow_feature;
    cfg_rsa.norm=norm;
    % define sliding window
    cfg_rsa.toi=toi;
    cfg_rsa.win_pow=win_pow; % in sec, power estimated for every win
    cfg_rsa.win=win; % duration of item
    cfg_rsa.slide=slide;       
    rsa=mcf_timeslidepowrsa(cfg_rsa,data)
    
    
    cfg_rsa.contrast_mat=contrast_mat;
    cfg_rsa.sortind=contrast_def.sortind_org2usedtrlinfo;
    sortind=cfg_rsa.sortind;
    contrast_vec=reshape(cfg_rsa.contrast_mat,[],1);
    
    % sort rsa to match contrast_mat
    %     rsa.trialinfo=rsa.trialinfo(sortind,:);
    rsa.rsa_mat=rsa.rsa_mat(sortind,:,:,:);
    rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:);
    rsa.rsa_mat=reshape(rsa.rsa_mat,[],n_bins,n_bins);
    
    % select only non nan trials
    sel_ind=~isnan(contrast_vec);
    contrast_vec=contrast_vec(sel_ind);
    rsa.contrast_vec=contrast_vec;
    rsa.rsa_mat=rsa.rsa_mat(sel_ind,:,:);
    %%%%%%%%%%%%
    
    % combine data to real avg and random avg
    cond_rsa(sub,1,:,:)=nanmean(rsa.rsa_mat(rsa.contrast_vec==1,:,:));
    cond_rsa(sub,2,:,:)=nanmean(rsa.rsa_mat(rsa.contrast_vec==0,:,:));
    
    
    [~,rand_ind]=sort(rand(numel(rsa.contrast_vec),nrand),1);
    rand_vec=repmat(rsa.contrast_vec,1,nrand);
    rand_vec=rand_vec(rand_ind);
    clear rand_ind
    for i=1:nrand
        rand_rsa(sub,1,i,:,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==1,:,:));
        rand_rsa(sub,2,i,:,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==0,:,:));
    end
    clear rand_vec
    %      save(fullfile(folder_out,strcat(sel_sub,'_rsastat')),'all_stat')
    %  clear all_stat
end
rsa.t1=cfg_rsa.t1;
rsa.t2=cfg_rsa.t2;
rsa.time=(cfg_rsa.t1+cfg_rsa.t2)./2;

% run data stats (using ft_freqstats)

data_dummy.label={sel_roi};
data_dummy.freq=rsa.time;
data_dummy.time=rsa.time;
data_dummy.dimord='subj_chan_freq_time';
data1=data_dummy;
data1.powspctrm=cond_rsa(:,1,:,:);

data2=data_dummy;
data2.powspctrm=cond_rsa(:,2,:,:);

% define freqstats
cfg=[];
cfg.avgoverchan =  'yes';
cfg.avgovertime =  'no';
cfg.avgoverfreq =  'no';

% first level
cfg.method           = 'montecarlo';
cfg.numrandomization = nrand;
cfg.correctm         =  'cluster';
cfg.correcttail      = 'prob';
cfg.statistic ='depsamplesT'
% for within-subjects (depsamplesT)
Nsub = size(data1.powspctrm,1);                                       %# of subjects?
design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

cfg.uvar     = 2;
cfg.ivar     = 1;
cfg.design = design;
stat_data=ft_freqstatistics (cfg,data1,data2);
clear design


% parfor?
parfor i=1:nrand
    % run permutation stats (using_ft_freqstats)
    Nsub = size(rand_rsa,1);
    data1=data_dummy;
    data1.powspctrm=reshape(squeeze(rand_rsa(:,1,i,:,:)),Nsub,1,n_bins,n_bins);
    
    data2=data_dummy;
    data2.powspctrm=reshape(squeeze(rand_rsa(:,2,i,:,:)),Nsub,1,n_bins,n_bins);
    
    cfg=[];
    cfg.avgoverchan =  'yes';
    cfg.avgovertime =  'no';
    cfg.avgoverfreq =  'no';
    
    % first level
    cfg.method           = 'montecarlo';
    cfg.numrandomization = 2;
    cfg.correctm         =  'cluster';
    cfg.correcttail      = 'prob';
    cfg.statistic ='depsamplesT'
    % for within-subjects (depsamplesT)
    Nsub = size(data1.powspctrm,1);
    %# of subjects?
    design=[];
    design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
    design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
    
    cfg.uvar     = 2;
    cfg.ivar     = 1;
    cfg.design = design;
    stat_rand=ft_freqstatistics (cfg,data1,data2);
    % get pos/neg clusters for each random
    rand_pos(i)=0;
    rand_neg(i)=0;
    if isfield(stat_rand,'posclusters')
        if ~isempty(stat_rand.posclusters)
            rand_pos(i)=stat_rand.posclusters(1).clusterstat;
        end
    end
    if isfield(stat_rand,'negclusters')
        if ~isempty(stat_rand.negclusters)
            rand_neg(i)=stat_rand.negclusters(1).clusterstat;
        end
    end
end


% check significance

rand_neg=sort(rand_neg,'ascend');
rand_pos=sort(rand_pos,'descend');

data_pos=0;
data_neg=0;
p_pos=1;
p_neg=1;
if isfield(stat_data,'posclusters')
    if ~isempty(stat_data.posclusters)
        data_pos=[stat_data.posclusters(:).clusterstat];
        for i=1:numel(data_pos)
            p_pos(i)=(nearest(rand_pos,data_pos(i))./nrand).*0.5;
        end
    end
end
if isfield(stat_data,'negclusters')
    if ~isempty(stat_data.negclusters)
        data_neg=[stat_data.negclusters(:).clusterstat];
        for i=1:numel(data_neg)
            p_neg(i)=(nearest(rand_neg,data_neg(i))./nrand).*0.5;
        end
    end
end


if isfield(stat_data,'posclusterslabelmat')
    maskpos=(stat_data.posclusterslabelmat<=sum(p_pos<0.05)&stat_data.posclusterslabelmat>0);
else
    maskpos=zeros(size(stat_data.stat));
end
if isfield(stat_data,'negclusterslabelmat')
    maskneg=(stat_data.negclusterslabelmat<=sum(p_neg<0.05)&stat_data.negclusterslabelmat>0);
else
    maskneg=zeros(size(stat_data.stat));
end

mask=maskpos+maskneg;

fig=figure
imagesc(stat_data.time,stat_data.time,squeeze(stat_data.stat),[-5 5])
hold on
colormap(jet_grey)
colorbar
title({[contrast,' in ',sel_roi];['pos tsum:',num2str(data_pos(1)),'p=',num2str(p_pos(1))];['neg tsum:',num2str(data_neg(1)),'p=',num2str(p_neg(1))]})
ylabel('t in s')
xlabel('t in s')
contour(stat_data.time,stat_data.time,squeeze(mask),1,'k')
set(gca,'YDir','normal')
% plot results
path_fig=fullfile( folder_out,'fig');
savefig(fig,[path_fig,'\',contrast,'_in_',sel_roi,],'compact')

stat_data.trial_rand.p_pos=p_pos;
stat_data.trial_rand.p_neg=p_neg;

stat_data.trial_rand.rand_pos=rand_pos;
stat_data.trial_rand.rand_neg=rand_neg;
save([path_fig,'\',contrast,'_in_',sel_roi,'.mat'],'stat_data')
close all
clear stat_data stat_rand rand_pos rand_neg p_pos p_neg data_neg data_pos


% first collect all channel for all subjects in a roi

% get rsa in each subject (spatial+power feature) sliding temp
% generalization



