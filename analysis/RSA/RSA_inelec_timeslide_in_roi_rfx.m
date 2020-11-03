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

nrand=1000;

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
load('D:\matlab_tools\jet_grey.mat')
    
    
% for every sub check for elecs in the roi (then run rsa/stats seperately
% for each roi)
for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    all_labels=datainfo.elec_info.bipolar.montage_withoutartichan.labelnew;
     [~,~,ind]=intersect(all_labels,datainfo.elec_info.bipolar.elec_mni.label,'stable');     
     all_elec_label=[datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{ind,1}];

      for r=1:numel(rois)
            sel_rois=getfield(roi,rois{r});
            sel_rois={sel_rois{:}};
            [~,~,label_ind]=intersect( all_elec_label,sel_rois,'stable');  
            sel_roi_ind=[];
            for l=1:numel(label_ind)
           sel_roi_ind=[sel_roi_ind,find(strcmp(all_elec_label,sel_rois{label_ind(l)}))];
            end
          all_roi_labels{sub,r}=all_labels(sel_roi_ind);
      end
end
clear sel_rois all_label all_elec_label sel_roi_in label_ind





for c=9%1:numel(contrasts)
        contrast=contrasts{c};
       folder_out=fullfile(path_out,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
        mkdir(folder_out)

for r=1:numel(rois)
    sel_roi=rois{r};
    
    % subs with electrodes in this roi
    sel_subs=allsubs(cellfun(@isempty, all_roi_labels(:,r))==0);
    rand_rsa=zeros(numel(sel_subs),2,nrand,cfg_rsa.n_bins,cfg_rsa.n_bins);
    cond_rsa=zeros(numel(sel_subs),2,cfg_rsa.n_bins,cfg_rsa.n_bins);

for sub=1:length(sel_subs)
    sel_sub=sel_subs{sub};
    sub_ind=find(strcmp(sel_sub,allsubs));
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat_sym')))
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
    cfg.channel     = all_roi_labels{sub_ind,r};
    data=ft_selectdata(cfg,data);
    
  data.trlinfo=trlinfo;
  all_labels=data.label;
  
  
  % rsa not for every chan, but for every roi
  
  % need to have roi rsa for every sub with electrode in workspace
  
  trlinfo=data.trlinfo;


  % define features
            freq=ft_freqanalysis(cfg_freq,data);
            clear data
    
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
    
    %%%%%% here: create time bin vectors
    % reshape in feature vec
    
    % tmp_vec: chan_trialXtimebin
    % average for each bin
    n_bins=cfg_rsa.n_bins;
    num_trl=size(freq.powspctrm,1);
    bin_mat=zeros(num_trl,numel(freq.label),numel(freq.freq),n_bins);
    for t=1:cfg_rsa.n_bins
        bin_mat(:,:,:,t)=nanmean(freq.powspctrm(:,:,:,cfg_rsa.ind_t1(t):cfg_rsa.ind_t2(t)),4);
    end
    tmp_vec=permute(bin_mat,[1,4,2,3]);% trial x bin x chan x freq
    tmp_vec=reshape(tmp_vec,[],numel(freq.label)*numel(freq.freq));
    clear bin_mat
    
    rsa.rsa_mat=zeros(num_trl,num_trl,n_bins,n_bins,'single');
        tmp=corr(tmp_vec(:,:)','Type','Spearman');
        % fisher z
        %rsa.rsa_mat=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
        tmp=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
        tmp=reshape(tmp,num_trl,n_bins,num_trl,n_bins);
        rsa.rsa_mat(:,:,:,:)=single(permute(tmp,[1,3,2,4]));
        clear tmp
    clear tmp_vec
    
    % rsa_mat is symmetric (timextime), set half of the matrix to Nan to
    % avoid spurious large clusters in cluster stat
    nan_mat=repmat(reshape(triu(ones(n_bins),1),1,1,n_bins,n_bins),num_trl,num_trl,1,1);
    rsa.rsa_mat(logical(nan_mat))=NaN;
    
    %rsa.dim='trial_trial_time_time';
    rsa.dim='trial_time_time';
    
    
    
%     rsa.cfg.corr='Spearman';
%     rsa.cfg.norm=norm;
%     rsa.cfg.pow=cfg_freq;
%     rsa.trlinfo=trlinfo;
%     rsa.label=freq.label;
%     rsa.t1=cfg_rsa.t1;
%     rsa.t2=cfg_rsa.t2;
%     rsa.time=(cfg_rsa.t1+cfg_rsa.t2)./2;
  
   sortind=cfg_rsa.sortind;
        contrast_vec=reshape(cfg_rsa.contrast_mat,[],1);
        % select only non nan trials
        sel_ind=~isnan(contrast_vec);
        contrast_vec=contrast_vec(sel_ind);
        rsa.contrast_vec=contrast_vec;
        
        % sort rsa to match contrast_mat
   %     rsa.trlinfo=rsa.trlinfo(sortind,:);
        rsa.rsa_mat=rsa.rsa_mat(sortind,:,:,:);
        rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:);
        rsa.rsa_mat=reshape(rsa.rsa_mat,[],n_bins,n_bins);
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
end
end 

% first collect all channel for all subjects in a roi

% get rsa in each subject (spatial+power feature) sliding temp
% generalization



