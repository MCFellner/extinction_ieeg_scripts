addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')




%% plot all stats results
% path_in='D:\Extinction\iEEG\analysis\pow\rfx\';
%
% % all_contrasts={'Aplusplus_halfs','Aplusminus_halfs', 'Aminusminus_halfs',...
% %     'Bplusplus_halfs','Bplusminus_halfs','Bminusminus_halfs',...
% % 'int_Aplusminus_Aminusminus_halfs','int_Bplusminus_Bminusminus_halfs',...
% % 'int_Aplusminus_Aplusplus_halfs','int_Bplusminus_Bplusplus_halfs',...
% % 'int_Aminusminus_Aplusplus_halfs','int_Bminusminus_Bplusplus_halfs'};
%
% %all_contrasts={'A_cs_2half_vs_A_nocs_2half'};
% %all_contrasts={'Aus_nous'};
% %all_contrasts={'Aplus_halfs','int_A_plus_minus_halfs'};
%
% all_contrasts={'Aplus_minus',...
% 'Bplusplus_Bplusminus',...
% 'Bplusminus_Bminusminus',...
% 'Bplusplus_Bminusminus',...
% 'Aplus_minus_2half',...
% 'Bplusplus_Bplusminus_2half',...
% 'Bplusminus_Bminusminus_2half',...
% 'Bplusplus_Bminusminus_2half',...
% 'Aplus_minus_1half',...
% 'Bplusplus_Bplusminus_1half',...
% 'Bplusminus_Bminusminus_1half',...
% 'Bplusplus_Bminusminus_1half'};
%
%
% for con=1:numel(all_contrasts)
% sel_con=all_contrasts{con};
%
%
% path_roi=fullfile(path_in,sel_con)
% path_fig=fullfile(path_roi,'fig');
% mkdir(path_fig)
% all_roi=dir(strcat(path_roi,'\'));
% all_roi={all_roi(:).name};
% all_roi=all_roi(cellfun(@numel,all_roi)>2)
% fig_ind=find(strcmp(all_roi,'fig'));
% all_roi(fig_ind)=[];
%
% load('D:\matlab_tools\jet_grey2.mat')
% for r=1:numel(all_roi)
%     sel_roi=all_roi{r};
%     sel_path=fullfile(path_roi,sel_roi)
%
%     all_stats=dir([fullfile(sel_path,'lf'),'\clusterstat*']);
%     all_stats={all_stats(:).name};
%     for s=1:numel(all_stats)
%         % load stat
%         load(fullfile(sel_path,'hf',all_stats{s}))
%         stats{1}=stat;
%         load(fullfile(sel_path,'lf',all_stats{s}))
%         stats{2}=stat;
%
%         % plot stat.stat, add contours for sig clusters (black) and white contours
%         % for non sig
%         % plot results
%         fig=figure
%         for i=1:numel(stats)
%                 mask_pos=zeros(size(squeeze(stats{i}.stat)));
%                 mask_pos_sig=zeros(size(squeeze(stats{i}.stat)));
%                 mask_pos_trend=zeros(size(squeeze(stats{i}.stat)));
%                 mask_neg=zeros(size(squeeze(stats{i}.stat)));
%                 mask_neg_sig=zeros(size(squeeze(stats{i}.stat)));
%                 mask_neg_trend=zeros(size(squeeze(stats{i}.stat)));
%             % get sig cluster
%             if isfield(stats{i},'posclusters')
%                 if ~isempty(stats{i}.posclusters)
%                     if any([stats{i}.posclusters(:).prob]<=0.1)
%                         sig_ind=find([stats{i}.posclusters(:).prob]<=0.05);
%                         if ~isempty(sig_ind)
%                         mask_pos_sig=squeeze(stats{i}.posclusterslabelmat<=sig_ind(end)&stats{i}.posclusterslabelmat>0);
%                         else
%                             sig_ind=0;
%                         end
%                         trend_ind=find([stats{i}.posclusters(:).prob]<=0.1);
%                                 mask_pos_trend=squeeze(stats{i}.posclusterslabelmat<=trend_ind(end) ...
%                                                      & stats{i}.posclusterslabelmat>sig_ind(end)...
%                                                      & stats{i}.posclusterslabelmat>0);
%                         mask_pos=squeeze(stats{i}.posclusterslabelmat>sig_ind(end));
%                     else
%                     end
%                 else
%                 end
%             else
%             end
%
%
%            if isfield(stats{i},'negclusters')
%                 if ~isempty(stats{i}.negclusters)
%                     if any([stats{i}.negclusters(:).prob]<=0.1)
%                         sig_ind=find([stats{i}.negclusters(:).prob]<=0.05);
%                         if ~isempty(sig_ind)
%                         mask_neg_sig=squeeze(stats{i}.negclusterslabelmat<=sig_ind(end)&stats{i}.negclusterslabelmat>0);
%                         else
%                             sig_ind=0;
%                         end
%                         trend_ind=find([stats{i}.negclusters(:).prob]<=0.1);
%                                 mask_neg_trend=squeeze(stats{i}.negclusterslabelmat<=trend_ind(end) ...
%                                                      & stats{i}.negclusterslabelmat>sig_ind(end)...
%                                                      & stats{i}.negclusterslabelmat>0);
%                         mask_neg=squeeze(stats{i}.negclusterslabelmat>sig_ind(end));
%                     else
%                     end
%                 else
%                 end
%             else
%             end
%
%
%             % contour sig
%             mask_alpha=mask_pos_sig+mask_neg_sig;
%             % contour no sig
%             nosig_alpha=mask_pos+mask_neg;
%             trend_alpha=mask_pos_trend+mask_neg_trend;
%             subplot(2,1,i)
%
%             H= imagesc(stats{i}.time,stats{i}.freq,squeeze(stats{i}.stat),[-5 5])
%             colormap(jet_grey2)
%             set(gca,'YDir','normal')
%             %set(H,'AlphaData',mask_alpha)
%             hold on
%             contour(stats{i}.time,stats{i}.freq,mask_alpha,1,'LineColor','k','LineWidth',1)
%             contour(stats{i}.time,stats{i}.freq,nosig_alpha,1,'LineColor','w','LineWidth',1)
%             contour(stats{i}.time,stats{i}.freq,trend_alpha,1,':','LineColor','k','LineWidth',1)
%
%             title([sel_roi, ':', sel_con])
%
%             clear mask_alpha nosig_alpha mask_neg_sig mask_neg mask_pos mask_pos_sig mask_pos_trend mask_neg_trend
%
%         end
%         savefig(fig,fullfile(path_fig,[sel_roi,'_',all_stats{s}(1:end-4),'.fig']),'compact')
%
%     end
% end
% end


%% plot pow time courses in sig clusters

% get vp with elec in roi
%
% roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
% roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};
% roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
% roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
% roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% roi.amy_r={'Right-Amygdala'};
% roi.amy_l={'Left-Amygdala'};
% roi.hip_l={'Left-Hippocampus'};
% roi.hip_r={'Right-Hippocampus'};

% roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
% roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
%roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
%roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};
roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};

rois=fieldnames(roi);


% contrasts of interest
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\pow\rfx\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};


% define all conditions to run

% all_contrasts.Aplus_minus.conditions={'Aplus','Aminus'};
% all_contrasts.Aplus_minus.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aplus_minus.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'}];
% all_contrasts.Aplus_minus.roi='amy_l';
% all_contrasts.Aplus_minus.freq='hf';
% all_contrasts.Aplus_minus.stat_windows=[2 3];
% all_contrasts.Aplus_minus.cluster='posclusters';
%
% all_contrasts.Aplus_minus.conditions={'Aplus','Aminus'};
% all_contrasts.Aplus_minus.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aplus_minus.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'}];
% all_contrasts.Aplus_minus.roi='amy_l';
% all_contrasts.Aplus_minus.freq='hf';
% all_contrasts.Aplus_minus.stat_windows=[4 5];
% all_contrasts.Aplus_minus.cluster='posclusters';

all_contrasts.Aplus_minus.conditions={'Aplus','Aminus'};
all_contrasts.Aplus_minus.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Aplus_minus.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'}];
all_contrasts.Aplus_minus.roi='ventraltempocci';
all_contrasts.Aplus_minus.freq='lf';
all_contrasts.Aplus_minus.stat_windows=[2 3];
all_contrasts.Aplus_minus.cluster='negclusters';
%
% all_contrasts.Bplusminus_Bminusminus.conditions={'Bplusminus','Bminusminus'};
% all_contrasts.Bplusminus_Bminusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bplusminus_Bminusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'}];
% all_contrasts.Bplusminus_Bminusminus.roi='hip_l';
% all_contrasts.Bplusminus_Bminusminus.freq='hf';
% all_contrasts.Bplusminus_Bminusminus.stat_windows=[4 5];
% all_contrasts.Bplusminus_Bminusminus.cluster='negclusters';
%
% all_contrasts.Bplusminus_Bminusminus.conditions={'Bplusminus','Bminusminus'};
% all_contrasts.Bplusminus_Bminusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bplusminus_Bminusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'}];
% all_contrasts.Bplusminus_Bminusminus.roi='hip_r';
% all_contrasts.Bplusminus_Bminusminus.freq='lf';
% all_contrasts.Bplusminus_Bminusminus.stat_windows=[4 5];
% all_contrasts.Bplusminus_Bminusminus.cluster='negclusters';

% all_contrasts.Bplusminus_Bminusminus.conditions={'Bplusminus','Bminusminus'};
% all_contrasts.Bplusminus_Bminusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bplusminus_Bminusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'}];
% all_contrasts.Bplusminus_Bminusminus.roi='ventraltempocci';
% all_contrasts.Bplusminus_Bminusminus.freq='hf';
% all_contrasts.Bplusminus_Bminusminus.stat_windows=[4 5];
% all_contrasts.Bplusminus_Bminusminus.cluster='negclusters';


all_cons=fieldnames(all_contrasts);

pre_item=-1;
post_item=5.5;
%stat_windows=[4 4.5];
% downsample to smallest sr
sr=1000;

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



load('D:\matlab_tools\jet_grey2.mat')

for con=1:numel(all_cons)
    % seperate in condition 1 & 2 & average over freqs
    sel_con=all_cons{con};
    sel_conditions=getfield(all_contrasts,sel_con);
    conditions=getfield(sel_conditions,'conditions');
    cond_def=getfield(sel_conditions,'cond_def');
    sel_freq=getfield(sel_conditions,'freq');
    sel_roi=getfield(sel_conditions,'roi');
    sel_window=getfield(sel_conditions,'stat_windows');
    conditions=getfield(sel_conditions,'conditions');
    sel_cluster=getfield(sel_conditions,'cluster');
    
    % get electrodes for each sub/roi
    
    % subs with electrodes in this roi
    sel_subs=allsubs(cellfun(@isempty, all_roi_labels(:,r))==0);
    load(fullfile(path_out,strcat(sel_con),sel_roi,sel_freq,strcat('clusterstat_',num2str(sel_window(1).*1000),'to', num2str(sel_window(2).*1000))),'stat')
    % select freq of interest based on clusterstats
    cluster_prob=getfield(stat,sel_cluster);
    sel_cluster_num=find([cluster_prob(:).prob]<0.05);
    clustermat=getfield(stat,[sel_cluster,'labelmat']);
    cluster_def=squeeze(clustermat>0&clustermat<=sel_cluster_num(end));
    sel_freq_band=stat.freq(sum(cluster_def,2)>0);
    cluster_def=cluster_def(sum(cluster_def,2)>0,:);
    for sub=1:length(sel_subs)
        sel_sub=sel_subs{sub};
        sub_ind=find(strcmp(sel_sub,allsubs));
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
        
        
        cfg=[];
        cfg.output='pow';
        cfg.toi =pre_item:0.05:post_item;
        cfg.keeptrials  = 'yes';
        cfg.pad='nextpow2';
        switch sel_freq
            case 'lf'
                cfg.foi     =sel_freq_band;
                cfg.method='wavelet';
                cfg.width = 7;
            case 'hf'
                cfg.method      = 'mtmconvol';
                cfg.taper       = 'dpss';
                cfg.foi         = sel_freq_band;
                cfg.t_ftimwin   = ones(1,length(cfg.foi))*0.3;
                cfg.tapsmofrq   = ones(1,length(cfg.foi))*10;
        end
        freq=ft_freqanalysis(cfg,data);
        % z trans
        cfg=[];
        cfg.time=[pre_item+0.5,post_item-0.5];
        freq=z_trans_TF_seltime(cfg, freq);
        
        % select freq for same time as stat
        
        t_stat1=nearest(freq.time,stat.time(1));
        t_stat2=nearest(freq.time,stat.time(end));
        % select time window and average across channels
        sel_pow=squeeze(nanmean(freq.powspctrm(:,:,:,t_stat1:t_stat2),2));
        cluster_mat=repmat(reshape(cluster_def,[1,size(cluster_def)]),size(sel_pow,1),1,1);
        freq_all_trials{sub}=sum(sum(sel_pow.*cluster_mat,2),3)./sum(sum(cluster_mat,2),3);
        freq_all_trlinfo{sub}=trlinfo;
        %   freq_all_trials{c,sub}=squeeze(nanmean(nanmean(nanmean(freq.powspctrm(sel_trials,:,:,:),1),2),3));
        
        
        % define two conditions
        for c=1:numel(cond_def)
            sel_cond_def=cond_def{c};
            for d=1:size(sel_cond_def,1)
                eval(strcat('tmp(:,d)=trlinfo(:,',sel_cond_def{d,1},')',sel_cond_def{d,2}));
            end
            sel_trials=find(sum(tmp,2)==size(sel_cond_def,1));
            clear tmp
            
            freq_avg_cond(c,sub,:)=squeeze(nanmean(nanmean(nanmean(freq.powspctrm(sel_trials,:,:,:),1),2),3));
        end
    end
    
    clear  data datainfo montage sel_elec sel_elec_tmp sel_ind sel_trials  trl trlinfo
    
    % get default colors for boundedline
    f=  figure
    fig=subplot(1,2,1)
    cmap_default=fig.ColorOrder;
    close all
    % plot power curve in trial
    fig=figure
    subplot(5,6,[1,9])
    hold on
    
    for c=1:numel(conditions)
        x1=freq.time;
        y1=squeeze(nanmean(freq_avg_cond(c,:,:),2));
        b1=squeeze(nanstd(freq_avg_cond(c,:,:),[],2))./sqrt(size(freq_avg_cond(c,:,:),2));
        boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
    end
    title([sel_con,sel_roi,sel_freq,'clusterstat_',num2str(sel_window(1).*1000),'to', num2str(sel_window(2).*1000)])
    clear x1 y1 b1
    
    subplot(5,6,[26])
    hold on
    
    plot(1:11)
    plot(2:12)
    legend(conditions)
    subplot(5,6,[4,12])
    
    % contour sig
    mask_alpha=squeeze(stat.mask);
    imagesc(stat.time,stat.freq,squeeze(stat.stat),[-5 5])
    colormap(jet_grey2)
    set(gca,'YDir','normal')
    hold on
    contour(stat.time,stat.freq,mask_alpha,1,'LineColor','k','LineWidth',1)
    
    
    % plot: trial curves
    trl_pow=nan(size(freq_all_trials,2),3,64);
    
    for sub=1:size(freq_all_trials,2)
        freq= freq_all_trials{sub};
        trlinfo=freq_all_trlinfo{sub};
        for c=1:3
            sel_trl=  trlinfo(:,6)==c;
            sel_trlinfo=trlinfo(sel_trl,:);
            
            % fix indices of trial reps
            sel_ind=[sel_trlinfo(sel_trlinfo(:,2)==1,15);...
                sel_trlinfo(sel_trlinfo(:,2)==2,15)+24;...
                sel_trlinfo(sel_trlinfo(:,2)==3,15)+48];
            sel_pow=freq(sel_trl);
            trl_pow(sub,c,sel_ind)=sel_pow;
        end
    end
    trl_pow=smoothdata(trl_pow,3,'movmean',5);
    subplot(5,6,[13,21])
    hold on
    
    
    % add time clustering permstat
    % cluster t test
    
    % dimord trl_pow 'subj_freq_rpt_cond_timewindow'
    data1.dimord='subj_chan_time';
    data1.label={sel_roi};
    data1.time=1:64;
    data2=data1;
    data1.avg(:,1,:)=squeeze(nanmean(trl_pow(:,1:2,:),2));
    data2.avg(:,1,:)=squeeze(trl_pow(:,3,:));
    
    cfg=[];
    cfg.latency     = [data1.time(1),data1.time(end)];
    cfg.avgoverchan =  'yes';
    cfg.avgovertime =  'no';
    % first level
    cfg.method           = 'montecarlo';
    cfg.numrandomization = 1000;
    cfg.correctm         =  'cluster';
    cfg.correcttail      = 'no';
    cfg.alpha            =0.05;
    cfg.statistic ='depsamplesT';
    Nsub = size(data1.avg,1);                                       %# of subjects?
    design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub)];
    design(2,:)  = [1:Nsub 1:Nsub];
    cfg.clusteralpha     = 0.05;
    cfg.uvar     = 2;
    cfg.ivar     = 1;
    cfg.design = design;
    cfg.tail=1
    stat= ft_timelockstatistics(cfg,data1, data2);
    clear design data1 data2
    
    if any(stat.mask)
        scatter(find(stat.mask),zeros(size(find(stat.mask))),'*')
    end
    for c=1:3
        % tmp= smoothdata(avg_pow(f,:,c,t),'movmean',3);
        tmp=squeeze(nanmean( trl_pow(:,c,:)));
        x1=1:64;
        y1=tmp;
        b1=squeeze(nanstd(trl_pow(:,c,:))./sqrt(size( trl_pow,1)));
        boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
        
        % plot(tmp)
    end
    plot([24,24],[0,0.5],'k:')
    plot([48,48],[0,0.5],'k:')
    
    subplot(5,6,[30])
    hold on
    plot(1:11)
    plot(2:12)
    plot(3:13)
    legend({'plusplus','plusminus','minusminus'})
    
    savefig(fig,fullfile(path_out,strcat(sel_con),sel_roi,sel_freq,strcat('clusterstat_',num2str(sel_window(1).*1000),'to', num2str(sel_window(2).*1000),'power_timecourses')))
    close all
    clear trl_pow x1 y1 b1  freq_all_trialsfreq_all_trlinfo freq_avg_cond
end



