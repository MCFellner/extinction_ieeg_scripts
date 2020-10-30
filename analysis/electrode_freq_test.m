addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')


%% freqs: check for significant freq electrodes: difference during item/
% difference anticipating us
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\pow\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
stat_windows=[2 2.5; 2.5 3; 3 3.5; 3.5 4]; % add rows for more windows
fois=[2 4;5 7;8 12;13 20;20 30; 35  60; 70 100];
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};

% downsample to smallest sr
sr=1000;

conditions={'B_switch','B_plus'};
cond1_def=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
cond2_def=[{'2'},{'==2'};{'6'},{'==1'}];
% conditions={'A_cs_2half','A_nocs_2half'};
% cond1_def=[{'2'},{'==1'};{'8'},{'==1'};{'15'},{'>=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% cond2_def=[{'2'},{'==1'};{'8'},{'==0'};{'15'},{'>=12'}];

for sub=8%18:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);
    
    
    % cut trials
    trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
    trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
    cfg=[];
    cfg.trl=trl;
    data=ft_redefinetrial(cfg,data);
    
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
    
    
    % define two conditions
    for d=1:size(cond1_def,1)
        eval(strcat('tmp(:,d)=trlinfo(:,',cond1_def{d,1},')',cond1_def{d,2}));
    end
    sel_trials1=find(sum(tmp,2)==size(cond1_def,1));
    clear tmp
    
    for d=1:size(cond2_def,1)
        eval(strcat('tmp(:,d)=trlinfo(:,',cond2_def{d,1},')',cond2_def{d,2}));
    end
    sel_trials2=find(sum(tmp,2)==size(cond2_def,1));
    clear tmp
    
    
    for f=1:length(fois)
        foi=fois(f,:);
        f_name=foi_names{f};
                
        cfg=[];
        cfg.output='pow';
        cfg.toi =pre_item:0.05:post_item;
        cfg.keeptrials  = 'yes';
        cfg.pad='nextpow2';
        switch f_name
            case {'theta1','theta2','alpha','beta1','beta2'}
                cfg.foilim     = foi;
                cfg.method='wavelet';
                cfg.width = 5;
            case {'gamma1','gamma2'}
                cfg.method      = 'mtmconvol';
                cfg.taper       = 'dpss';
                cfg.foi         = foi(1):5:foi(2);
                cfg.t_ftimwin   = ones(1,length(cfg.foi))*0.3;
                cfg.tapsmofrq   = ones(1,length(cfg.foi))*10;
        end
        freq=ft_freqanalysis(cfg,data);
        
        % z trans
        cfg=[];
        cfg.time=[pre_item+0.5,post_item-0.5];
        freq=z_trans_TF_seltime(cfg, freq);
        
        % seperate in condition 1 & 2 & average over freqs
        freq1=freq;
        freq1.powspctrm=nanmean(freq.powspctrm(sel_trials1,:,:,:),3);
        freq1.cumtapcnt=nanmean(freq.cumtapcnt(sel_trials1,:),2);
        freq1.freq=nanmean(freq.freq);
        freq2=freq;
        freq2.powspctrm=nanmean(freq.powspctrm(sel_trials2,:,:,:),3);
        freq2.cumtapcnt=nanmean(freq.cumtapcnt(sel_trials2,:),2);
        freq2.freq=nanmean(freq.freq);
        
        % loop over channels
        for t=1:size(stat_windows,1)
            sel_window=   stat_windows(t,:);
            parfor e=1:numel(freq1.label)
                sel_elec=freq1.label{e};
                               
                % contrast
                cfg=[];
                cfg.latency     = sel_window;
                cfg.channel   = freq1.label{e};
                cfg.avgoverchan =  'yes';
                cfg.avgovertime =  'no';
                cfg.avgoverfreq =  'yes';

                % first level
                cfg.method           = 'montecarlo';
                cfg.numrandomization = 10000;
                cfg.correctm         =  'cluster';
                cfg.correcttail      = 'prob';
                cfg.statistic ='indepsamplesT'
                design = zeros(1,size(freq1.powspctrm,1) + size(freq2.powspctrm,1));
                design(1,1:size(freq1.powspctrm,1)) = 1;
                design(1,(size(freq1.powspctrm,1)+1):(size(freq1.powspctrm,1) + size(freq2.powspctrm,1)))= 2;
                cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
                cfg.design = design;
                stat=ft_freqstatistics (cfg,freq1,freq2)
                
                if isfield(stat,'posclusters')
                    if ~isempty(stat.posclusters)
                        pos_prob(e)=stat.posclusters(1).prob;
                    else
                        pos_prob(e)=1;
                    end
                else
                    pos_prob(e)=1;
                end
                
                if isfield(stat, 'negclusters')
                    if ~isempty(stat.negclusters)
                        neg_prob(e)=stat.negclusters(1).prob;
                    else
                        neg_prob(e)=1;
                    end
                else
                    neg_prob(e)=1;
                end
            end
            freqstat.label=freq1.label;
            freqstat.cfg=cfg;
            freqstat.prob_pos=pos_prob;
            freqstat.prob_neg=neg_prob;
            freqstat.time_window= sel_window;
            freqstat.conditions=conditions;
            freqstat.cond1_def=cond1_def;
            freqstat.cond2_def=cond2_def;
            freqstat.sel_trials1=sel_trials1;
            freqstat.sel_trials2=sel_trials2;
            freqstat.foi=foi;
            freqstat.f_name=f_name;
            
            % save in specific folder
            folder_out=fullfile(path_out,strcat(conditions{1},'_vs_',conditions{2}),strcat(conditions{1},'_vs_',conditions{2},'_',f_name,num2str(foi(1)),'to',num2str(foi(2)),'Hz','_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));
            mkdir(folder_out)
            save(fullfile(folder_out,strcat(sel_sub,'freqstat')),'freqstat')
            clear freqstat pos_prob neg_prob stat trl
        end
    end
    clear sel_trials1 sel_trials2 freq freq1 freq2 trlinfo 
end

%%
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_pow='D:\Extinction\iEEG\analysis\pow\';


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};


stat_windows=[2 2.5; 2.5 3; 3 3.5; 3.5 4];
fois=[2 4;5 7;8 12;13 20;20 30; 35  60; 70 100];
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};


alpha_def=0.05;

%conditions={'A_us','A_nous'};
%conditions={'B_switch','B_plus'};
conditions={'A_cs_2half','A_nocs_2half'};

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


 for f=1:length(fois)
        foi=fois(f,:);
        f_name=foi_names{f};
  for t=1:size(stat_windows,1)
            sel_window=   stat_windows(t,:);       
        
  path_in=fullfile(path_pow,strcat(conditions{1},'_vs_',conditions{2}),strcat(conditions{1},'_vs_',conditions{2},'_',f_name,num2str(foi(1)),'to',num2str(foi(2)),'Hz','_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));


sig_ind=[];
all_label=[];
all_pos=[];
for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    load(fullfile(path_in,strcat(sel_sub,'freqstat')))
    
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
       
    stat=freqstat;
    sig_tmp=zeros(size(stat.prob_pos));
    sig_tmp=sig_tmp+(stat.prob_pos<alpha_def);
    sig_tmp=sig_tmp-(stat.prob_neg<alpha_def);

sig_def{sub}=sig_tmp;
all_elec{sub}=stat.label;
        
% get positions        
[~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');     
all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,1);

sig_ind=[sig_ind;sig_tmp'];
all_label=[all_label;stat.label'];
all_pos=[all_pos;all_elec_pos{sub}];
end

[region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def,roi)
save(fullfile(path_in,'results_table'),'region_count','subject_count','roi_count')

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

cfg.path_fig=path_in;
mcf_ieegelectrodeplotter(cfg)

end
 end
  
 
 %% combined results across frequencies
 
 
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_pow='D:\Extinction\iEEG\analysis\pow\';
 stat_windows=[2 2.5; 2.5 3; 3 3.5; 3.5 4];
%stat_windows=[4 4.5];
fois=[2 4;5 7;8 12;13 20;20 30; 35  60; 70 100];
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};

%conditions={'A_us','A_nous'};
conditions={'B_switch','B_plus'};
%conditions={'A_cs_2half','A_nocs_2half'};

load('D:\matlab_tools\jet_grey.mat')
colormap_pos=[ones(2,3);jet_grey(33:end,:)];
tmp_ind=32:-1:1;
colormap_neg=[ones(2,3);jet_grey(tmp_ind,:)];
figure 
for t=1:size(stat_windows,1)
            sel_window=   stat_windows(t,:);     
           
            
 for f=1:length(fois)
        foi=fois(f,:);
        f_name=foi_names{f};
           
path_in=fullfile(path_pow,strcat(conditions{1},'_vs_',conditions{2}),strcat(conditions{1},'_vs_',conditions{2},'_',f_name,num2str(foi(1)),'to',num2str(foi(2)),'Hz','_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));
load(fullfile(path_in,'results_table'))

% roi rel_sigelec_pos
tmp=roi_count.relsigelec_pos;
roi_pos(:,f)=[tmp{:}];

tmp=roi_count.relsigelec_neg;
roi_neg(:,f)=[tmp{:}];
% plot data using imgesc
% ,'region_count','subject_count','roi_count'
 end
 
 rois=roi_count.roilabel;
ax(1)=subplot(size(stat_windows,1),2,(t*2)-1)
imagesc(roi_pos,[0 1])
xticks(1:numel(fois))
xticklabels(foi_names)
yticks(1:numel(rois))
yticklabels(rois)
title(strcat('pos diff:',conditions{1},' vs ',conditions{2},' toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'))
colormap(ax(1),colormap_pos)
colorbar
ax(2)=subplot(size(stat_windows,1),2,t*2)
imagesc(roi_neg,[0 1])
xticks(1:numel(fois))
xticklabels(foi_names)
yticks(1:numel(rois))
yticklabels(rois)
title(strcat('neg diff:',conditions{1},'vs ',conditions{2},' toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'))
colormap(ax(2),colormap_neg)
colorbar

end