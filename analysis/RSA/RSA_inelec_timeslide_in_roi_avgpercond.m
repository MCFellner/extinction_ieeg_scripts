addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\extinction_ieeg_scripts\additional_functions')

%% average rsa values for contrasts

% get rsa value frm sig clusters
% or get values from timeXtime window

% define cluster

% extract averages for different conditions (same as the contrast mats?)

% different

% plot time course of contrasts using sliding avg definition
% plot bars using half definition
% plot bars using block definition


path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};



%%%%%%%%%%%% % define cluster/ time window to select data from
all_contrasts=[];

all_contrasts{1}.contrast= 'type1to2_vs_type2to3_mask_block1_interaction_type1to2_vs_type2to3_mask_block2';
all_contrasts{1}.roi='amy';
all_contrasts{1}.toi=[2 3];

all_contrasts{2}.contrast= 'type1to2_vs_type2to3_mask_block1_interaction_type1to2_vs_type2to3_mask_block2';
all_contrasts{2}.roi='amy';
all_contrasts{2}.toi=[2 3.8];

all_contrasts{3}.contrast= 'item_specific_mask_block1_interaction_item_specific_mask_block3';
all_contrasts{3}.roi='ifg';
all_contrasts{3}.toi=[2 3];

all_contrasts{4}.contrast= 'item_specific_mask_block1_interaction_item_specific_mask_block3';
all_contrasts{4}.roi='ifg';
all_contrasts{4}.toi=[2 3.8];

all_contrasts{5}.contrast= 'item_specific_mask_block1_interaction_item_specific_mask_block2';
all_contrasts{5}.roi='dm_pfc';
all_contrasts{5}.toi=[2 3.8];

all_contrasts{6}.contrast= 'type1to2_vs_type2to3_mask_block1_interaction_type1to2_vs_type2to3_mask_block2';
all_contrasts{6}.roi='dm_pfc';
all_contrasts{6}.toi=[2 3.8];

all_contrasts{7}.contrast= 'type1to2_vs_type2to3_mask_block1_interaction_type1to2_vs_type2to3_mask_block2';
all_contrasts{7}.roi='dm_pfc';
all_contrasts{7}.toi=[2 3];

all_contrasts{8}.contrast= 'type1to2_vs_type2to3_mask_block1_interaction_type1to2_vs_type2to3_mask_block2';
all_contrasts{8}.roi='hip_r';
all_contrasts{8}.toi=[2 3.8];

all_contrasts{9}.contrast= 'item_specific';
all_contrasts{9}.roi='ventraltempocci';
all_contrasts{9}.toi=[2 3];

all_contrasts{10}.contrast= 'item_specific';
all_contrasts{10}.roi='ventraltempocci';
all_contrasts{10}.toi=[2 3.8];
  %%
for con_clus=1:numel(all_contrasts)
    contrast= all_contrasts{con_clus}.contrast;
    roi= all_contrasts{con_clus}.roi;
    toi=all_contrasts{con_clus}.toi;
    
  
    
    %%%%%%%%%rsa definitions
    % segment data in different trial parts
    % item window: -1 to 4 (imp
    pre_item=-1.5;
    post_item=5;
    win_pow=0.05; % in sec, power estimated for every win
    win=0.200; % duration of item
    slide=0.05;
    
    pow_feature='powlogscale';
    norm='z_crosstrials';
    
    % downsample to smallest sr
    sr=1000;
    
    
    folder_stats=fullfile(path_out,[pow_feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)],'stats',contrast);
    
    load(fullfile(folder_stats,'\fig\',[contrast,'_in_',roi,'.mat']));
    cluster_def.type='stat';
    cluster_def.stat_mask=stats.trial_rand.mask;
    cluster_def.stat_tx=stats.time; % define time
    cluster_def.stat_ty=stats.freq; % define time to ensure correct selection
    % timewindow
    % cluster_def.type='window';
    % cluster_def.tx=[2 3]; % defined in sec start/end
    % cluster_def.ty=[2 3];
    
    %%%%%%%% contrasts to plots (single plot for each group)
    % every contrast gets trial course, barplots for blocks (subblocks)
    
    % item specfic contrast
    % valence contrast
    
    contrasts={'item_specific','type1to2_vs_type2to3'};
    mask_def={'block','halfs','trialcourse'};
    contrast_masks={    {'block1','block2','block3'},...
        {'first_half_block1','second_half_block1',...
        'first_half_block2','second_half_block2',...
        'first_half_block3','second_half_block3'},...
        {'trial_slidingavg_def'}};
    %%%%%%%%%%%%%%%
    
    % all rois definitions
    all_roi.hip_l={'Left-Hippocampus'};
    all_roi.hip_r={'Right-Hippocampus'};
    all_roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
    all_roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
    all_roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
    all_roi.amy_r={'Right-Amygdala'};
    all_roi.amy_l={'Left-Amygdala'};
    all_roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
    all_roi.amy={'Left-Amygdala','Right-Amygdala'};
    all_roi.hip={'Right-Hippocampus','Left-Hippocampus'};
    
    
    rois=fieldnames(all_roi);
    sel_rois=getfield(all_roi,roi);
    
    
    %%%%%%%%%%%%%%
    load('D:\matlab_tools\jet_grey.mat')
    
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
    
    n_bins=numel(toi(1):slide:(toi(2)-win));
    
    % initialize result matrix
    
    all_contrast_rsa=[];
    for con_m=1:numel(contrasts)
        sel_con=contrasts{con_m};
        for m=1:numel(mask_def)
            sel_def=mask_def{m};
            sel_masks=contrast_masks{m};
            all_contrast_rsa=setfield(all_contrast_rsa,[sel_con,'_x_',sel_def],[]);
        end
    end
    
    for sub=1:length(sel_subs)
        sel_sub=sel_subs{sub};
        sub_ind=find(strcmp(sel_sub,allsubs));
        
        
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
        clear trl
        
        % defintion for features
        cfg_rsa.freqdef=pow_feature;
        cfg_rsa.norm=norm;
        % define sliding window
        cfg_rsa.toi=toi;
        cfg_rsa.win_pow=win_pow; % in sec, power estimated for every win
        cfg_rsa.win=win; % duration of item
        cfg_rsa.slide=slide;
        cfg_rsa.roi=roi;
        rsa=mcf_timeslidepowrsa(cfg_rsa,data)
        
        % select data from defined window/cluster
        % average data for choosen masks
        cfg_sel=cluster_def;
        rsa=mcf_rsadataselect(cfg_sel,rsa);
        
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat_sym')))
        
        % sort trials to same order as in contrat_mats
        sortind=contrast_def.sortind_org2usedtrlinfo;
        rsa.rsa_mat=rsa.rsa_mat(sortind,:);
        rsa.rsa_mat=rsa.rsa_mat(:,sortind);
        %contrast_mat=getfield(contrast_def,contrast);
        
        % vectorize rsa_mat
        rsa_vec=reshape(rsa.rsa_mat,numel(rsa.rsa_mat),1);
        
        all_contrast_mat=[];
        % selcontrast is contrast x mask
        for con_m=1:numel(contrasts)
            sel_con=contrasts{con_m};
            
            for m=1:numel(mask_def)
                sel_def=mask_def{m};
                sel_masks=contrast_masks{m};
                for s=1:numel(sel_masks)
                    sel_mask=sel_masks{s};
                    
                    sel_contrast =[sel_con,'_mask_',sel_mask];
                    tmp=mcf_contrastmatdef(contrast_def,sel_contrast);
                    contrast_mat_tmp{s}=reshape(tmp,[],size(rsa_vec,1));
                    clear tmp
                end
                contrast_mat_tmp=vertcat(contrast_mat_tmp{:});
                all_contrast_mat=setfield(all_contrast_mat,[sel_con,'_x_',sel_def],contrast_mat_tmp);
                clear contrast_mat_tmp sel_contrast
            end
        end
        clear con con_m tmp sel_con sel_def sel_contrast sel_mask m
        %     %%%%%%%%% get average values for each contrast_vec x mask_vec
        %
        all_con=fieldnames(all_contrast_mat);
        for con=1:numel(all_con)
            sel_con=all_con{con};
            sel_mat= getfield(all_contrast_mat,sel_con);
            tmp=getfield(all_contrast_rsa,sel_con);
            for c=1:size(sel_mat,1)
                tmp(sub,c,1)= mean( rsa_vec(sel_mat(c,:)==1));
                tmp(sub,c,2)= mean( rsa_vec(sel_mat(c,:)==2));
            end
            all_contrast_rsa=setfield(all_contrast_rsa,sel_con,tmp);
            clear tmp
        end
        % plot result
        
        
        %
        %     fig=figure
        %     imagesc(stats.time,stats.time,squeeze(stats.stat),[-5 5])
        %     hold on
        %     colormap(jet_grey)
        %     colorbar
        %     title({[contrast];[roi];['pos tsum:',num2str(stats.trial_rand.data_pos(1)),'p=',num2str(stats.trial_rand.p_pos(1))];['neg tsum:',num2str(stats.trial_rand.data_neg(1)),'p=',num2str(stats.trial_rand.p_neg(1))]})
        %     ylabel('t in s')
        %     xlabel('t in s')
        %     contour(stats.time,stats.time,squeeze(stats.trial_rand.mask),1,'k')
        %     set(gca,'YDir','normal')
        %     %         path_fig=fullfile( folder_out,'fig');
        %     %         mkdir(path_fig)
        %     %         savefig(fig,[path_fig,'\',contrast,'_in_',roi],'compact')
        %     %
        %     %         save([path_fig,'\',contrast,'_in_',roi,'.mat'],'stats')
        %     close all
    end
    
    all_con=fieldnames(all_contrast_rsa);
    % item specific figure
    
    addpath('D:\matlab_tools\Violinplot-Matlab-master')
    addpath(genpath('D:\matlab_tools\boundedline\kakearney-boundedline-pkg-8179f9a'))
    
    folder_fig=fullfile(folder_stats,roi,'cluster_extracted_values');
    mkdir(folder_fig)
    %% figure itemspecific
    fig=figure
    
    rsa_mat=reshape(permute(all_contrast_rsa.item_specific_x_block,[1,3,2]),numel(sel_subs),[]);
    factors={'type','block'};
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',2,1),[],1);
    type=reshape(repmat([{'wi'},{'bi'}],3,1)',[],1);
    cat_mat=[type,block];
    colorscheme=[1 0 0; 0 0 1];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('itemspecific x block')
    
    savefig(fig,fullfile(folder_fig,'itemspecific_effects.fig'))
    close all
    
    %% figure itemspecific half
    fig=figure
    
    rsa_mat=reshape(permute(all_contrast_rsa.item_specific_x_halfs,[1,3,2]),numel(sel_subs),[]);
    factors={'type','block','half'};
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',4,1,1),[],1);
    type=reshape(repmat([{'wi'},{'bi'}]',6,1,1),[],1);
    half=reshape(repmat([{'first'},{'second'}],2,3,1),[],1);
    cat_mat=[type,block,half];
    colorscheme=[1 0 0; 0 0 1];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('itemspecific x block x halfs')
    
    savefig(fig,fullfile(folder_fig,'itemspecifichalfs_effects.fig'))
    close all
    
    
    %% trialcourse item specific
    colorscheme=[1 0 0; 0 0 1];
    
    co = colorscheme;
    set(groot,'defaultAxesColorOrder',co)
    cmap_default=co;
    fig=figure
    hold on
    subplot(2,2,1:2)
    for ty=1:2
        x1=1:size(all_contrast_rsa.item_specific_x_trialcourse,2)
        y1=squeeze(nanmean(all_contrast_rsa.item_specific_x_trialcourse(:,:,ty)));
        b1=squeeze(nanstd(all_contrast_rsa.item_specific_x_trialcourse(:,:,ty)))./sqrt(numel(sel_subs));
        
        boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
    end
    ylabel('similarity')
    xlabel('trialnumber')
    title('trial course  itemspecific')
    
    % add block markers
    ax_def=gca;
    hold on
    plot([24 24],ax_def.YLim,'k:')
    plot([48 48],ax_def.YLim,'k:')
    
    % add cluster permutation
    data1_mat=permute(all_contrast_rsa.item_specific_x_trialcourse(:,:,1),[1,3,2]);
    data2_mat=permute(all_contrast_rsa.item_specific_x_trialcourse(:,:,2),[1,3,2]);
    data1.dimord='subj_chan_time';
    data1.label={contrast};
    data1.time=1:size(data2_mat,3);
    data2=data1;
    data1.avg=data1_mat;
    data2.avg=data2_mat;
    
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
    cfg.correcttail='prob';
    stat= ft_timelockstatistics(cfg,data1, data2);
    
    % add stars at significant clusters
    if any(stat.mask)
        hold on
        scatter(find(stat.mask),ones(size(find(stat.mask))).*ax_def.YLim(2),'*')
    end
    legend([{'wi'},{'stde'},{'bi'},{'stde'}],'Location','northeastoutside')
    
    savefig(fig,fullfile(folder_fig,'itemspecifictrialcourse_effects.fig'))
    close all
    set(groot,'defaultAxesColorOrder','remove')
    
    %% figure item specific wi-bi difference
    fig=figure
    
    rsa_mat=all_contrast_rsa.item_specific_x_block(:,:,1)-all_contrast_rsa.item_specific_x_block(:,:,2);
    factors={'block'};
    block=[{'block1'};{'block2'};{'block3'}];
    cat_mat=[block];
    colorscheme=[0.5 0 0.5];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('itemspecific x block')
    
    savefig(fig,fullfile(folder_fig,'itemspecific_block_effects.fig'))
    close all
    
    %% figure itemspecific half differences
    fig=figure
    
    rsa_mat=all_contrast_rsa.item_specific_x_halfs(:,:,1)-all_contrast_rsa.item_specific_x_halfs(:,:,2);
    factors={'block','half'};
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',2,1),[],1);
    half=reshape(repmat([{'first'},{'second'}],1,3,1),[],1);
    cat_mat=[block,half];
    colorscheme=[1 0 0; 0 0 1];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('itemspecific x halfs')
    
    savefig(fig,fullfile(folder_fig,'itemspecific_halfs_effects.fig'))
    close all
    
    
    %% figure valence
    fig=figure
    
    rsa_mat=reshape(permute(all_contrast_rsa.type1to2_vs_type2to3_x_block,[1,3,2]),numel(sel_subs),[]);
    factors={'type','block'};
    
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',2,1),[],1);
    type=reshape(repmat([{'1to2'},{'2to3'}],3,1)',[],1);
    cat_mat=[type,block];
    colorscheme=[0 0.5 0.5; 0.5 0.5 1];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('valence specific x block')
    
    savefig(fig,fullfile(folder_fig,'valencespecific_effects.fig'))
    close all
    
    %% figure valence x half
    fig=figure
    
    rsa_mat=reshape(permute(all_contrast_rsa.type1to2_vs_type2to3_x_halfs,[1,3,2]),numel(sel_subs),[]);
    factors={'type','block','half'};
    
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',4,1),[],1);
    type=reshape(repmat([{'1to2'},{'2to3'}],6,1)',[],1);
    half=reshape(repmat([{'first'},{'second'}],2,3,1),[],1);
    
    cat_mat=[type,block,half];
    colorscheme=[0 0.5 0.5; 0.5 0.5 1];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('valence specific x half x block')
    
    savefig(fig,fullfile(folder_fig,'valencespecifichalf_effects.fig'))
    close all
    
    %% trialcourse valence specific
    colorscheme=[0 0.5 0.5; 0.5 0.5 1];
    
    co = colorscheme;
    set(groot,'defaultAxesColorOrder',co)
    cmap_default=co;
    fig=figure
    hold on
    subplot(2,2,1:2)
    for ty=1:2
        x1=1:size(all_contrast_rsa.type1to2_vs_type2to3_x_trialcourse,2)
        y1=squeeze(nanmean(all_contrast_rsa.type1to2_vs_type2to3_x_trialcourse(:,:,ty)));
        b1=squeeze(nanstd(all_contrast_rsa.type1to2_vs_type2to3_x_trialcourse(:,:,ty)))./sqrt(numel(sel_subs));
        
        boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
    end
    ylabel('similarity')
    xlabel('trialnumber')
    title('valence: trial course')
    % add block markers
    ax_def=gca;
        hold on

    plot([24 24],ax_def.YLim,'k:')
    plot([48 48],ax_def.YLim,'k:')
        hold on

    % add cluster permutation
    data1_mat=permute(all_contrast_rsa.type1to2_vs_type2to3_x_trialcourse(:,:,1),[1,3,2]);
    data2_mat=permute(all_contrast_rsa.type1to2_vs_type2to3_x_trialcourse(:,:,2),[1,3,2]);
    data1.dimord='subj_chan_time';
    data1.label={contrast};
    data1.time=1:size(data2_mat,3);
    data2=data1;
    data1.avg=data1_mat;
    data2.avg=data2_mat;
    
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
    cfg.correcttail='prob';
    stat= ft_timelockstatistics(cfg,data1, data2);
    
    % add stars at significant clusters
    hold on
    if any(stat.mask)
            hold on

        scatter(find(stat.mask),ones(size(find(stat.mask))).*ax_def.YLim(2),'*')
    end
    legend([{'type1to2'},{'stde'},{'type2to3'},{'stde'}],'Location','northeastoutside')
    
    savefig(fig,fullfile(folder_fig,'valencespecifictrialcourse_effects.fig'))
    close all
    set(groot,'defaultAxesColorOrder','remove')
    %% figure valence diff
    
    fig=figure
    
    rsa_mat=all_contrast_rsa.type1to2_vs_type2to3_x_block(:,:,1)-all_contrast_rsa.type1to2_vs_type2to3_x_block(:,:,2);
    factors={'block'};
    
    block=[{'block1'};{'block2'};{'block3'}];
    cat_mat=[block];
    colorscheme=[1 0 0.5];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('valence specific x block')
    
    savefig(fig,fullfile(folder_fig,'valencespecific_block_effects.fig'))
    close all
    
    %% figure valence diff halfs
    
    fig=figure
    
    rsa_mat=all_contrast_rsa.type1to2_vs_type2to3_x_halfs(:,:,1)-all_contrast_rsa.type1to2_vs_type2to3_x_halfs(:,:,2);
    factors={'block','half'};
    block=reshape(repmat([{'block1'};{'block2'};{'block3'}]',2,1,1),[],1);
    half=reshape(repmat([{'first'},{'second'}],1,3,1),[],1);cat_mat=[block];
    colorscheme=[1 0 0.5];
    cat_mat=[block,half];
    
    rmtable=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme);
    sgtitle('valence specific: halfs x block')
    
    savefig(fig,fullfile(folder_fig,'valencespecific_halfs_effects.fig'))
    close all
    keep path_info path_preproc path_out path_designmat allsubs  all_contrasts
end