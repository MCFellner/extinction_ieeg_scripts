addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\extinction_ieeg_scripts\additional_functions')
addpath('D:\matlab_tools\CircStat')
%% rsa for iEEG

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};





% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1.5;
post_item=5;
toi=[0 3];
win_pow=0.001; % in sec, power estimated for every win
win=3; % duration of item
slide=3;

pow_feature='phase';
norm='no_norm';

% downsample to smallest sr
sr=1000;
nrand=1000;

%
all_roi.hip_l={'Left-Hippocampus'};
all_roi.hip_r={'Right-Hippocampus'};
all_roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
all_roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% %roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
all_roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% %roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
all_roi.amy_r={'Right-Amygdala'};
all_roi.amy_l={'Left-Amygdala'};
all_roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
all_roi.amy={'Left-Amygdala','Right-Amygdala'};
all_roi.hip={'Right-Hippocampus','Left-Hippocampus'};

rois=fieldnames(all_roi);

for r=1:numel(rois)
    roi=rois{r};
    sel_rois=getfield(all_roi,roi);
    
    
    %%%%%%%%%%%%%%
    load('D:\matlab_tools\jet_grey.mat')
    
    folder_out=fullfile(path_out,[pow_feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)],roi,'rsa_mat');
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
    rsa_ga.rand_rsa=zeros(numel(sel_subs),2,nrand,n_bins,n_bins);
    rsa_ga.cond_rsa=zeros(numel(sel_subs),2,n_bins,n_bins);
    
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
        clear cfg_preproc trl
        
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
        
        save(fullfile(folder_out,strcat(sel_sub,'_rsa')),'rsa')
        clear rsa
    end
end

%%
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';

toi=[0 3];
pow_feature='phase';
norm='no_norm';
n_freq=50;
nrand=1000;

contrasts={'video_specific_mask_block1','video_specific_mask_block2','video_specific_mask_block3',...
    'video_specific_mask_block1_interaction_video_specific_mask_block2',...
    'video_specific_mask_block1_interaction_video_specific_mask_block3',...
    'video_specific_mask_block2_interaction_video_specific_mask_block3'};
%

all_roi.hip_l={'Left-Hippocampus'};
all_roi.hip_r={'Right-Hippocampus'};
all_roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
all_roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% %roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
all_roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% %roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
all_roi.amy_r={'Right-Amygdala'};
all_roi.amy_l={'Left-Amygdala'};
all_roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
all_roi.amy={'Left-Amygdala','Right-Amygdala'};
all_roi.hip={'Right-Hippocampus','Left-Hippocampus'};

rois=fieldnames(all_roi);
%%%%%%%%%%%%%%
load('D:\matlab_tools\jet_grey.mat')
for r=1:numel(rois)
    roi=rois{r};
    
    
    folder_out=fullfile(path_out,[pow_feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)],roi);
    folder_rsa=fullfile(folder_out,'rsa_mat');
    all_subs=dir([folder_rsa,'\*.mat']);
    all_subs={all_subs(:).name};
    for cons=1:numel(contrasts)
        contrast=contrasts{cons};
        rsa_ga=[];
        for sub=1:numel(all_subs)
            sel_sub=all_subs{sub}(1:end-8);
            load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat_sym')))
            load(fullfile(folder_rsa,[sel_sub,'_rsa.mat']))
            
            %contrast_mat=getfield(contrast_def,contrast);
            contrast_mat=mcf_contrastmatdef(contrast_def,contrast);
            
            % get data for contrast & rand distribution
            % get condition and rand contrasts
            cfg_con.generate_rand='yes';
            cfg_con.nrand=nrand;
            cfg_con.contrast_mat=contrast_mat;
            cfg_con.sortind=contrast_def.sortind_org2usedtrlinfo;
            [rsa_cond]=mcf_rsacontrasts(cfg_con,rsa);
            
            rsa_ga.cond_rsa(sub,:,:)=rsa_cond.cond_rsa;
            rsa_ga.rand_rsa(sub,:,:,:)=rsa_cond.rand_rsa;
        end
        rsa_ga.freq=rsa_cond.freq;
        rsa_ga.dim_cond='subj_cond_freq';
        rsa_ga.dim_rand='subj_cond_rand_freq';
        rsa_ga.roi=rsa_cond.roi;
        
        % run data stats (using ft_freqstats)
        cfg_stats.nrand=nrand;
        cfg_stats.permutation='yes';
        cfg_stats.twosidedtest='yes'
        
        stats=mcf_rsacondstats(cfg_stats,rsa_ga)
        
        addpath(genpath('D:\matlab_tools\boundedline\kakearney-boundedline-pkg-8179f9a'))
        %         % plot result
        fig=figure
        colorscheme=[1 0 0; 0 0 1];
        hold on
        for cond=1:2
            x1=1:numel(rsa_ga.freq);
            y1=squeeze(nanmean(rsa_ga.cond_rsa(:,cond,:)));
            b1=squeeze(nanstd(rsa_ga.cond_rsa(:,cond,:)))./sqrt(numel(all_subs));
            
            boundedline(x1, y1, b1, 'cmap',colorscheme(cond,:),'transparency',0.1,'alpha');
        end
        sel_xticks=0:5:numel(rsa_ga.freq);
        sel_xticklabels=round([0,rsa_ga.freq(sel_xticks(2:end))]);
        xticks(sel_xticks)
        xticklabels(sel_xticklabels)
        % plot(rsa_ga.freq,stats.stat,'Color',colorscheme(cond,:))
        ylabel('similarity')
        
        xlabel('freq')
        title_str=strrep([roi,': ',contrast],'_',' ');
        title({title_str;;['pos tsum:',num2str(stats.trial_rand.data_pos(1)),'p=',num2str(stats.trial_rand.p_pos(1))];['neg tsum:',num2str(stats.trial_rand.data_neg(1)),'p=',num2str(stats.trial_rand.p_neg(1))]})
        
        ax_def=gca;
        hold on
        if any(stats.trial_rand.mask)
            scatter(find(stats.trial_rand.mask),ones(size(find(stats.trial_rand.mask))).*ax_def.YLim(2),'*')
        end
        legend([{'cond1'},{'stde'},{'cond2'},{'stde'}],'Location','northeastoutside')
        
        
        path_fig=fullfile( folder_out,'fig');
        mkdir(path_fig)
        savefig(fig,[path_fig,'\',contrast,'_in_',roi],'compact')
        
        save([path_fig,'\',contrast,'_in_',roi,'.mat'],'stats')
        close all
        clear rsa_cond rsa_ga stats
    end
end



