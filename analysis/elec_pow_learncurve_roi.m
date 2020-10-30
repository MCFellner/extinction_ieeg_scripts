addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')


%% roi specific rfx regression analysis (combining sub specific t-maps)
% get vp with elec in roi
%
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

rois=fieldnames(roi);


% contrasts of interest
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_in='D:\Extinction\iEEG\analysis\pow\regression\';
path_out='D:\Extinction\iEEG\analysis\pow\regression\roi\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};


stat_windows=[0 4.5;2 3;3 4]; % add rows for more windows

% define all conditions to run
all_con{1}.conditions={'type1_A'};
all_con{1}.type='vs0';

all_con{2}.conditions={'type2_A'};
all_con{2}.type='vs0';

all_con{3}.conditions={'type3_A'};
all_con{3}.type='vs0';

all_con{4}.conditions={'type1_B'};
all_con{4}.type='vs0';

all_con{5}.conditions={'type2_B'};
all_con{5}.type='vs0';

all_con{6}.conditions={'type3_B'};
all_con{6}.type='vs0';


all_con{7}.conditions={'type1_A','type3_A'};
all_con{7}.type='simple_contrast';

all_con{8}.conditions={'type1_A','type2_A'};
all_con{8}.type='simple_contrast';

all_con{9}.conditions={'type2_A','type3_A'};
all_con{9}.type='simple_contrast';

all_con{10}.conditions={'type1_B','type3_B'};
all_con{10}.type='simple_contrast';

all_con{11}.conditions={'type1_B','type2_B'};
all_con{11}.type='simple_contrast';

all_con{12}.conditions={'type2_B','type3_B'};
all_con{12}.type='simple_contrast';

for con=1:numel(all_con)
    conditions=all_con{con}.conditions;
    type=all_con{con}.type;

    freqs={'hf'}%,'hf'}
    for f=1:numel(freqs)
        freq_def=freqs{f};

        % get electrodes for each sub/roi
        for r=1:numel(rois)
            sel_roi=rois{r};
            roi_def=getfield(roi,sel_roi);
            sub_num=0;
            for sub=1:length(allsubs)
                sel_sub=allsubs{sub};
                % electrodeinfo
                info_file=strcat(path_info,sel_sub,'_datainfo');
                load(info_file)
                sel_elec_tmp=datainfo.elec_info.bipolar.elec_ct_mr.label(ismember([datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{:}],[roi_def]));
                % only select after preproc elecs
                sel_elec=intersect(sel_elec_tmp,datainfo.artifact_info.rejectvisual_bip.elecsin);

                count_elec(sub,r)=numel(sel_elec);

                if numel(sel_elec)>0
                    sub_num=sub_num+1;


                    for c=1:numel(conditions)
                        sel_cond=conditions{c};
                        sel_path=fullfile(path_in,sel_cond);
                        load(fullfile(sel_path,[sel_sub,'_',freq_def,'_reg_tmap']))

                        [~,elec_ind]=intersect(freq_reg.label,sel_elec,'stable');
                        freq_all{c,sub_num}=freq_reg;
                        freq_all{c,sub_num}.powspctrm=nanmean(freq_reg.powspctrm(elec_ind,:,:),1);
                        freq_all{c,sub_num}.dimord='chan_freq_time';
                        freq_all{c,sub_num}.label={sel_roi};



                    end
                else
                end
                clear freq data datainfo montage sel_elec sel_elec_tmp sel_ind sel_trials  trl trlinfo

            end

            for c=1:numel(conditions)
                sel_cond=conditions{c};

                cfg=[];
                cfg.keepindividual='yes';
                ga{c}=ft_freqgrandaverage(cfg,freq_all{c,:});
            end

            switch type
                case 'interaction'
                    ga1=ga{1};
                    ga2=ga{1};
                    ga1.powspctrm=ga{1}.powspctrm-ga{2}.powspctrm;
                    ga2.powspctrm=ga{3}.powspctrm-ga{4}.powspctrm;
                    cond_name=[conditions{1},'-',conditions{2},'_vs_',conditions{3},'-',conditions{4}];

                case 'simple_contrast'
                    ga1=ga{1};
                    ga2=ga{2};
                    cond_name=[conditions{1},'_vs_',conditions{2}];

                case 'vs0'
                    ga1=ga{1};
                    ga2=ga1;
                    ga2.powspctrm=zeros(size(ga1.powspctrm));
                    cond_name=conditions{1};
            end
            clear ga

            for w=1:size(stat_windows)
                sel_window=stat_windows(w,:);

                cfg=[];
                cfg.latency     = sel_window;
                cfg.avgoverchan =  'yes';
                cfg.avgovertime =  'no';
                cfg.avgoverfreq =  'no';

                % first level
                cfg.method           = 'montecarlo';
                cfg.numrandomization = 1000;
                cfg.correctm         =  'cluster';
                cfg.correcttail      = 'prob';
                cfg.statistic ='depsamplesT'
                % for within-subjects (depsamplesT)
                Nsub = size(ga1.powspctrm,1);                                       %# of subjects?
                design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
                design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

                cfg.uvar     = 2;
                cfg.ivar     = 1;
                cfg.design = design;
                stat=ft_freqstatistics (cfg,ga1,ga2)

                sel_path=fullfile(path_out,cond_name,sel_roi,freq_def);
                mkdir(sel_path)
                save(fullfile(sel_path,strcat('clusterstat_',num2str(sel_window(1).*1000),'to', num2str(sel_window(2).*1000))),'stat')
                clear design

            end
            clear freq_all

        end
    end
end

%% plot all stats results
path_in='D:\Extinction\iEEG\analysis\pow\regression\roi\';

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

rois=fieldnames(roi);


all_con{1}.conditions={'type1_A'};
all_con{1}.type='vs0';

all_con{2}.conditions={'type2_A'};
all_con{2}.type='vs0';

all_con{3}.conditions={'type3_A'};
all_con{3}.type='vs0';

all_con{4}.conditions={'type1_B'};
all_con{4}.type='vs0';

all_con{5}.conditions={'type2_B'};
all_con{5}.type='vs0';

all_con{6}.conditions={'type3_B'};
all_con{6}.type='vs0';


all_con{7}.conditions={'type1_A','type3_A'};
all_con{7}.type='simple_contrast';

all_con{8}.conditions={'type1_A','type2_A'};
all_con{8}.type='simple_contrast';

all_con{9}.conditions={'type2_A','type3_A'};
all_con{9}.type='simple_contrast';

all_con{10}.conditions={'type1_B','type3_B'};
all_con{10}.type='simple_contrast';

all_con{11}.conditions={'type1_B','type2_B'};
all_con{11}.type='simple_contrast';

all_con{12}.conditions={'type2_B','type3_B'};
all_con{12}.type='simple_contrast';

for con=1:numel(all_con)
    conditions=all_con{con}.conditions;
    type=all_con{con}.type;
    
    switch type
        case 'interaction'
            cond_name=[conditions{1},'-',conditions{2},'_vs_',conditions{3},'-',conditions{4}];
        case 'simple_contrast'
            cond_name=[conditions{1},'_vs_',conditions{2}];
        case 'vs0'
            cond_name=conditions{1};
    end
    
    
    
    
    all_roi=rois;
    
    
    load('D:\matlab_tools\jet_grey2.mat')
    for r=1:numel(all_roi)
        sel_roi=all_roi{r};
        path_roi=fullfile(path_in,cond_name);
        path_fig=fullfile(path_roi,'fig');
        mkdir(path_fig)
        sel_path=fullfile(path_roi,sel_roi)
        
        all_stats=dir([fullfile(sel_path,'lf'),'\clusterstat*']);
        all_stats={all_stats(:).name};
        for s=1:numel(all_stats)
            % load stat
            load(fullfile(sel_path,'hf',all_stats{s}))
            stats{1}=stat;
            load(fullfile(sel_path,'lf',all_stats{s}))
            stats{2}=stat;
            
            % plot stat.stat, add contours for sig clusters (black) and white contours
            % for non sig
            % plot results
            fig=figure

                for i=1:numel(stats)
                 mask_pos=zeros(size(squeeze(stats{i}.stat)));
                mask_pos_sig=zeros(size(squeeze(stats{i}.stat)));
                mask_pos_trend=zeros(size(squeeze(stats{i}.stat)));
                mask_neg=zeros(size(squeeze(stats{i}.stat)));
                mask_neg_sig=zeros(size(squeeze(stats{i}.stat)));
                mask_neg_trend=zeros(size(squeeze(stats{i}.stat)));
                    
                    % get sig cluster
                  %  [stats{i}.posclusters(:).prob]
                    if isfield(stats{i},'posclusters')
                        if ~isempty(stats{i}.posclusters)
                            if any([stats{i}.posclusters(:).prob]<=0.1)
                                sig_ind=find([stats{i}.posclusters(:).prob]<=0.05);
                                if ~isempty(sig_ind)
                                    mask_pos_sig=squeeze(stats{i}.posclusterslabelmat<=sig_ind(end)&stats{i}.posclusterslabelmat>0);
                                else
                                    sig_ind=0;
                                end
                                trend_ind=find([stats{i}.posclusters(:).prob]<=0.1);
                                mask_pos_trend=squeeze(stats{i}.posclusterslabelmat<=trend_ind(end)&stats{i}.posclusterslabelmat>sig_ind&stats{i}.posclusterslabelmat>0);
                                
                                mask_pos=squeeze(stats{i}.posclusterslabelmat>sig_ind(end));
                            else
                            end
                        else
                        end
                    else
                    end
                    
                    
                    
                  %  [stats{i}.negclusters(:).prob
                    if isfield(stats{i},'negclusters')
                        if ~isempty(stats{i}.negclusters)
                            if any([stats{i}.negclusters(:).prob]<=0.1)
                                sig_ind=find([stats{i}.negclusters(:).prob]<=0.05);
                                if ~isempty(sig_ind)
                                    mask_neg_sig=squeeze(stats{i}.negclusterslabelmat<=sig_ind(end) & stats{i}.negclusterslabelmat>0);
                                else
                                    sig_ind=0;
                                end
                                trend_ind=find([stats{i}.negclusters(:).prob]<=0.1);
                                mask_neg_trend=squeeze(stats{i}.negclusterslabelmat<=trend_ind(end) ...
                                                     & stats{i}.negclusterslabelmat>sig_ind(end)...
                                                     & stats{i}.negclusterslabelmat>0);
                                
                                mask_neg=squeeze(stats{i}.negclusterslabelmat>sig_ind(end));
                            else
                            end
                        else
                        end
                    else
                    end
                    
                    % contour sig
                    mask_alpha=mask_pos_sig+mask_neg_sig;
                    % contour no sig
                    nosig_alpha=mask_pos+mask_neg;
                    trend_alpha=mask_pos_trend+mask_neg_trend;
                    subplot(2,1,i)
                    
                    H= imagesc(stats{i}.time,stats{i}.freq,squeeze(stats{i}.stat),[-5 5])
                    colormap(jet_grey2)
                    set(gca,'YDir','normal')
                    %set(H,'AlphaData',mask_alpha)
                    hold on
                    contour(stats{i}.time,stats{i}.freq,mask_alpha,1,'LineColor','k','LineWidth',1)
                    contour(stats{i}.time,stats{i}.freq,nosig_alpha,1,'LineColor','w','LineWidth',1)
                    contour(stats{i}.time,stats{i}.freq,trend_alpha,1,':','LineColor','k','LineWidth',1)
                    
                    title([sel_roi, ':', cond_name])
                    
                    clear mask_alpha nosig_alpha mask_neg_sig mask_neg mask_pos mask_pos_sig mask_pos_trend mask_neg_trend
                    
                end
                savefig(fig,fullfile(path_fig,[sel_roi,'_',all_stats{s}(1:end-4),'.fig']),'compact')
                close all
            end
        end
    end
