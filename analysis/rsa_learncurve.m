% rsa learning curves

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')
%%
%
% path_rsa='D:\Extinction\iEEG\analysis\rsa\';
% path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% toi=[3 4];
% win=0.1; % in sec, power estimated for every win
% feature='powlogscale';
% norm='z_crosstrials';
% sel_folder=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
% nrand=1000;
%
% contrasts={'learn_cond1_block1','learn_cond1_block2','learn_cond1_block3',...
%             'learn_cond2_block1','learn_cond2_block2','learn_cond2_block3',...
%             'learn_cond3_block1','learn_cond3_block2','learn_cond3_block3'};
%
% for c=1:numel(contrasts)
% contrast=contrasts{c};
%
% for sub=1:length(allsubs)
%     sel_sub=allsubs{sub};
% % load rsa matrix
% load(fullfile(sel_folder,strcat(sel_sub,'_rsa_')))
% load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
%
% contrast_mat=getfield(contrast_def,contrast);
%
% % sort rsa to match contrast_mat
% sortind=contrast_def.sortind_org2usedtrlinfo;
% rsa.trlinfo=rsa.trlinfo(sortind,:);
% rsa.rsa_mat=rsa.rsa_mat(:,sortind,:);
% rsa.rsa_mat=rsa.rsa_mat(:,:,sortind);
%
%
% num_trial=size(rsa.trlinfo,1);
% num_chan=numel(rsa.label);
%
% contrast_vec=reshape(contrast_mat,[],1);
%
%
% rsa_tmp=permute(rsa.rsa_mat,[2,3,1]);
% rsa_tmp=reshape(rsa_tmp,[],num_chan);
%
% sel_ind=find(contrast_vec>0&~isnan(contrast_vec));
% sel_rsa=rsa_tmp(sel_ind,:);
% sel_reg=contrast_vec(sel_ind);
%
%
% sel_data.time=[1 1.1];
% sel_data.dimord='rpt_chan';
% %sel_data.avg=repmat(sel_rsa,1,1,2);
% sel_data.avg=zscore(sel_rsa);
% sel_data.label=rsa.label;
%
%     cfg=[];
%
%                 cfg.avgoverchan =  'no';
%                 cfg.method           = 'analytic';
%                 cfg.statistic ='indepsamplesregrT';
% %                 cfg.numrandomization = 1000;
% %                 cfg.correctm         =  'cluster';
% %                 cfg.correcttail      = 'prob';
%                 cfg.design = sel_reg;
%                 cfg.ivar=1;
%                 cfg.cvar=[];
%                 stat=ft_timelockstatistics (cfg,sel_data);
%
%
% stat.rsa=zscore(sel_rsa);
% stat.reg=sel_reg;
%
% folder_out=fullfile(sel_folder,contrast);
% mkdir(folder_out)
% save(fullfile(folder_out,strcat(sel_sub,'_rsa_learnstat')),'stat')
% clear stat
% end
% end
%
% %%
% path_rsa='D:\Extinction\iEEG\analysis\rsa\';
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% toi=[3 4];
% win=0.1; % in sec, power estimated for every win
% feature='powlogscale';
% norm='z_crosstrials';
% nrand=1000;
%
% contrasts={'learn_cond1_block1','learn_cond1_block2','learn_cond1_block3',...
%     'learn_cond2_block1','learn_cond2_block2','learn_cond2_block3',...
%     'learn_cond3_block1','learn_cond3_block2','learn_cond3_block3'};
% roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
% roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};
% roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
% roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
% roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
%
% roi.amy_r={'Right-Amygdala'};
% roi.amy_l={'Left-Amygdala'};
% roi.hip_l={'Left-Hippocampus'};
% roi.hip_r={'Right-Hippocampus'};
%
% roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
% roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%
% path_in=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
% alpha_def=0.05;
% for c=1:numel(contrasts)
%     contrast=contrasts{c};
%     sel_path=fullfile(path_in,contrast);
%     sig_ind=[];
%     all_label=[];
%     all_pos=[];
%     for sub=1:length(allsubs)
%         sel_sub=allsubs{sub};
%         load(fullfile(sel_path,strcat(sel_sub,'_rsa_learnstat')))
%
%         info_file=strcat(path_info,sel_sub,'_datainfo');
%         load(info_file)
%
%         sig_tmp=(stat.prob<alpha_def).*sign(stat.stat);
%
%         sig_def{sub}=sig_tmp';
%         all_elec{sub}=stat.label;
%
%         % get positions
%         [~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');
%         all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
%         all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,1);
%
%         sig_ind=[sig_ind;sig_tmp];
%         all_label=[all_label;stat.label];
%         all_pos=[all_pos;all_elec_pos{sub}];
%     end
%
%     [region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def,roi)
% save(fullfile(sel_path,'results_table'),'region_count','subject_count','roi_count')
%
% end
% clear
%
%% combine data from different rois for rfx analysis
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
% 
% rois=fieldnames(roi);
% 
% path_rsa='D:\Extinction\iEEG\analysis\rsa\';
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% 
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% toi=[2 3];
% win=0.1; % in sec, power estimated for every win
% feature='powlogscale';
% norm='z_crosstrials';
% nrand=1000;
% 
% contrasts={'learn_cond1_block1','learn_cond1_block2','learn_cond1_block3',...
%     'learn_cond2_block1','learn_cond2_block2','learn_cond2_block3',...
%     'learn_cond3_block1','learn_cond3_block2','learn_cond3_block3'};
% 
% % two conditions to contrast
% %contrasts={'learn_cond1_block1','learn_cond3_block1'};
% 
% 
% path_in=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
% alpha_def=0.05;
% % get electrodes for each sub/roi
% 
% % get electrodes for each sub/roi
% for r=1:numel(rois)
%     sel_roi=rois{r};
%     roi_def=getfield(roi,sel_roi);
%     sub_num=0;
%     for sub=1:length(allsubs)
%         sel_sub=allsubs{sub};
%         % electrodeinfo
%         info_file=strcat(path_info,sel_sub,'_datainfo');
%         load(info_file)
%         sel_elec_tmp=datainfo.elec_info.bipolar.elec_ct_mr.label(ismember([datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{:}],[roi_def]));
%         % only select after preproc elecs
%         sel_elec=intersect(sel_elec_tmp,datainfo.artifact_info.rejectvisual_bip.elecsin);
% 
%         roi_sub(r,sub)=numel(sel_elec)>0;
%         roi_sub_elec{r,sub}=sel_elec;
%     end
% end
% 
% for r=1:numel(rois)
%     sel_roi=rois{r};
% 
%     sel_allsubs={allsubs{roi_sub(r,:)}};
%     elec_def=[roi_sub_elec(r,roi_sub(r,:))];
% 
% 
% for c=1:numel(contrasts)
%     contrast=contrasts{c};
%     sel_path=fullfile(path_in,contrast);
% 
%     for sub=1:length(sel_allsubs)
%         sel_sub=sel_allsubs{sub};
%         load(fullfile(sel_path,strcat(sel_sub,'_rsa_learnstat')))
%         sel_elec=elec_def{sub};
%         [~,elec_ind]=intersect(stat.label,sel_elec,'stable');
%         avg_t(sub,c)=nanmean(stat.stat(elec_ind));
%     end
% 
% end
% 
% 
% 
% [h(r,:),p(r,:),~,ts]=ttest(avg_t)
% t_stat(r,:)=ts.tstat;
% 
% end
% figure
% ax1=subplot(1,2,1)
% imagesc(p,[0 0.1])
% colormap(hot)
% tmp=colormap;
% colormap(ax1,tmp(65-(1:size(tmp,1)),:))
% xticks(1:numel(contrasts))
% xticklabels(contrasts)
% xtickangle(60)
% yticks(1:numel(rois))
% yticklabels(rois)
% title('p-value, sig regression in roi, rfx')
% colorbar
% ax2=subplot(1,2,2)
% imagesc(t_stat,[-5 5])
% load('D:\matlab_tools\jet_grey.mat')
% colormap(ax2,jet_grey)
% xticks(1:numel(contrasts))
% xticklabels(contrasts)
% xtickangle(60)
% yticks(1:numel(rois))
% yticklabels(rois)
% title('t-value, sig regression in roi, rfx')
% hold on
% contour(1:numel(contrasts),1:numel(rois),p<0.05,1,'LineColor','k','LineWidth',1)
% colorbar
%% plot average curve in roi


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

path_rsa='D:\Extinction\iEEG\analysis\rsa\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

toi=[2 3];
win=0.1; % in sec, power estimated for every win
feature='powlogscale';
norm='z_crosstrials';
nrand=1000;

contrasts={'learn_cond1_block1','learn_cond1_block2','learn_cond1_block3',...
    'learn_cond2_block1','learn_cond2_block2','learn_cond2_block3',...
    'learn_cond3_block1','learn_cond3_block2','learn_cond3_block3'};

% two conditions to contrast
%contrasts={'learn_cond1_block1','learn_cond3_block1'};


path_in=fullfile(path_rsa,strcat(feature,'_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)));
alpha_def=0.05;
% get electrodes for each sub/roi

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
        
        roi_sub(r,sub)=numel(sel_elec)>0;
        roi_sub_elec{r,sub}=sel_elec;
    end
end

contrasts={'learn_cond1_block1','learn_cond1_block2','learn_cond1_block3';...
    'learn_cond2_block1','learn_cond2_block2','learn_cond2_block3';...
    'learn_cond3_block1','learn_cond3_block2','learn_cond3_block3'};

for r=1:numel(rois)
    sel_roi=rois{r};
    
    sel_allsubs={allsubs{roi_sub(r,:)}};
    elec_def=[roi_sub_elec(r,roi_sub(r,:))];
    
    all_rsa=nan(numel(sel_allsubs),3,61);
    
    for sub=1:length(sel_allsubs)
        sel_sub=sel_allsubs{sub};
        sel_elec=elec_def{sub};
        % load rsa matrix
        load(fullfile(path_in,strcat(sel_sub,'_rsa_')))
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
        [~,elec_ind]=intersect(rsa.label,sel_elec,'stable');
        sortind=contrast_def.sortind_org2usedtrlinfo;
        rsa.trlinfo=rsa.trlinfo(sortind,:);
        rsa_mat=rsa.rsa_mat(:,sortind,:);
        rsa_mat=rsa_mat(:,:,sortind);
        rsa_mat=reshape(rsa_mat,numel(rsa.label),[]);
        
        for c=1:size(contrasts,1)
            contrast=contrasts(c,:);
            
            
            for b=1:numel(contrast)
                sel_contrast=contrast{b};
                contrast_mat=getfield(contrast_def,sel_contrast);
                
                
                [~,col]=find(contrast_mat>0);
                trl_ind=rsa.trlinfo(col,15);
                if b==2
                    trl_ind=trl_ind+23;
                elseif b==3
                    trl_ind=trl_ind+46
                else
                end
                sel_ind=find(reshape(contrast_mat,[],1)>0);
                sel_rsa=nanmean(rsa_mat(elec_ind,sel_ind));
                
                all_rsa(sub,c,trl_ind)=sel_rsa;
                
                
            end
        end
    end
    f=  figure
    fig=subplot(1,2,1)
    cmap_default=fig.ColorOrder;
    
    close all
    fig=figure
    subplot(6,5,1:25)
    for c=1:3
        x1=1:61;
        y1=squeeze(nanmean(all_rsa(:,c,:)));
        b1=squeeze(nanstd(all_rsa(:,c,:))./sqrt(size(all_rsa,1)));
        boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
        
    end
    title([sel_roi,': rsa trial to trial sim'])
    plot([23,23],[0,0.1],'k:')
    plot([46,46],[0,0.1],'k:')
    
    subplot(6,5,26)
    hold on
    plot(squeeze(nanmean(all_rsa(:,:,:)))')
    axis off
    legend({'plusplus','plusminus','minusminus'});
    
    savefig(fig,fullfile(path_in,[sel_roi,'unsmoothed_learncurve']),'compact')
    close all
    
    all_rsa=smoothdata(all_rsa,3,'movmean',5)
    
           % dimord trl_pow 'subj_freq_rpt_cond_timewindow'
            data1.dimord='subj_chan_time';
            data1.label={sel_roi};
            data1.time=1:61;
            data2=data1;
            data3=data1;
            data1.avg(:,1,:)=squeeze(all_rsa(:,1,:));
            data2.avg(:,1,:)=squeeze(all_rsa(:,2,:));
            data3.avg(:,1,:)=squeeze(all_rsa(:,3,:));
    
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
                cfg.statistic ='depsamplesFmultivariate';
                 Nsub = size(data1.avg,1);                                       %# of subjects?
                design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub) 3*ones(1,Nsub)];
                design(2,:)  = [1:Nsub 1:Nsub 1:Nsub];
   cfg.clusteralpha     = 0.05;
                cfg.uvar     = 2;
                cfg.ivar     = 1;
                cfg.design = design;
                cfg.tail=1
    stat= ft_timelockstatistics(cfg,data1, data2,data3);
    
            clear design data1 data2 data3
    
    
    f=  figure
    fig=subplot(1,2,1)
    cmap_default=fig.ColorOrder;
    
    close all
    fig=figure
    subplot(6,5,1:25)
    scatter(find(stat.mask),zeros(size(find(stat.mask))),'*')
    hold on
    for c=1:3
        x1=1:61;
        y1=squeeze(nanmean(all_rsa(:,c,:)));
        b1=squeeze(nanstd(all_rsa(:,c,:))./sqrt(size(all_rsa,1)));
        boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
        
    end
    title([sel_roi,': rsa trial to trial sim'])
    plot([23,23],[0,0.1],'k:')
    plot([46,46],[0,0.1],'k:')
    
    subplot(6,5,26)
    hold on
    plot(squeeze(nanmean(all_rsa(:,:,:)))')
    axis off
    legend({'plusplus','plusminus','minusminus'});
    
    savefig(fig,fullfile(path_in,[sel_roi,'smoothed_learncurve']),'compact')
end