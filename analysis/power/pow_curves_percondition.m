addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')



% plot power curves for each conditon in each roi, timewindow & freq


%% first step: freq per trial and timewindow for every sub
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
% path_out='D:\Extinction\iEEG\analysis\pow\avg_pow\';
% mkdir(path_out)
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% % segment data in different trial parts
% % item window: -1 to 4 (imp
% pre_item=-1;
% post_item=5.5;
% toi=2:0.5:5;
% fois=[2 4;5 7;8 12;13 20;20 30; 35  60; 70 100];
% foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};
%
% % downsample to smallest sr
% sr=1000;
%
%
% for sub=27%1:length(allsubs)
%     sel_sub=allsubs{sub};
%     % electrodeinfo
%     info_file=strcat(path_info,sel_sub,'_datainfo');
%     load(info_file)
%
%     load(strcat(path_preproc,sel_sub,'_data.mat'))
%
%     % apply bipolar montage
%     montage=datainfo.elec_info.bipolar.montage_withoutartichan;
%     data = ft_apply_montage(data,montage);
%
%     % cut trials
%     trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
%     trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
%     trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
%     cfg=[];
%     cfg.trl=trl;
%     data=ft_redefinetrial(cfg,data);
% clear trl
%     % downsample to common sampling rate
%     cfg=[];
%     cfg.resamplefs      = sr;
%     cfg.detrend='yes';
%     data=ft_resampledata(cfg,data);
%
%     % only select artifree trials in trialinfo
%     trlinfo=datainfo.trialinfo;
%     trlinfo=trlinfo(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree,:);
%
%     cfg=[]
%     cfg.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
%     data=ft_selectdata(cfg,data);
%
%     for f=1:length(fois)
%         foi=fois(f,:);
%         f_name=foi_names{f};
%          sel_path=fullfile(path_out,f_name);
%         mkdir(sel_path)
%
%         cfg=[];
%         cfg.output='pow';
%         cfg.toi =pre_item:0.05:post_item;
%         cfg.keeptrials  = 'yes';
%         cfg.pad='nextpow2';
%         switch f_name
%             case {'theta1','theta2','alpha','beta1','beta2'}
%                 cfg.foilim     = foi;
%                 cfg.method='wavelet';
%                 cfg.width = 5;
%             case {'gamma1','gamma2'}
%                 cfg.method      = 'mtmconvol';
%                 cfg.taper       = 'dpss';
%                 cfg.foi         = foi(1):5:foi(2);
%                 cfg.t_ftimwin   = ones(1,length(cfg.foi))*0.3;
%                 cfg.tapsmofrq   = ones(1,length(cfg.foi))*10;
%         end
%         freq=ft_freqanalysis(cfg,data);
%
%         % z trans
%         cfg=[];
%         cfg.time=[pre_item+0.5,post_item-0.5];
%         freq=z_trans_TF_seltime(cfg, freq);
%         freq.trialinfo=trlinfo;
%
%         % average across freq
%         freq.powspctrm=nanmean(freq.powspctrm,3);
%         freq.freq=nanmean(freq.freq);
%
%         % average for tois
%         tmp_pow=zeros(size(trlinfo,1),numel(freq.label),1,numel(toi)-1);
%         for t=1:(numel(toi)-1)
%             t1=nearest(freq.time,toi(t));
%             t2=nearest(freq.time,toi(t+1))-1;
%             tmp_pow(:,:,:,t)=nanmean(freq.powspctrm(:,:,:,t1:t2),4);
%             tmp_time(t)=mean([toi(t) toi(t+1)]);
%         end
%         freq.powspctrm=tmp_pow;
%         freq.time=tmp_time;
%
%         % save (folder for each freq)
%
%         save(fullfile(sel_path,[sel_sub,'_avgpow']),'freq')
%
%     end
%
% end

%% load freq avg files and plot power curves

path_in='D:\Extinction\iEEG\analysis\pow\avg_pow\';
path_fig='D:\Extinction\iEEG\analysis\pow\avg_pow\fig\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

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
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};
conditions={'plusplus','plusminus','minusminus'};
%conditions={'us','nous'};


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



for r=1:numel(rois)
    sel_roi=rois{r};
    
    sel_allsubs={allsubs{roi_sub(r,:)}};
    sel_path=fullfile(path_in,foi_names{1});
    elec_def=[roi_sub_elec(r,roi_sub(r,:))];
    load(fullfile(sel_path,[sel_allsubs{1},'_avgpow']));
    
    trl_pow=nan(numel(sel_allsubs),numel(foi_names),64,3,numel(freq.time));
    for f=1:numel(foi_names)
        sel_freq=foi_names{f};
        sel_path=fullfile(path_in,sel_freq);
        
        
        for sub=1:length(sel_allsubs)
            sel_sub=sel_allsubs{sub};
            sel_elec=elec_def{sub};
            load(fullfile(sel_path,[sel_sub,'_avgpow']));
            [~,elec_ind]=intersect(freq.label,sel_elec,'stable');
            
            for c=1:numel(conditions)
                sel_trl=  freq.trialinfo(:,6)==c;
                sel_trlinfo=freq.trialinfo(sel_trl,:);
                
                % fix indices of trial reps
                sel_ind=[sel_trlinfo(sel_trlinfo(:,2)==1,15);...
                    sel_trlinfo(sel_trlinfo(:,2)==2,15)+24;...
                    sel_trlinfo(sel_trlinfo(:,2)==3,15)+48];
                sel_pow=squeeze(nanmean(freq.powspctrm(sel_trl,elec_ind,:,:),2));
                trl_pow(sub,f,sel_ind,c,:)=sel_pow;
            end
        end
        
        
    end

    trl_pow_sm=smoothdata(trl_pow,3,'movmean',7);
    
    avg_pow=squeeze(nanmean(trl_pow_sm));
    f=  figure
    fig=subplot(1,2,1)
    cmap_default=fig.ColorOrder;
    close (f)
   fig= figure
    num_f=numel(foi_names);
    num_t=size(trl_pow,5);
    
    for t=1:num_t
    subplot(num_f+2,num_t+1,t+1)
        axis off
        text(0,0, [num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'])
    end
    for f=1:num_f
        subplot(num_f+2,num_t+1,((f)*(num_t+1))+1)
        axis off
        text(0,0, foi_names{f})
    end
    
    for f=1:num_f

        for t=1:num_t
            % cluster f test
            
                % dimord trl_pow 'subj_freq_rpt_cond_timewindow'
            data1.dimord='subj_chan_time';
            data1.label={sel_roi};
            data1.time=1:64;
            data2=data1;
            data3=data1;
            data1.avg(:,1,:)=squeeze(trl_pow_sm(:,f,:,1,t));
            data2.avg(:,1,:)=squeeze(trl_pow_sm(:,f,:,2,t));
            data3.avg(:,1,:)=squeeze(trl_pow_sm(:,f,:,3,t));
    
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
            
            subplot(num_f+2,num_t+1,((f).*(num_t+1))+t+1)
            hold on
            if any(stat.mask)
            scatter(find(stat.mask),zeros(size(find(stat.mask))),'*')
            end
            for c=1:numel(conditions)
                % tmp= smoothdata(avg_pow(f,:,c,t),'movmean',3);
                tmp= avg_pow(f,:,c,t);
                x1=1:64;
                y1=tmp;
                b1=squeeze(nanstd(trl_pow_sm(:,f,:,c,t))./sqrt(size(trl_pow,1)));
                boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
                
                % plot(tmp)
            end
            plot([24,24],[0,0.5],'k:')
            plot([48,48],[0,0.5],'k:')
         %   title([sel_roi,': ',foi_names{f},' ',num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'])
        end
    end
    subplot(num_f+1,num_t+1,(num_f.*(num_t+1))+1)
    hold on
    for c=1:numel(conditions)
        tmp=squeeze(avg_pow(f,:,c,t));
        plot(tmp)
    end
    legend(conditions)
    sel_path=fullfile(path_fig,sel_roi);
    mkdir(sel_path)
    savefig(fig,[sel_path,'\smoothed_timecourse_sm7'],'compact')
    close all
    
    %barplot: blockab x type
    
    fig=figure
    num_f=numel(foi_names);
    num_t=size(trl_pow,5);
    num_sub=size(trl_pow,1);
    for f=1:num_f
        for t=1:num_t
            subplot(num_f+1,num_t,((f-1).*num_t)+t)            
            hold on            
            sel_pow=squeeze(trl_pow(:,f,:,:,t));                  
        % for block x type block a and block b
            sel_pow_block(:,:,1)=nanmean(sel_pow(:,1:24,:),2);
            sel_pow_block(:,:,2)=nanmean(sel_pow(:,25:48,:),2);           
            tmp=reshape(sel_pow_block,size(sel_pow,1),[]);
            data = table(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),tmp(:,5),tmp(:,6),...
                'VariableNames',{'V1','V2','V3','V4','V5','V6'});
            design(:,1) = {'A','A','A','B','B','B'};
            design(:,2) = {'1','2','3','1','2','3'};      
            
            within=table(design(:,1),design(:,2),'VariableNames',{'block','type'})
            rm = fitrm(data,'V1-V6~1','WithinDesign',within)
            
            ranovatbl = ranova(rm,'WithinModel','block*type')
            hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3])
           xticklabels(conditions)
            scatter([0.85.*ones(num_sub,1);1.15.*ones(num_sub,1);...
            1.85.*ones(num_sub,1);2.15.*ones(num_sub,1);...
            2.85.*ones(num_sub,1);3.15.*ones(num_sub,1)],...
            reshape(tmp,1,[]),3,'ko','filled');
            
            %xlabel(['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})])
                 
            title({[sel_roi,': ',foi_names{f},' ',num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'];['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})]})
            
        end
    end
    subplot(num_f+1,num_t,((f).*num_t)+1)            
      hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3])
           xticklabels(conditions)
            legend({'blockA','blockB'}) 
            
    sel_path=fullfile(path_fig,sel_roi);
    mkdir(sel_path)
    savefig(fig,[sel_path,'\blockAB_type'],'compact')
    clear design sel_pow_block
    close all
    
    
    
     %barplot: block a: halfs x type
    
    fig=figure
    num_f=numel(foi_names);
    num_t=size(trl_pow,5);
    num_sub=size(trl_pow,1);
    for f=1:num_f
        for t=1:num_t
            subplot(num_f+1,num_t,((f-1).*num_t)+t)            
            hold on            
            sel_pow=squeeze(trl_pow(:,f,:,:,t));                  
        % for type x halfes
            sel_pow_block(:,:,1)=nanmean(sel_pow(:,1:12,:),2);
            sel_pow_block(:,:,2)=nanmean(sel_pow(:,13:24,:),2);           
            tmp=reshape(sel_pow_block,size(sel_pow,1),[]);
            data = table(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),tmp(:,5),tmp(:,6),...
                'VariableNames',{'V1','V2','V3','V4','V5','V6'});
            design(:,1) = {'first','first','first','second','second','second'};
            design(:,2) = {'1','2','3','1','2','3'};      
            
            within=table(design(:,1),design(:,2),'VariableNames',{'halfs','type'})
            rm = fitrm(data,'V1-V6~1','WithinDesign',within)
            
            ranovatbl = ranova(rm,'WithinModel','halfs*type')
            hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3])
           xticklabels(conditions)
            scatter([0.85.*ones(num_sub,1);1.15.*ones(num_sub,1);...
            1.85.*ones(num_sub,1);2.15.*ones(num_sub,1);...
            2.85.*ones(num_sub,1);3.15.*ones(num_sub,1)],...
            reshape(tmp,1,[]),3,'ko','filled');
            
           % xlabel(['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})])
                  
            title({[sel_roi,': ',foi_names{f},' ',num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'];['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})]})
            
        end
    end
    subplot(num_f+1,num_t,((f).*num_t)+1)            
      hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3])
           xticklabels(conditions)
            legend({'firsthalf','secondhalf'}) 
            
    sel_path=fullfile(path_fig,sel_roi);
    mkdir(sel_path)
    savefig(fig,[sel_path,'\blockA_halfs_type'],'compact')
    clear design sel_pow_block
    close all   
   
      %barplot: block b: halfs x type
    
    fig=figure
    num_f=numel(foi_names);
    num_t=size(trl_pow,5);
    num_sub=size(trl_pow,1);
    for f=1:num_f
        for t=1:num_t
            subplot(num_f+1,num_t,((f-1).*num_t)+t)            
            hold on            
            sel_pow=squeeze(trl_pow(:,f,:,:,t));                  
        % for type x halfes
            sel_pow_block(:,:,1)=nanmean(sel_pow(:,25:36,:),2);
            sel_pow_block(:,:,2)=nanmean(sel_pow(:,37:48,:),2);           
            tmp=reshape(sel_pow_block,size(sel_pow,1),[]);
            data = table(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),tmp(:,5),tmp(:,6),...
                'VariableNames',{'V1','V2','V3','V4','V5','V6'});
            design(:,1) = {'first','first','first','second','second','second'};
            design(:,2) = {'1','2','3','1','2','3'};      
            
            within=table(design(:,1),design(:,2),'VariableNames',{'halfs','type'})
            rm = fitrm(data,'V1-V6~1','WithinDesign',within)
            
            ranovatbl = ranova(rm,'WithinModel','halfs*type')
            hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3])
           xticklabels(conditions)
            scatter([0.85.*ones(num_sub,1);1.15.*ones(num_sub,1);...
            1.85.*ones(num_sub,1);2.15.*ones(num_sub,1);...
            2.85.*ones(num_sub,1);3.15.*ones(num_sub,1)],...
            reshape(tmp,1,[]),3,'ko','filled');
            
         %   xlabel(['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})])
                  
            title({[sel_roi,': ',foi_names{f},' ',num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'];['block p=',num2str(ranovatbl{3,5}),' type p=',num2str(ranovatbl{5,5}),'int p=',num2str(ranovatbl{7,5})]})
            
        end
    end
    subplot(num_f+1,num_t,((f).*num_t)+1)            
      hold on
            bar(reshape(nanmean(tmp),[],2))
            xticks([1 2 3]);
           xticklabels(conditions);
            legend({'firsthalf','secondhalf'}) 
            
    sel_path=fullfile(path_fig,sel_roi);
    mkdir(sel_path)
    savefig(fig,[sel_path,'\blockB_halfs_type'],'compact')
    clear design sel_pow_block
    close all   
end

%% sanity check: contrast us vs no us

path_in='D:\Extinction\iEEG\analysis\pow\avg_pow\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

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
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};
conditions={'us','nous'};


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



for r=1:numel(rois)
    sel_roi=rois{r};
    
    sel_allsubs={allsubs{roi_sub(r,:)}};
    sel_path=fullfile(path_in,foi_names{1});
    elec_def=[roi_sub_elec(r,roi_sub(r,:))];
    load(fullfile(sel_path,[sel_allsubs{1},'_avgpow']));
    
    trl_pow=nan(numel(sel_allsubs),numel(foi_names),2,numel(freq.time));
    for f=1:numel(foi_names)
        sel_freq=foi_names{f};
        sel_path=fullfile(path_in,sel_freq);
        
        
        for sub=1:length(sel_allsubs)
            sel_sub=sel_allsubs{sub};
            sel_elec=elec_def{sub};
            load(fullfile(sel_path,[sel_sub,'_avgpow']));
            [~,elec_ind]=intersect(freq.label,sel_elec,'stable');
            
            for c=0:1
                sel_trl=  freq.trialinfo(:,2)==1& freq.trialinfo(:,9)==c;
                sel_trlinfo=freq.trialinfo(sel_trl,:);
                
                sel_pow=squeeze(nanmean(nanmean(freq.powspctrm(sel_trl,elec_ind,:,:),1),2));
                trl_pow(sub,f,c+1,:)=sel_pow;
            end
        end
        
        
    end
    
    figure
    num_f=numel(foi_names);
    num_t=size(trl_pow,4);
    for f=1:num_f
        for t=1:num_t
            subplot(num_f+1,num_t,((f-1).*num_t)+t)
            
            
            tmp= squeeze(trl_pow(:,f,1:2,t));
            [h,p,~]=ttest(tmp(:,1),tmp(:,2));
            bar(nanmean(tmp))
            hold on
            scatter([ones(size(tmp,1),1);ones(size(tmp,1),1).*2],reshape(tmp,1,[]))
            xticklabels({'us','nous'})
            if h
                text(1.5,0.1,'*')
            end
            title([sel_roi,': ',foi_names{f},' ',num2str(freq.time(t).*1000-250),'to',num2str(freq.time(t).*1000+250),'ms'])
        end
    end
end