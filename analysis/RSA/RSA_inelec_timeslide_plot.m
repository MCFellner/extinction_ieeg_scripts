addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% plot timeslide result in each channel

% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_rsa='D:\Extinction\iEEG\analysis\rsa\';
% path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
% 
% 
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% 
% contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
%     'cs_specific','cs_specific_block1','cs_specific_block2',...
%     'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};
% feature='powlogscale';
% norm='z_crosstrials';
% toi=[2 4];
% win_pow=0.05;
% for c=1:numel(contrasts)
%     contrast=contrasts{c};
%     folder_in=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
%     path_out=fullfile(folder_in,'fig');
%     mkdir(path_out)
%     for sub=1:length(allsubs)
%         sel_sub=allsubs{sub};
%         sel_folder=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
%         load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
%         cfg_rsa.contrast_mat=getfield(contrast_def,contrast);
%         cfg_rsa.sortind=contrast_def.sortind_org2usedtrlinfo;
%         % electrodeinfo
%         info_file=strcat(path_info,sel_sub,'_datainfo');
%         load(info_file)
%         load(fullfile(folder_in,strcat(sel_sub,'_rsastat')),'all_stat')
%         
%         num_chan=numel(all_stat);
%         fig_row=ceil(sqrt(num_chan));
%        fig= figure
%         for chan=1:num_chan
%             subplot(fig_row,fig_row,chan)
%             imagesc(all_stat{chan}.time,all_stat{chan}.freq,squeeze(all_stat{chan}.stat),[-5 5])
%             sel_label=all_stat{chan}.label{1};
%             [~,label_ind]= intersect(datainfo.elec_info.bipolar.elec_ct_mr.label,sel_label,'stable');
%             hold on
%             contour(all_stat{chan}.time,all_stat{chan}.freq,squeeze(all_stat{chan}.mask),1,'k')
%             set(gca,'YDir','normal')
%             title([sel_label,datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{label_ind}]);
%         end
%         savefig(fig,[path_out,'\',sel_sub],'compact')
%         close all
%     end
% end


%% get sig elec summary table


% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_rsa='D:\Extinction\iEEG\analysis\rsa\';
% path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
% 
% 
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
% roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% %roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% %roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% roi.amy_r={'Right-Amygdala'};
% roi.amy_l={'Left-Amygdala'};
% roi.hip_l={'Left-Hippocampus'};
% roi.hip_r={'Right-Hippocampus'};
% roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
% %roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
% 
% rois=fieldnames(roi);
% 
% 
% contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
%     'cs_specific','cs_specific_block1','cs_specific_block2',...
%     'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};
% feature='powlogscale';
% norm='z_crosstrials';
% toi=[2 4];
% win_pow=0.05;
% for c=1:numel(contrasts)
%     contrast=contrasts{c};
%     folder_in=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
%     path_out=fullfile(folder_in,'fig');
%     mkdir(path_out)
%     
% sig_ind=[];
% all_label=[];
% all_pos=[];
%     for sub=1:length(allsubs)
%         sel_sub=allsubs{sub};
%         sel_folder=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
%         load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
%         cfg_rsa.contrast_mat=getfield(contrast_def,contrast);
%         cfg_rsa.sortind=contrast_def.sortind_org2usedtrlinfo;
%         % electrodeinfo
%         info_file=strcat(path_info,sel_sub,'_datainfo');
%         load(info_file)
%         load(fullfile(folder_in,strcat(sel_sub,'_rsastat')),'all_stat')
%         
% %%%%%%%%%%%%%%%%%%%%     
% 
%         sig_tmp_pos=zeros(numel(all_stat),1); %vec 1/0/-1 for sig electrodes
%         sig_tmp_neg=zeros(numel(all_stat),1); %vec 1/0/-1 for sig electrodes
%         sig_tmp=zeros(numel(all_stat),1);
%         sel_label=cell(numel(all_stat),1);
%         numel(all_stat)
%         for chan=1:numel(all_stat)
%             sel_stat=all_stat{chan};
%             sel_label{chan}=sel_stat.label{1};
%             if isfield(sel_stat,'posclusters')
%                 if ~isempty(sel_stat.posclusters)
%                 sig_tmp_pos(chan)=sel_stat.posclusters(1).prob<0.05;
%                 else
%                 sig_tmp_pos(chan)=0;
%                 end
%             else
%                 sig_tmp_pos(chan)=0;
%             end
%                 if isfield(sel_stat,'negclusters')
%                 if ~isempty(sel_stat.negclusters)
%                 sig_tmp_neg(chan)=sel_stat.negclusters(1).prob<0.05;
%                 else
%                 sig_tmp_neg(chan)=0;
%                 end
%             else
%                 sig_tmp_neg(chan)=0;
%             end
%             
%         end
% 
%     sig_tmp(sig_tmp_neg==1)=-1;
%     sig_tmp(sig_tmp_pos==1)=1;
% 
%     sig_def{sub}=sig_tmp';
%     all_elec{sub}=sel_label;
%         clear sel_label
% % get positions        
% [~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');     
% all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
% all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,1);
% 
% sig_ind=[sig_ind;sig_tmp];
% all_label=[all_label;all_elec{sub}];
% all_pos=[all_pos;all_elec_pos{sub}];
% end
% 
% [region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def,roi)
% save(fullfile(folder_in,'results_table'),'region_count','subject_count','roi_count')
% 
% %
% num_all=sum([region_count.absnumelecinregion{:}])
% num_sig=sum([region_count.absnumsigelec_pos{:}])
% rel_sig=sum([region_count.absnumsigelec_pos{:}])/sum([region_count.absnumelecinregion{:}])
% 
% clear num_all num_sig rel_sig region_count 
% % plot electrodes and count electrodes per region/pat
% % sort elecs for subfunction in sig and no sig group
% sort_ind=sig_ind;
% all_ind=unique(sort_ind);
% for i=1:numel(all_ind)
%  elec_pos{i}=all_pos(sig_ind==all_ind(i),:);
%  elec_label{i}=all_label(sig_ind==all_ind(i));
% end
% 
% cfg.elec_pos =elec_pos;
% cfg.elec_label=elec_label;
% cfg.col_def=[0,0,1;0.5 0.5 0.5;1,0,0];
% cfg.trans_def=[1 0.2 1];
% views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
% views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
% cfg.views=views;
% cfg.hemispheres={'left','right'};
% cfg.hemisphere_surf={'D:\matlab_tools\fieldtrip-20200130\template\anatomy\surface_pial_left.mat',...
%     'D:\matlab_tools\fieldtrip-20200130\template\anatomy\surface_pial_right.mat'};
% 
% cfg.path_fig=folder_in;
% mcf_ieegelectrodeplotter(cfg)
% clear cfg elec_pos elec_label all_ind sort_ind sig_def all_elec_pos all_elec_label
% end

%% get average t-map for each roi

% plot average
% plot all electrodes in roi

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_rsa='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';


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


contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
    'cs_specific','cs_specific_block1','cs_specific_block2',...
    'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};
feature='powlogscale';
norm='z_crosstrials';
toi=[2 4];
win_pow=0.05;

load('D:\matlab_tools\jet_grey.mat')
for c=1:numel(contrasts)
    contrast=contrasts{c};
    folder_in=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
    path_out=fullfile(folder_in,'fig');
    mkdir(path_out)
    
sig_ind=[];
all_label=[];
all_pos=[];
    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        sel_folder=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
        cfg_rsa.contrast_mat=getfield(contrast_def,contrast);
        cfg_rsa.sortind=contrast_def.sortind_org2usedtrlinfo;
        % electrodeinfo
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
        load(fullfile(folder_in,strcat(sel_sub,'_rsastat')),'all_stat')
        
sel_label=[];
   for chan=1:numel(all_stat)
             sel_label{chan}=all_stat{chan}.label{1};
   end
 
        [~,~,ind]=intersect(sel_label,datainfo.elec_info.bipolar.elec_mni.label,'stable');     
        all_elec_label=[datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{ind,1}];
        
   % for every roi check for electrodes and average t-map     
        for r=1:numel(rois)
            sel_rois=getfield(roi,rois{r});
            sel_rois={sel_rois{:}};
            [~,~,label_ind]=intersect( all_elec_label,sel_rois,'stable');  
            sel_roi_elec=[];
            for l=1:numel(label_ind)
           sel_roi_elec=[sel_roi_elec,find(strcmp(all_elec_label,sel_rois{label_ind(l)}))];
            end
            
            %all_stats{sub,r}=all_stat(sel_roi_elec);
            sel_tmaps=[];
            if ~isempty(sel_roi_elec)
            num_elec=numel(sel_roi_elec)
            sel_tmaps=zeros([num_elec,size(squeeze(all_stat{1}.stat))]);
             for e=1:num_elec
             sel_tmaps(e,:,:)=all_stat{e}.stat;
             %sel_masks(e,:,:)=all_stat{e}.mask;

             end
             avg_tmap{sub,r}=squeeze(nanmean(sel_tmaps,1));
            end
            
        end
    end
    
 t_num=size(squeeze(all_stat{1}.stat),1);
   figure

% plot average tmap
for r=1:numel(rois)
   sel_rois=rois{r};
   sel_tmap=reshape([avg_tmap{:,r}],t_num,t_num,[]);
   subplot(3,3,r)
   imagesc(all_stat{1}.time,all_stat{1}.time,squeeze(nanmean(sel_tmap,3)),[-4 4])
   colormap(jet_grey)
   colorbar
   title([contrast,':',sel_rois,': mean t'])
   ylabel('t in s')
   xlabel('t in s')
end
end