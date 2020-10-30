%% check location of remaining electrodes: bip processed data checking bip labels (freesurferDK, nearest GM)
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};

roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.ifg=[roi.ifg_r,roi.ifg_l];

roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
roi.dm_pfc=[roi.dm_pfc_r,roi.dm_pfc_l];

roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};

roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
roi.ventraltempocci=[roi.ventraltempocci_l,roi.ventraltempocci_r];


% load a freesurfer atlas to get all labels
load('D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub01\fieldtrip\c_sub01\atlasDK.mat')
freesurfer_label=atlasDK.aparclabel;
for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        % electrodeinfo
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
              
        % channels surviving preprocessing
        clean_chan=datainfo.artifact_info.rejectvisual_bip.elecsin';
        
        % get indices of channel in analabel
        [chan,~,ind]=intersect(clean_chan,datainfo.elec_info.bipolar.elec_mni.label,'stable');
        % check whether each channel is matched
        if numel(chan)~=numel(clean_chan)
            error('channel labels do not match')
        end       
   
        all_labels=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,:);
        all_labels=[all_labels{:}]';
            
        % goal: count of how many subjects have at least one contact in an area    
        for r=1:numel(freesurfer_label)
           count_region(r,sub)=sum(strcmp(freesurfer_label{r},all_labels)); 
        end

end
 
all_ind= find(sum(count_region,2));
count_region=count_region(all_ind,:);
all_labels=freesurfer_label(all_ind);

% put data in table
all_labels_count=table(all_labels, count_region,sum(count_region>0,2),'VariableNames',{'regions','numelecinpat','numpatperregion'});
% add roi count
all_rois=fieldnames(roi);
all_labels=[all_labels(:)];

for r=1:numel(all_rois)
    sel_def=all_rois{r};
    sel_label=getfield(roi,sel_def)
[tmp,ind]=intersect(all_labels,sel_label','stable' )
elec_in_roi(r,:)=sum(count_region(ind,:),1);

end

roi_count=table(all_rois,elec_in_roi,sum(elec_in_roi>0,2),'VariableNames',{'regions','numelecpat','numpatperregion'});

file_out=fullfile(path_info,strcat('countregion_bipafterpreproc'));
save(file_out,'roi_count','roi','all_labels_count')