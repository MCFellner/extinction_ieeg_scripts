%% query datainfo for different regions
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

% setup definition

sel_atlas='aal';
sel_regions={'Amygdala_L','Hippocampus_L','Rectus_L','Frontal_Sup_Orb_L','Frontal_Med_Orb_L','Cingulum_Ant_L',...
            'Amygdala_R','Hippocampus_R','Rectus_R','Frontal_Sup_Orb_R','Frontal_Med_Orb_R','Cingulum_Ant_R'};
distance_region=3; % distance of electrodes from region

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
          'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};

% load datainfo for each subject
for n=1:numel(allsubs)
sel_sub=allsubs{n};  

% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

% query datainfo 
cfg.atlas=sel_atlas;
cfg.region= sel_regions;
cfg.distance_region= distance_region;
elec_selection=mcf_elec_region_selector(cfg,datainfo);

% sum elec in each region
sum_per_region(n,:)=elec_selection.count;
clear elec_selection elec_region_mat 
end
