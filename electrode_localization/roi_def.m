%% regions of interest for analysis

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

% plot individual freesurfer atlas
load('D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub01\fieldtrip\c_sub01\atlasDK.mat')
mri=ft_read_mri('D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub01\fieldtrip\c_sub01\T1_fs.nii')
mri.coordsys='acpc';
colorcube=colorcube;
mri.aparc=atlasDK.aparc
atlasDK.coordsys='acpc';
cfg.method        = 'ortho';
cfg.anaparameter  = 'aparc';
cfg.atlas=atlasDK;
cfg.coordsys  = 'acpc';
cfg.funparameter ='aparc';
cfg.funcolormap=[colorcube(1:2:end,:);colorcube(2:2:end,:)];
ft_sourceplot(cfg,mri)