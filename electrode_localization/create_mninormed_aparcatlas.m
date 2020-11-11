%% create a aparc atlas normalized to mni

path_data='D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub03\fieldtrip\c_sub03'
load(fullfile(path_data,'atlasDK.mat'))

mr=ft_read_mri(fullfile(path_data,'T1_fs.nii'))
cfg=[];
cfg.nonlinear ='no';
norm_mr=ft_volumenormalise(cfg,mr)
atlasDK.coordsys='ras';
atlasDK.anatomy=atlasDK.aparc;

cfg=[];
cfg.spmparams        = norm_mr.params;
norm_atlasDK=ft_volumenormalise(cfg,atlasDK)
norm_atlasDK.anatomy=round(norm_atlasDK.anatomy);
norm_atlasDK.aparc=round(norm_atlasDK.aparc);

atlas_2mr=rmfield(norm_atlasDK,'aparc')

ft_write_mri('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.nii',norm_atlasDK,'dataformat','nifti_spm')
norm_atlasDK.aparclabel=atlasDK.aparclabel;
save('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.mat','norm_atlasDK')