%% create a aparc atlas normalized to mni

path_data='D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub01\fieldtrip\c_sub01'
load(fullfile(path_data,'atlasDK.mat'))

mr=ft_read_mri(fullfile(path_data,'T1_fs.nii'))
cfg=[];
norm_mr=ft_volumenormalise(cfg,mr)
atlasDK.coordsys='ras';
atlasDK.anatomy=atlasDK.aparc;

cfg=[];
cfg.spmparams        = norm_mr.params;
norm_atlasDK=ft_volumenormalise(cfg,atlasDK)

atlas_2mr=rmfield(norm_atlasDK,'aparc')
ft_write_mri('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.nii',norm_atlasDK,'dataformat','nifti_spm')

norm_atlasDK.aparclabel=atlasDK.aparclabel;
save('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.mat','norm_atlasDK')