addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
%% create a aparc atlas normalized to mni

path_data='D:\Extinction\iEEG\data\freesurfer_anat\output\c_sub03\fieldtrip\c_sub03'
load(fullfile(path_data,'atlasDK.mat'))

mr=ft_read_mri(fullfile(path_data,'T1_fs.nii'))

mr.coordsys='acpc';
cfg=[];
norm_mr=ft_volumenormalise(cfg,mr)

atlasDK.coordsys='acpc';
atlasDK.anatomy=atlasDK.aparc;
atlasDK.anatomylabel=atlasDK.aparclabel;
%%% every region needs to be normalized on its own...
% loop trough every label, normalized mask, define the region using value >1
atlasDK_tmp=atlasDK;
atlasDK_tmp=rmfield(atlasDK_tmp, 'aparc');
atlasDK_tmp=rmfield(atlasDK_tmp, 'aparclabel');

tmp_anatomy=zeros(size(norm_mr.anatomy));
for r=1:numel(atlasDK.anatomylabel)
       cfg=[];
    cfg.atlas=atlasDK;
    cfg.inputcoord='acpc';
    cfg.roi=atlasDK.anatomylabel{r};
    tmp_mask=ft_volumelookup(cfg,atlasDK);
    atlasDK.anatomylabel{r}
    
atlasDK_tmp.anatomy=tmp_mask;
atlasDK_tmp.anatomylabel=atlasDK.anatomylabel(r);

%  create a mask for this region    

cfg=[];
cfg.spmparams        = norm_mr.params;
tmp_norm=ft_volumenormalise(cfg,atlasDK_tmp)
% combine masks to anatomy
tmp_anatomy(tmp_norm.anatomy>0)=r;

end
% atlasDK.anatomy=tmp_anatomy;
% atlasDK.aparc=tmp_anatomy;

atlasDK.anatomy=atlasDK.aparc;
atlasDK=rmfield(atlasDK,'aparc');
cfg=[];
cfg.spmparams        = norm_mr.params;
norm_atlasDK=ft_volumenormalise(cfg,atlasDK)

%%%%%%% round is not appropiate(nonsense values)
% norm_atlasDK.anatomy=round(norm_atlasDK.anatomy);
% norm_atlasDK.aparc=round(norm_atlasDK.aparc);

norm_atlasDK.anatomy=tmp_anatomy;
norm_atlasDK.aparc=tmp_anatomy;
norm_atlasDK.aparclabel=atlasDK.aparclabel;
norm_atlasDK.anatomylabel=atlasDK.aparclabel;



ft_write_mri('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.nii',norm_atlasDK,'dataformat','nifti_spm')
norm_atlasDK.aparclabel=atlasDK.aparclabel;
save('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.mat','norm_atlasDK')


%% check fsaverage

%%% normalizing fsaverage may be the better solution:
path_in='D:\Extinction\iEEG\rawdata\China\data_bids_forfmriprep\derivatives\freesurfer\fsaverage\mri'

filename=fullfile(path_in,'aparc+aseg.mgz')

standard_atlasDK = ft_read_mri(filename, 'dataformat','freesurfer_mgz')
filename_out='D:\Extinction\iEEG\rawdata\China\data_bids_forfmriprep\derivatives\freesurfer\fsaverage\mri\test.nii'
ft_write_mri(filename_out,standard_atlasDK,'dataformat', 'nifti')
%%%% missing here: norm fs mri us transform on atlas for aparc in mni space
% maybe coregister enough?
atlasDK=ft_read_atlas(filename, 'dataformat','freesurfer_mgz')
