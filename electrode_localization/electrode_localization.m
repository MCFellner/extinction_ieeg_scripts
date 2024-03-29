
% fieldtrip ieeg electrode localization pipeline
%addpath('D:\matlab_tools\fieldtrip-20190611')
addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% read in mri
path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

for sub=1:numel(allsubs)

sel_sub=allsubs{sub};   

sel_folder=strcat(path_data,sel_sub,'\anat\');
cd (sel_folder)
mri = ft_read_mri('T1.nii'); 

% description from tutorial:
    % check the coordinate system of the mri
        % important left-right definition:
            % if left to right increases to the right: scan left-to-right
            % orientation

        %     1 Visualize the coordinate system of the MRI or CT by inputting the 
        %     following command: ft_determine_coordsys(mri) 
        %     2 Determine which of the three axes, x, y, or z, runs through or 
        %     along the left�right axis of the subject�s head. This axis is the 
        %     left�right axis for this anatomical volume. 
        %     3 Determine the orientation of the left�right axis. If the values on 
        %     the left�right axis increase to the right (indicated by a �+� sign), 
        %     then the scan has a left-to-right orientation. If the values on the 
        %     left�right axis increase to the left, then the scan has a right-to-left
        %     orientation. 
        %     4 Write down the orientation of the scan�s left�right axis.
%%%%% our mr seem flipped (only kind of a good match in the coregistration
%%%%% later)

ft_determine_coordsys(mri);

% align mri to acpc (for freesurfer): remember orientation: L= l+
display('unflipped version: mark here right according to the orientation marked before')
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);

% here: mark r on the "wrong" side
display('flipped version: mark here right on the other side compared to before')
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'acpc';
mri_acpc_flipped = ft_volumerealign(cfg, mri)

% for some reason this mri has only negative values, shift values to
% positive ones
% check for negative values in data

min_anatomy=min(min(min(mri_acpc.anatomy)));
if min_anatomy<0
    mri_acpc.anatomy=mri_acpc.anatomy+abs(min_anatomy);
    mri_acpc_flipped.anatomy=mri_acpc_flipped.anatomy+abs(min_anatomy);
end

% write it to file
cfg           = [];
cfg.filename  = [sel_sub '_MR_acpc'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy'
ft_volumewrite(cfg, mri_acpc);

cfg           = [];
cfg.filename  = [sel_sub '_MR_acpc_flipped'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy'
ft_volumewrite(cfg, mri_acpc_flipped);

close all 
keep allsubs sub path_data
end
%% 
path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};


for sub=5:numel(allsubs)
sel_sub=allsubs{sub}   

sel_folder=strcat(path_data,sel_sub,'\anat\');
cd (sel_folder)

% import ct
file_ct=dir('CT2.*');
ct = ft_read_mri(file_ct.name); % we used the dcm series
%
%ft_determine_coordsys(ct); % here check with ppt in info folder

% realign to ctf (nasion, lpa, rpa and interhemispheric location)
cfg           = [];
cfg.method    = 'interactive';
%cfg.coordsys  = 'ctf';
cfg.coordsys  = 'acpc';
ct_ctf = ft_volumerealign(cfg, ct);
ct_acpc=ct_ctf;
% same coordinates as the mri
%ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');

% cfg = [];
% cfg.anaparameter = 'anatomy';
% cfg.funparameter = 'anatomy';
% cfg.location = [0 0 60];
% ft_sourceplot(cfg, ct_acpc)

%import  processed mri
filename=strcat(sel_sub,'_MR_acpc.nii') % change here to nii from freesurfer
fsmri_acpc = ft_read_mri(filename); 
fsmri_acpc.coordsys = 'acpc';

% fusion of ct and mri: unflipped
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'acpc';
cfg.viewresult  = 'yes';
cfg.parameter='anatomy';
ct_acpc_f = ft_volumerealign(cfg,ct_acpc,fsmri_acpc);

% verify fusion!
% for verifying fusion: also use mricron (better scrollable)
% if fusion sucks: redo all steps, think about whether there might be a
% l/r flip in the data

% write fused mri to file
filename=strcat(sel_sub,'_ct_acpc_f') 
cfg           = [];
cfg.filename  = filename;
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);

% fusion ct and mri: flipped
filename=strcat(sel_sub,'_MR_acpc_flipped.nii') % change here to nii from freesurfer
fsmri_acpc_flipped = ft_read_mri(filename); 
fsmri_acpc_flipped.coordsys = 'acpc';

cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'acpc';
cfg.viewresult  = 'yes';
cfg.parameter='anatomy';
ct_acpc_f = ft_volumerealign(cfg,ct_acpc,fsmri_acpc_flipped);

% write fused mri to file
filename=strcat(sel_sub,'_ct_acpc_f_toflippedMR') % change here to nii from freesurfer
cfg           = [];
cfg.filename  = filename;
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);
close all 
keep allsubs sub path_data
end

%%%%%%%%%%%% % verify fusion!
% for verifying fusion:  use mricron (better scrollable)
% decide whether mr flipped or unflipped

%%%%%%%%%%%%%%%

%% define correct T1: flipped or unflipped and move to freesurfer exchange folder
path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

for i=1:numel(allsubs)
sel_sub=allsubs{i}    
load(strcat(path_info,sel_sub,'_datainfo.mat'))

flip_def=input('Is original T1 flipped or not? write ''flipped'' or ''no_flip''  ');
datainfo.anat.mri.filename='T1.nii';
datainfo.anat.mri.flipped=flip_def;
save(strcat(path_info,sel_sub,'_datainfo.mat'),'datainfo')

%fs_path=strcat('I:\Extinction\iEEG\data\eeg\freesurfer_anat\mri\',sel_sub,'\');
fs_path=strcat('D:\Extinction\iEEG\data\freesurfer_anat\mri\',sel_sub,'\');
mkdir(fs_path)
if strcmp(flip_def,'flipped')
filename=strcat(path_data,sel_sub,'\anat\',sel_sub,'_MR_acpc_flipped.nii');
elseif strcmp(flip_def,'no_flip')
filename=strcat(path_data,sel_sub,'\anat\',sel_sub,'_MR_acpc.nii');
else
    error('flip not defined')
end

file_out=strcat(fs_path,sel_sub,'_MR_acpc.nii');
copyfile(filename,file_out)
end
      
%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start linux virtual box 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % cortical surface extraction: freesurfer
% % 
% % % % needs to be done using linux
% % % sel_sub='c_sub01';
% % % path_data='/media/freesurfer/Elements/Extinction/iEEG/data/eeg/freesurfer_anat/';
% % % 
% % % 
% % % fshome     = '/usr/local/freesurfer';
% % % subdir     = strcat(path_data,'output/',sel_sub,'/');
% % % mrfile     = strcat(path_data,'mri/',sel_sub,'/',sel_sub,'_MR_acpc.nii');
% % % mkdir(subdir)
% % % 
% % % system(['export FREESURFER_HOME=' fshome '; ' ...
% % % 'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
% % % 'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
% % % 'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])
% % % 
% % % %% export mgz and atlas to windows readable files
% % % 
% % % addpath('/media/freesurfer/Elements/matlab_tools/fieldtrip-20190108')
% % % ft_defaults
% % % 
% % % sel_sub='c_sub01';
% % % 
% % % path_data='/media/freesurfer/Elements/Extinction/iEEG/data/eeg/freesurfer_anat/output/';
% % % 
% % % path_in=strcat(path_data,sel_sub,'/freesurfer/mri/');
% % % path_out=strcat(path_data,sel_sub,'/fieldtrip/',sel_sub,'/');
% % % 
% % % mkdir(path_out)
% % % cd(path_in)
% % % 
% % % % save atlas as matlab file for readability in windows
% % % atlasDK=ft_read_atlas('aparc+aseg.mgz');
% % % save(strcat(path_out,'atlasDK.mat'),'atlasDK')
% % % clear atlasDK
% % % atlasDKT40=ft_read_atlas('aparc.DKTatlas+aseg.mgz');
% % % save(strcat(path_out,'atlasDKT40.mat'),'atlasDKT40')
% % % clear atlasDKT40
% % % atlasDest=ft_read_atlas('aparc.a2009s+aseg.mgz');
% % % save(strcat(path_out,'atlasDest.mat'),'atlasDest')
% % % clear atlasDest
% % % 
% % % %save fs mri (slightly changes dimensions...)
% % % mri=ft_read_mri('T1.mgz')
% % % cfg=[];
% % % cfg.parameter     = 'anatomy';
% % % cfg.filename      = 'T1_fs';
% % % cfg.filetype      = 'nifti'
% % % ft_volumewrite(cfg,mri)
% % % clear mri
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linux/virtual box end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%% move files back to folders
% path_data='D:\Extinction\iEEG\data\freesurfer_anat\';
% fs_path=strcat('I:\Extinction\iEEG\data\eeg\freesurfer_anat\');
% 
% copyfile(fs_path,path_data)
% 

%for new files
fs_path='D:\Extinction\iEEG\rawdata\China\data_bids_forfmriprep\derivatives\freesurfer'
path_data='D:\Extinction\iEEG\data\freesurfer_anat\output';
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};
bidssubs = {'sub-c19','sub-c21','sub-c22','sub-c23','sub-c24','sub-c25','sub-c26','sub-c29'};

for sub=2:numel(allsubs)
 sel_bids=fullfile(fs_path,bidssubs{sub});
 sel_out=fullfile(path_data,allsubs{sub},'freesurfer');
 mkdir(sel_out)
 copyfile(sel_bids,sel_out)   
end

%% export fs data to windows friendly files (for new bids data)
path_data='D:\Extinction\iEEG\data\freesurfer_anat\output';
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};


for sub=1:numel(allsubs)
sel_sub=allsubs{sub};

path_in=fullfile(path_data,sel_sub,'freesurfer','mri');
path_out=fullfile(path_data,sel_sub,'fieldtrip',sel_sub);

mkdir(path_out)
cd(path_in)

% save atlas as matlab file for readability in windows
atlasDK=ft_read_atlas('aparc+aseg.mgz');
save(fullfile(path_out,'atlasDK.mat'),'atlasDK')
clear atlasDK
atlasDKT40=ft_read_atlas('aparc.DKTatlas+aseg.mgz');
save(fullfile(path_out,'atlasDKT40.mat'),'atlasDKT40')
clear atlasDKT40
atlasDest=ft_read_atlas('aparc.a2009s+aseg.mgz');
save(fullfile(path_out,'atlasDest.mat'),'atlasDest')
clear atlasDest
%save fs mri (slightly changes dimensions...)
mri=ft_read_mri('T1.mgz')
cfg=[];
cfg.parameter     = 'anatomy';
cfg.filename      = fullfile(path_out,'T1_fs');
cfg.filetype      = 'nifti'
ft_volumewrite(cfg,mri)
clear mri
end

%% coregister freesurfer mr to ct
path_data='D:\Extinction\iEEG\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

for i=1:numel(allsubs)
sel_sub=allsubs{i};    
path_mr=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub);
path_ct=strcat(path_data,'rawdata\extinction_ieeg\',sel_sub,'\anat\');

cd (path_ct)

%import ct
file_ct=dir('CT2.*');
ct = ft_read_mri(file_ct.name);
%
%ft_determine_coordsys(ct); % here check with ppt in info folder, are electrodes in the right hemisphere? (pos=right)

% realign to ctf (nasion, lpa, rpa and interhemispheric location)
cfg           = [];
cfg.method    = 'interactive';
%cfg.coordsys  = 'ctf';
%ct_ctf = ft_volumerealign(cfg, ct);
cfg.coordsys='acpc';
ct_acpc=ft_volumerealign(cfg, ct);


% % same coordinates as the mri
% ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');
% 
% % 
% cfg = [];
% cfg.anaparameter = 'anatomy';
% cfg.funparameter = 'anatomy';
% cfg.location = [0 0 0];
% ft_sourceplot(cfg, ct_acpc)

% cfg = [];
% cfg.anaparameter = 'anatomy';
% cfg.funparameter = 'anatomy';
% %cfg.location = [0 0 60];
% ft_sourceplot(cfg, fsmri_acpc)


%import  processed mri
cd(path_mr)
fsmri_acpc = ft_read_mri('T1_fs.nii'); 
fsmri_acpc.coordsys = 'acpc';

% fusion of ct and mri: unflipped
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'acpc';
cfg.viewresult  = 'yes';
cfg.parameter='anatomy';
ct_acpc_f = ft_volumerealign(cfg,ct_acpc,fsmri_acpc);

% write fused mri to file
filename=strcat(sel_sub,'_ct_acpc_f') % change here to nii from freesurfer
cfg           = [];
cfg.filename  = filename;
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);
end

%% write electrodeinfo files; check names and numbers

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

for i=4%:numel(allsubs)
sel_sub=allsubs{i};  
% electrode names from ieeg files
 eeg_file=strcat(path_data,'data\preproc\ieeg\readin\', sel_sub,'_data.mat');
load(eeg_file)
 info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

elec_info.labels_orderindata=data.label;
% resort labels for easier labeling
% define electrodes
ekg_count=0;
for i=1:numel(data.label)
    sel_lab=data.label{i};
    if strncmp(sel_lab, 'EKG',3)
    channel_type{i}='EKG';
    electrode{i,1}='EKG';
    ekg_count=ekg_count+1;
    electrode{i,2}=ekg_count;
    else
    channel_type{i}='iEEG';
    for j=1:numel(sel_lab)
    num_ind(j)=str2double(sel_lab(j));
    end
    electrode{i,1}=sel_lab(isnan(num_ind));
    electrode{i,2}=str2double(sel_lab(~isnan(num_ind)));
    clear num_ind
    end
end

cell_ieeg=sortrows(electrode(strcmp(channel_type,'iEEG'),:));

for i=1:size(cell_ieeg,1)
   sorted_ieeg{i,1}=strcat(cell_ieeg{i,1},num2str(cell_ieeg{i,2}));
end
elec_info.orderindata.label=data.label;
elec_info.orderindata.channel_type=channel_type;
elec_info.reordered.cell_ieeg=cell_ieeg;
elec_info.reordered.sorted_ieeg=sorted_ieeg;

datainfo.elec_info=elec_info;
save(info_file,'datainfo')
keep path_data path_info allsubs

end
% %% copy relevant files to Q (dataexchange)
% 
%% datainfo, ct, mr
% path_data='D:\Extinction\iEEG\';
% path_dataout='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\electrode_localization\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
% 
% for i=1:numel(allsubs)
% sel_sub=allsubs{i};  
% path_mr=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
% path_info=strcat(path_data,'data\preproc\ieeg\datainfo\');
% path_out=strcat(path_dataout,sel_sub,'\');
% mkdir(path_out)
% 
% info_file=strcat(path_info,sel_sub,'_datainfo.mat');
% ct_file=strcat(path_mr,sel_sub,'_ct_acpc_f.nii')
% mr_file=strcat(path_mr,'T1_fs.nii');
% 
% info_out=strcat(path_out,sel_sub,'_datainfo.mat');
% ct_out=strcat(path_out,sel_sub,'_ct_acpc_f.nii')
% mr_out=strcat(path_out,sel_sub,'_T1_fs.nii');
% 
% copyfile(info_file,info_out)
% copyfile(ct_file,ct_out)
% copyfile(mr_file,mr_out)
% 
% end



%% mark electrode in coregistered ct files

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

i=8%:numel(allsubs)
sel_sub=allsubs{i};  

% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

%load coregistered ct &mr
path_mr=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub);
cd(path_mr)
ct_file=strcat(sel_sub,'_ct_acpc_f.nii');
ct = ft_read_mri(ct_file);
mr_file=strcat('T1_fs.nii');
mr = ft_read_mri(mr_file);
ct.coordsys='acpc';
mr.coordsys='acpc';

cfg         = [];
cfg.channel = datainfo.elec_info.reordered.sorted_ieeg;
if isfield(datainfo.elec_info,'elec_ct_mr')
cfg.elec           =datainfo.elec_info.elec_ct_mr;
end
elec_ct_mr = ft_electrodeplacement(cfg, ct, mr);
datainfo.elec_info.elec_ct_mr=elec_ct_mr;
save(info_file,'datainfo')
keep path_data path_info allsubs
close all
delete(gcf)


% 



%% get electrode information for paris data & move necessary files in localization folder
% path_data='D:\Extinction\iEEG\';
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% 
% allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% for i=8:numel(allsubs)
% sel_sub=allsubs{i};  
% path_anat=strcat(path_data,'\rawdata\extinction_ieeg\',sel_sub,'\anat\epiloc\');
% path_fsout=strcat(path_data,'\data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');
% mkdir(path_fsout)
% 
% path_fs=strcat(path_anat,'segmentations\ref_t1mri\freesurfer\');
% path_mr=strcat(path_anat,'nifti\');
% path_ct=strcat(path_anat,'reg_ima\')
% 
% % get mr scan 
% cd(path_mr)
% mr_file=dir('*_t1mri.nii.gz'); % this is the nifti on which the freesurfer segmentation is run
% mr = ft_read_mri(mr_file.name);
% 
% cd(path_ct)
% ct_file=dir('*_ct_post_2_t1mri.nii.gz'); % this is the ct coregistered to the mr
% ct = ft_read_mri(ct_file.name);
% 
% cd(path_fsout)
% cfg=[];
%  cfg.parameter     = 'anatomy';
%  cfg.filename      = 'T1_fs';
%  cfg.filetype      =  'nifti';
% ft_volumewrite(cfg,mr)
% cfg=[];
%  cfg.parameter     = 'anatomy';
%  cfg.filename      = strcat(sel_sub,'_ct_acpc_f.nii');
%  cfg.filetype      =  'nifti';
% ft_volumewrite(cfg,ct)
% 
% % get fs atlas
% cd(path_fs)
% atlasDK=ft_read_atlas('aparc+aseg.nii.gz');
% save(strcat(path_fsout,'atlasDK.mat'),'atlasDK')
% clear atlasDK
% 
% atlasDest=ft_read_atlas('aparc.a2009s+aseg.nii.gz');
% save(strcat(path_fsout,'atlasDest.mat'),'atlasDest')
% clear atlasDest
% 
% % write electrodeinfo in datainfo
% load(strcat(path_data,'data\preproc\ieeg\readin\',sel_sub,'_data.mat'))
% load(strcat(path_info,sel_sub,'_datainfo.mat'))
% elec_info.orderindata.label=data.label;
% elec_info.orderindata.channel_type=data.hdr.chantype;
% 
% cd(path_anat)
% 
% csv_file=dir('*tmp.csv');
% [file_elec{1},file_elec{2},...
%     file_mni{1},file_mni{2},file_mni{3},...
%     file_mrmm_ac{1},file_mrmm_ac{2},file_mrmm_ac{3},...
%     file_mrvox_ac{1},file_mrvox_ac{2},file_mrvox_ac{3},...
%     file_mrmm{1},file_mrmm{2},file_mrmm{3},...
%     file_ctmm_ac{1},file_ctmm_ac{2},file_ctmm_ac{3},...
%     file_ctvox_ac{1},file_ctvox_ac{2},file_ctvox_ac{3},...
%     ]=textread(csv_file.name,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*s %*s %*s %*s %*s %*s %*s %*s','delimiter',';','headerlines',2);
% 
% % there seem to be empty last rows:
% if isempty(file_elec{1}{end})
% file_elec{1}(end)=[];
% end
%         for i=1:numel(file_elec{1})    
%         elec{i,1}=strcat(file_elec{1}{i},'_',file_elec{2}{i});
%         ct_mr_pos(i,:)=[str2double(file_mrmm{1}{i}),str2double(file_mrmm{2}{i}),str2double(file_mrmm{3}{i})];
%         mni_pos(i,:)=[str2double(file_mni{1}{i}),str2double(file_mni{2}{i}),str2double(file_mni{3}{i})];
%         end
% clear file_elec file_mni file_mrmm_ac file_mrvox_ac file_mrmm file_ctmm_ac file_ctmm_ac file_ctvox
% 
% % check match of elecs and data.label
% if strcmp(sel_sub,'p_sub02')
%   for i=1:numel(elec)
%       if strncmp(elec{i},'H2ag',4)
%         elec{i}=strcat('Ha2g_',elec{i}(end));
%       end
%   end
% elseif strcmp(sel_sub,'p_sub03')
%   for i=1:numel(elec)
%       if strncmp(elec{i},'Col',3)
%         elec{i}=strcat('Col1_',elec{i}(end));
%       elseif strncmp(elec{i},'HaT',3)
%         elec{i}=strcat('HaT2_',elec{i}(end));
%       end
%   end
% elseif strcmp(sel_sub,'p_sub05')
%   for i=1:numel(elec)
%       if strncmp(elec{i},'AmT',3)
%         elec{i}=strcat('AmT2_',elec{i}(end));
%       elseif strncmp(elec{i},'TPol',4)
%         elec{i}=strcat('Tpol_',elec{i}(end));
%       end
%   end
%   elseif strcmp(sel_sub,'p_sub08')
%   for i=1:numel(data.label)
%       if strncmp(data.label{i},'SiPa',4)
%         data.label{i}=strcat('SiPa_',data.label{i}(end));
%       end
%       if strncmp(data.label{i},'T1p',3)& numel(data.label{i})==4
%         data.label{i}=strcat('T1p_',data.label{i}(end));
%       end
%       if strncmp(data.label{i},'T2as',4)
%         data.label{i}=strcat('T2as_',data.label{i}(end));
%       end      
%       if strncmp(data.label{i},'GyAn',4)
%         data.label{i}=strcat('GyAn_',data.label{i}(end));
%       end
%       if strncmp(data.label{i},'T2ai',4)
%         data.label{i}=strcat('T2ai_',data.label{i}(end));
%       end    
%       if strncmp(data.label{i},'T2mi',4)
%         data.label{i}=strcat('T2mi_',data.label{i}(end));
%       end 
%       if strncmp(data.label{i},'T3p',3)& numel(data.label{i})==4
%         data.label{i}=strcat('T3p_',data.label{i}(end));
%       end 
%       if strncmp(data.label{i},'T2ps',4)
%         data.label{i}=strcat('T2ps_',data.label{i}(end));
%       end
%       if strncmp(data.label{i},'T2pi',4)
%         data.label{i}=strcat('T2pi_',data.label{i}(end));
%       end     
%       if strncmp(data.label{i},'T2p',3)& numel(data.label{i})==4
%         data.label{i}=strcat('T2p_',data.label{i}(end));
%       end 
%   end
% end
% elecs_not_in_data=setdiff(elec,data.label)
% elecs_not_in_epiloc=setdiff(data.label,elec)
% check=input('missing labels ok?'); % quick check if labels are matching
% 
% 
% elec_info.elec_ct_mr.unit='mm';
% elec_info.elec_ct_mr.coordsys='fsaverage';
% elec_info.elec_ct_mr.label=elec;
% elec_info.elec_ct_mr.tra=diag(ones(numel(elec),1));
% elec_info.elec_ct_mr.elecpos=ct_mr_pos;
% elec_info.elec_ct_mr.chanpos=elec_info.elec_ct_mr.elecpos;
% elec_info.elec_ct_mr.cfg.method='epiloc';
% elec_info.elec_ct_mr.cfg.filename=csv_file;
% elec_info.elec_ct_mr.cfg.mrfile=mr_file;
% elec_info.elec_ct_mr.cfg.ctfile=ct_file;
% 
% elec_info.elec_mni.unit='mm';
% elec_info.elec_mni.coordsys='mni';
% elec_info.elec_mni.label=elec;
% elec_info.elec_mni.tra=diag(ones(numel(elec),1));
% elec_info.elec_mni.elecpos=mni_pos;
% elec_info.elec_mni.chanpos=elec_info.elec_mni.elecpos;
% elec_info.elec_mni.cfg.method='epiloc';
% elec_info.elec_mni.cfg.filename=csv_file;
% 
% datainfo.elec_info=elec_info;
% 
% % check electrode label using the ct/mr plotting
% ct.coordsys='fsaverage';
% mr.coordsys='fsaverage';
% cfg         = [];
% cfg.elec           =elec_info.elec_ct_mr;
% elec_tmp = ft_electrodeplacement(cfg, ct, mr);
% 
% check2=input('do position/labels match?')
% clear mr ct elec ct_mr_pos mni_pos
% 
% 
% if strcmp(sel_sub,'p_sub08')
% datainfo.elec_info.orderindata.label=data.label;
% datainfo.elec_info.orderindata.channel_type=data.hdr.chantype;%
% save(strcat(path_data,'data\preproc\ieeg\readin\',sel_sub,'_data.mat'),'data')
% end
% %save(strcat(path_info,sel_sub,'_datainfo.mat'),'datainfo')
% 
% end


 
%% normalized mrs for normalization parameters
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
% 
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

      
for i=1:numel(allsubs)
sel_sub=allsubs{i};  

path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
path_fs=strcat(path_data,'\data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');

cd(path_anat)
 mr_file=strcat('T1_fs.nii');
 mr = ft_read_mri(mr_file);
mr.coordsys='acpc';

% % volume based normalization (to get mni coordinates)
cfg            = [];
cfg.nonlinear  = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod  = 'new';
mr_mni = ft_volumenormalise(cfg, mr);

cfg=[];
 cfg.parameter     = 'anatomy';
 cfg.filename      = 'T1_fs_normed';
 cfg.filetype      =  'nifti';
ft_volumewrite(cfg,mr_mni) % save normed mr to check normalization

save(strcat(sel_sub,'_norm_mri.mat'),'mr_mni')
end
 
%% apply norm param to get mni coordinates
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};


allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};


for i=1:numel(allsubs)
sel_sub=allsubs{i};  

% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
path_fs=strcat(path_data,'\data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');

cd(path_anat)
load(strcat(sel_sub,'_norm_mri.mat'))

% % apply parameters to electrode postions
elec_tmp = datainfo.elec_info.elec_ct_mr;
elec_tmp.elecpos = ft_warp_apply(mr_mni.params,elec_tmp.elecpos, 'individual2sn');
elec_tmp.chanpos = ft_warp_apply(mr_mni.params, elec_tmp.chanpos, 'individual2sn');
elec_tmp.coordsys = 'mni';
elec_tmp.cfg=mr_mni.cfg;
elec_tmp.cfg.file=strcat(path_anat,sel_sub,'_norm_mri.mat');
% 
datainfo.elec_info.elec_mni=elec_tmp;

% save normalized coordinates in datainfo
info_file=strcat(path_info,sel_sub,'_datainfo');
save(info_file,'datainfo')
end
%% get labels using individualized freesurfer coordinates and normalized  coordinates
% 
% use individual freesurfer atlas
% use mni based atlas

% search with different query ranges: 1,3,5,7,9,11 (r=query/2-0.5),
% radius:0,1,2,3,4,5

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
%allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

%allsubs = {'p_sub08'};
allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

all_atlas={'freesurferDestrieux','afni','brainweb','freesurferDK','aal'}%,}
query_range=[1,3,5,7,9,11];

for i=1:numel(allsubs)
% save labels in file & as excel sheet
sel_sub=allsubs{i};
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

for a=1:numel(all_atlas)
%path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');
path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
path_templates='D:\matlab_tools\fieldtrip-20200130\template\atlas\';

sel_atlas=all_atlas{a};

if ~any(strcmp(datainfo.elec_info.elec_ct_mr.label,datainfo.elec_info.elec_mni.label))
error('electrode labels in elec_ct_mr and elect_mni do not match');
end

elec_labels.labels=datainfo.elec_info.elec_ct_mr.label;
elec_labels.query_ranges=query_range;
datainfo.elec_info.ana_labels.labels=datainfo.elec_info.elec_ct_mr.label;
datainfo.elec_info.ana_labels.query_ranges=query_range;

      switch sel_atlas
        case 'freesurferDK'
            elecinfo=datainfo.elec_info.elec_ct_mr;            
            cd(path_anat)
            load('atlasDK.mat')
            atlas=atlasDK;
            atlas.coordsys='fsaverage'; 
            inputcoord='fsaverage';
            datainfo.elec_info.ana_labels.freesurferDK_def.file=strcat(path_anat,'atlasDK.mat');
            datainfo.elec_info.ana_labels.freesurferDK_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.ana_labels.freesurferDK_def.atlas_unit=atlas.unit;
            datainfo.elec_info.ana_labels.freesurferDK_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.ana_labels.freesurferDK_def.elec_def=elecinfo;
        case 'freesurferDestrieux'
            elecinfo=datainfo.elec_info.elec_ct_mr;              
            cd(path_anat)
            load('atlasDest.mat')
            atlas=atlasDest;
            atlas.coordsys='fsaverage'; 
            inputcoord='fsaverage';
            datainfo.elec_info.ana_labels.freesurferDestrieux_def.file=strcat(path_anat,'atlasDest.mat');
            datainfo.elec_info.ana_labels.freesurferDestrieux_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.ana_labels.freesurferDestrieux_def.atlas_unit=atlas.unit;
            datainfo.elec_info.ana_labels.freesurferDestrieux_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.ana_labels.freesurferDestrieux_def.elec_def=elecinfo;
        case 'aal'
             elecinfo=datainfo.elec_info.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'aal\ROI_MNI_V4.nii');
             atlas=ft_read_atlas(atlas_file);
            inputcoord='mni';
            datainfo.elec_info.ana_labels.aal_def.file=atlas_file;
            datainfo.elec_info.ana_labels.aal_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.ana_labels.aal_def.atlas_unit=atlas.unit;
            datainfo.elec_info.ana_labels.aal_def.sphereradius_mm=(query_range-1);
            datainfo.elec_info.ana_labels.aal_def.elec_def=elecinfo;

        case 'afni'
             elecinfo=datainfo.elec_info.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'afni\TTatlas+tlrc.HEAD');
             atlas=ft_read_atlas(atlas_file);
             inputcoord='mni';
            datainfo.elec_info.ana_labels.afni_def.file=atlas_file;
            datainfo.elec_info.ana_labels.afni_def.atlas_resolution=[1 1 1];
            datainfo.elec_info.ana_labels.afni_def.atlas_unit=atlas.unit;
            datainfo.elec_info.ana_labels.afni_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.ana_labels.afni_def.elec_def=elecinfo;
        case 'brainweb'
             elecinfo=datainfo.elec_info.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'brainweb\brainweb_discrete.mat');
             load(atlas_file)
             inputcoord='mni';
            datainfo.elec_info.ana_labels.brainweb_def.file=atlas_file;
            datainfo.elec_info.ana_labels.brainweb_def.atlas_resolution=[1 1 1];
            datainfo.elec_info.ana_labels.brainweb_def.atlas_unit=atlas.unit;
            datainfo.elec_info.ana_labels.brainweb_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.ana_labels.brainweb_def.elec_def=elecinfo;
        otherwise
    end
      
for r=1:numel(query_range)
sel_range=query_range(r);

%parfor possible
parfor e=1:numel(elecinfo.label)
sel_pos=elecinfo.elecpos(e,:);
cfg            = [];
cfg.roi        = sel_pos;
cfg.atlas      = atlas;
cfg.output     = 'label';
cfg.minqueryrange = sel_range;
cfg.maxqueryrange = sel_range;
cfg.inputcoord=inputcoord;
labels = ft_volumelookup(cfg, atlas);

if sum(labels.count)>0
tmp_label(e,r)={labels.name(labels.count>0)'};
else
tmp_label(e,r)={'no label found'};
end

end
end

eval(strcat('datainfo.elec_info.ana_labels.',sel_atlas,'=tmp_label;'));
clear tmp_label
end

save(info_file,'datainfo')
end

%% add info to paris electrode info
% info about electrodes

% % info about the channeltype
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% for n=8%2:numel(allsubs)
% sel_sub=allsubs{n};  
% 
% % electrodeinfo
% info_file=strcat(path_info,sel_sub,'_datainfo');
% load(info_file)
% 
% % there are some electrode inconsistencies (electrodes in epiloc not in
% % data and vice versa), these issues are fixed here
% 
% iEEG_labels=datainfo.elec_info.elec_ct_mr.label;
% inctmr_not_indata=setdiff(datainfo.elec_info.elec_ct_mr.label, datainfo.elec_info.orderindata.label)
% % delete channels  inctmr_not_indata (no need for the elecpos then)
% if ~isempty(inctmr_not_indata)
%      %delete in elec_ct_mr
%      [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.elec_ct_mr.label)
%      datainfo.elec_info.elec_ct_mr.label(ind)=[];
%      datainfo.elec_info.elec_ct_mr.chanpos(ind,:)=[];
%      datainfo.elec_info.elec_ct_mr.elecpos(ind,:)=[];
%      datainfo.elec_info.elec_ct_mr.tra(ind,:)=[];
%      datainfo.elec_info.elec_ct_mr.tra(:,ind)=[];
%      
%     %delete in elec_mni
%      [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.elec_mni.label)
%      datainfo.elec_info.elec_mni.label(ind)=[];
%      datainfo.elec_info.elec_mni.chanpos(ind,:)=[];
%      datainfo.elec_info.elec_mni.elecpos(ind,:)=[];
%      datainfo.elec_info.elec_mni.tra(ind,:)=[];
%      datainfo.elec_info.elec_mni.tra(:,ind)=[];
%     %delete in ana_labels
%      [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.ana_labels.labels)
%      datainfo.elec_info.ana_labels.freesurferDK(ind,:)=[];
%      datainfo.elec_info.ana_labels.freesurferDestrieux(ind,:)=[];
%      datainfo.elec_info.ana_labels.aal(ind,:)=[];
%      datainfo.elec_info.ana_labels.afni(ind,:)=[];
%      datainfo.elec_info.ana_labels.brainweb(ind,:)=[];
%      
%       datainfo.elec_info.ana_labels.freesurferDK_def.elec_def=datainfo.elec_info.elec_ct_mr;
%       datainfo.elec_info.ana_labels.freesurferDestrieux_def.elec_def=datainfo.elec_info.elec_ct_mr;
%       datainfo.elec_info.ana_labels.aal_def.elec_def=datainfo.elec_info.elec_mni;
%       datainfo.elec_info.ana_labels.afni_def.elec_def=datainfo.elec_info.elec_mni;
%       datainfo.elec_info.ana_labels.brainweb_def.elec_def=datainfo.elec_info.elec_mni;
%     end
%   
% 
% indata_not_inctmr=setdiff(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label)
% % define channels indata_not_inctmr specifically in channeltype
% % 'nopos_iEEG'
%     % define datainfo.elec_info.orderindata.channel_type
% % get channel_type for indata_not_inctmr
%     [~,~,ind]=intersect(indata_not_inctmr,datainfo.elec_info.orderindata.label)
%     chantype_indata_not_inctmr=datainfo.elec_info.orderindata.channel_type(ind);
%     [num2cell(ind),indata_not_inctmr,chantype_indata_not_inctmr]
%     nopos_iEEG_ind=input('Which channels are iEEG?, type indices[]');
%     if ~isempty(nopos_iEEG_ind)
%        datainfo.elec_info.orderindata.channel_type(nopos_iEEG_ind)={'nopos_iEEG'};
%     end
% 
% % now channels are consistent, fix datainfo.elec_info.orderindata.channel_type
%     % set for all labels with elecpos (aka iEEG) the channel type label 'iEEG'
%   [~,ind,~]=intersect(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label);
%   datainfo.elec_info.orderindata.channel_type(ind)={'iEEG'};
%      [datainfo.elec_info.orderindata.label,datainfo.elec_info.orderindata.channel_type]
%     check=input('channel_types generally ok?')
%     
% % add datainfo.elec_info.reordered (only for iEEG channels)
%     ieeg_label=datainfo.elec_info.orderindata.label(strcmp('iEEG',datainfo.elec_info.orderindata.channel_type));
%     %add datainfo.elec_info.reordered.cell: first column name of electrode, 2 name of contact  
%     for e=1:numel(ieeg_label)
%     sel_label=ieeg_label{e};   
%     ind_=strfind(sel_label,'_');
%     datainfo.elec_info.reordered.cell_ieeg{e,1}=sel_label(1:ind_-1);
%     datainfo.elec_info.reordered.cell_ieeg{e,2}=str2num(sel_label(ind_+1:end));
%     end
%   [datainfo.elec_info.reordered.cell_ieeg,ind]=  sortrows( datainfo.elec_info.reordered.cell_ieeg);
%     
%     datainfo.elec_info.reordered.sorted_ieeg=ieeg_label(ind);
%     [datainfo.elec_info.reordered.sorted_ieeg,datainfo.elec_info.reordered.cell_ieeg]
%     check=input('are electrodes and numbers well seperated?')
%     
% save(info_file,'datainfo')
% clear datainfo ind sel_label ind_ chantype_indata_not inctmr check e ieeg_label iEEG_labels inctmr_not_indata indata_not_inctmr nopo_iEEG_ind
% end
%% plot all electrodes on mni brain

% load all datainfos and get normed electrode coordinates
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_figs='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\figure\';
mkdir(path_figs)
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
         'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
         'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
         'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};

     
     load('D:\Extinction\iEEG\scripts\additional_functions\sel_colorseries.mat')
for n=1:numel(allsubs)
sel_sub=allsubs{n};  

% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)
sub{n}=sel_sub;

% for all recorded electrodes:
% all_pos{n}=datainfo.elec_info.elec_mni.elecpos;
% all_label{n}=datainfo.elec_info.elec_mni.label;

% all clean electrodes, bipolar ref
% sel clean elec
sel_labels=datainfo.elec_info.bipolar.montage_withoutartichan2.labelnew;
[labels,ind]=intersect(datainfo.elec_info.bipolar.elec_mni.label,sel_labels,'stable');

all_pos{n}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
all_label{n}=labels';

end



views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
mesh.coordsys = 'mni';
hemispheres={'left','right'};
elec_def=[-1,1];
for h=1:numel(hemispheres)
    sel_hemi=hemispheres{h};
    sel_elec_def=elec_def(h);
    load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_',sel_hemi,'.mat')));
figure('units','pixels','position',[0 0 800 1000])
ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
    hold on
    % elec to plot
    elec_toplot.unit ='mm';
    elec_toplot.coordsys ='mni';
    
    for n=1:numel(sub)

     sel_elec=sign(all_pos{n}(:,1))==sel_elec_def;

    elec_toplot.label=all_label{n}(sel_elec,:);
    elec_toplot.elecpos=all_pos{n}(sel_elec,:);
    elec_toplot.chanpos=all_pos{n}(sel_elec,:);
    sel_color=sel_col{n};
    if numel(sel_elec)>0
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color);
    end
    end
    view([-90 20]);
    material dull;
    view(squeeze(views(h,1,:))');
    c1=camlight(0,0);
    set(c1, 'style', 'infinite');

    view(squeeze(views(h,2,:))');
    c2=camlight(0, 0);
    set(c2, 'style', 'infinite');

    view(squeeze(views(h,3,:))');
    print('-f1','-r600','-dtiff',fullfile(path_figs,strcat('bipolar',sel_hemi,'_lat.tiff'))) 

    view(squeeze(views(h,4,:))');
    print('-f1','-r600','-dtiff',fullfile(path_figs,strcat('bipolar',sel_hemi,'_med.tiff'))) 
    clear c1 c2 
cd (path_figs)
    % get movie
    azis=-180:0.5:180;
    for i=1:numel(azis)
        view([azis(i) 0]);
        F(i)=getframe(gcf);
        drawnow
    end
    writerObj=VideoWriter(strcat('video_','bipolar',sel_hemi,'.avi'));
   % writerObj.Height=1880;
    open(writerObj);
    for i=1:length(F)
        frame=F(i);
        writeVideo(writerObj,frame)
    end
    close(writerObj);
    
    close all
end
    
%% query datainfo for different regions
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

% setup definition

sel_atlas='aal';
sel_regions={'Amygdala_L','Hippocampus_L','Rectus_L','Frontal_Sup_Orb_L','Frontal_Med_Orb_L','Cingulum_Ant_L',...
            'Amygdala_R','Hippocampus_R','Rectus_R','Frontal_Sup_Orb_R','Frontal_Med_Orb_R','Cingulum_Ant_R'};
distance_region=3; % distance of electrodes from region

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
         'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
         'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
         'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};


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

%% plot a specific region and its electrodes
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_templates='D:\matlab_tools\fieldtrip-20200130\template\atlas\';
path_figs='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\figure\';

% aal region def
%sel_regions={'Amygdala_L','Hippocampus_L'};
%sel_regions={'Rectus_L','Frontal_Sup_Orb_L','Frontal_Med_Orb_L'};
%sel_regions={'Amygdala_R','Hippocampus_R','Rectus_R','Frontal_Sup_Orb_R','Frontal_Med_Orb_R'};
%sel_regions={'Amygdala_L','Hippocampus_L','Rectus_L','Frontal_Sup_Orb_L','Frontal_Med_Orb_L'};

% all_roi.hip_l={'Left-Hippocampus'};
% all_roi.hip_r={'Right-Hippocampus'};
% all_roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
% all_roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% %roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% all_roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% %roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% all_roi.amy_r={'Right-Amygdala'};
% all_roi.amy_l={'Left-Amygdala'};
% all_roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
% %roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%all_roi.amy={'Left-Amygdala','Right-Amygdala'};

% sel_regions={'Left-Amygdala','Left-Hippocampus','ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal'};
% sel_color= [1,0,0;0,1,0;0,0,1;0,0,1];%rgb value for each region

sel_regions={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital',...
    'ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole',...
    'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
sel_color= [1,1,0;1,1,0;1,1,0;...
    1,1,0;1,1,0;1,1,0;1,1,0;...
    0,1,1;0,1,1];
%     
%     



plot_elec='no';
whole_brain='yes';
distance_region=2; % distance of electrodes from region to be included in the plot
hemisphere='left';
sel_atlas='freesurferDK';
%sel_color= [1,0,0;0,1,0;0,0,1];%rgb value for each region
%sel_color= [0,0,1;0,0,1;0,0,1];%rgb value for each region

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
         'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
         'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
         'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25','c_sub26','c_sub29'};



views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
f1=figure
hold on

switch sel_atlas
    case 'aal'
        atlas_def=strcat(path_templates,'aal\ROI_MNI_V4.nii');
        atlas=ft_read_atlas(atlas_def);

    case 'freesurferDK'
       load('D:\Extinction\iEEG\data\template\atlasDK_norm2mni.mat') 
       atlas=norm_atlasDK;
       atlas.anatomylabel=atlas.aparclabel;
       atlas.coordsys='mni';
end


switch whole_brain
    case 'yes'
        switch hemisphere
            case 'left'
                load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_left.mat')));
                ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
                h=1;
            case 'right'
                load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_right.mat')));
                ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
                h=2;
        end
    otherwise
        switch hemisphere
            case 'left'
                h=1;
            case 'right'
                h=2;
        end
end

for r=1:numel(sel_regions)
    % get mesh for a region
    cfg=[];
    cfg.atlas=atlas;
    cfg.inputcoord='mni';
    cfg.roi=sel_regions{r};
    sel_mask{r}=ft_volumelookup(cfg,atlas);
    
    seg{r}=keepfields(atlas,{'dim','unit','coordsys','transform'});
    seg{r}.brain=sel_mask{r};
    
    cfg=[];
    cfg.method='iso2mesh';
    cfg.radbound=1;
    cfg.maxsurf=0;
    cfg.tissue='brain';
    cfg.numvertices=100000;
    %cfg.smooth=8;
    cfg.spmversion='spm12';
    sel_mesh{r}=ft_prepare_mesh(cfg,seg{r});
    
    ft_plot_mesh(sel_mesh{r},'facealpha',0.2,  'edgealpha',0.2,'facecolor', sel_color(r,:),'edgecolor','none');
end

% add electrode in this region
for n=1:numel(allsubs)
    sel_sub=allsubs{n};
    
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    % query datainfo
    cfg.atlas=sel_atlas;
    cfg.region= sel_regions;
    cfg.distance_region= distance_region;
    cfg.ref='bipolar';
    cfg.clean='clean';
    cfg.ana_labels=sel_atlas;
    elec_selection=mcf_elec_region_selector(cfg,datainfo);
    
    % sum elec in each region
    all_pos{n}=vertcat(elec_selection.elecpos_per_roi{:});
    sub{n}=sel_sub;
    all_label{n}=horzcat(elec_selection.labels_per_roi{:});
    
    clear elec_selection elec_region_mat
end
switch plot_elec
    case 'yes'
% elec to plot
elec_toplot.unit ='mm';
elec_toplot.coordsys ='mni';
load('D:\Extinction\iEEG\extinction_ieeg_scripts\additional_functions\sel_colorseries.mat')

for n=1:numel(sub)   
    elec_toplot.label=all_label{n};
    elec_toplot.elecpos=all_pos{n};
    elec_toplot.chanpos=all_pos{n};
    sel_color=sel_col{n};
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color,'elecsize',4);
end

    case 'no'
end
view([-90 20]);
material dull;
view(squeeze(views(h,1,:))');
c1=camlight(0,0);
set(c1, 'style', 'infinite');

view(squeeze(views(h,2,:))');
c2=camlight(0, 0);
set(c2, 'style', 'infinite');

view(squeeze(views(h,3,:))');
print(f1,'-r600','-dtiff',fullfile(path_figs,strcat(hemisphere,strcat(sel_regions{:}),'bipolar','brain',whole_brain,'_lat.tiff')))

view(squeeze(views(h,4,:))');
print(f1,'-r600','-dtiff',fullfile(path_figs,strcat(hemisphere,strcat(sel_regions{:}),'bipolar','brain',whole_brain,'_med.tiff')))
clear c1 c2



%% get labels using individualized freesurfer coordinates and normalized  coordinates for bipolar contacts (first need to construct bipolar referencencing scheme)
% 
% use individual freesurfer atlas
% use mni based atlas
% 
% search with different query ranges: 1,3,5,7,9,11 (r=query/2-0.5),
% radius:0,1,2,3,4,5
% 
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
 %allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
 %          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

%allsubs = {'p_sub08'};

%all_atlas={'freesurferDestrieux','afni','brainweb','freesurferDK','aal'}%,}
all_atlas={'freesurferDK','aal'};
query_range=[1,3,5,7];

for i=1:numel(allsubs)
% save labels in file & as excel sheet
sel_sub=allsubs{i};
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

for a=1:numel(all_atlas)
path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');
%path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
path_templates='D:\matlab_tools\fieldtrip-20200130\template\atlas\';

sel_atlas=all_atlas{a};

if ~any(strcmp(datainfo.elec_info.elec_ct_mr.label,datainfo.elec_info.elec_mni.label))
error('electrode labels in elec_ct_mr and elect_mni do not match');
end

elec_labels.labels=datainfo.elec_info.elec_ct_mr.label;
elec_labels.query_ranges=query_range;
datainfo.elec_info.ana_labels.labels=datainfo.elec_info.elec_ct_mr.label;
datainfo.elec_info.ana_labels.query_ranges=query_range;

      switch sel_atlas
        case 'freesurferDK'
            elecinfo=datainfo.elec_info.bipolar.elec_ct_mr;            
            cd(path_anat)
            load('atlasDK.mat')
            atlas=atlasDK;
            atlas.coordsys='fsaverage'; 
            inputcoord='fsaverage';
            datainfo.elec_info.bipolar.ana_labels.freesurferDK_def.file=strcat(path_anat,'atlasDK.mat');
            datainfo.elec_info.bipolar.ana_labels.freesurferDK_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.bipolar.ana_labels.freesurferDK_def.atlas_unit=atlas.unit;
            datainfo.elec_info.bipolar.ana_labels.freesurferDK_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.bipolar.ana_labels.freesurferDK_def.elec_def=elecinfo;
        case 'freesurferDestrieux'
            elecinfo=datainfo.elec_info.bipolar.elec_ct_mr;              
            cd(path_anat)
            load('atlasDest.mat')
            atlas=atlasDest;
            atlas.coordsys='fsaverage'; 
            inputcoord='fsaverage';
            datainfo.elec_info.bipolar.ana_labels.freesurferDestrieux_def.file=strcat(path_anat,'atlasDest.mat');
            datainfo.elec_info.bipolar.ana_labels.freesurferDestrieux_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.bipolar.ana_labels.freesurferDestrieux_def.atlas_unit=atlas.unit;
            datainfo.elec_info.bipolar.ana_labels.freesurferDestrieux_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.bipolar.ana_labels.freesurferDestrieux_def.elec_def=elecinfo;
        case 'aal'
             elecinfo=datainfo.elec_info.bipolar.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'aal\ROI_MNI_V4.nii');
             atlas=ft_read_atlas(atlas_file);
            inputcoord='mni';
            datainfo.elec_info.bipolar.ana_labels.aal_def.file=atlas_file;
            datainfo.elec_info.bipolar.ana_labels.aal_def.atlas_resolution=atlas.hdr.volres;
            datainfo.elec_info.bipolar.ana_labels.aal_def.atlas_unit=atlas.unit;
            datainfo.elec_info.bipolar.ana_labels.aal_def.sphereradius_mm=(query_range-1);
            datainfo.elec_info.bipolar.ana_labels.aal_def.elec_def=elecinfo;

        case 'afni'
             elecinfo=datainfo.elec_info.bipolar.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'afni\TTatlas+tlrc.HEAD');
             atlas=ft_read_atlas(atlas_file);
             inputcoord='mni';
            datainfo.elec_info.bipolar.ana_labels.afni_def.file=atlas_file;
            datainfo.elec_info.bipolar.ana_labels.afni_def.atlas_resolution=[1 1 1];
            datainfo.elec_info.bipolar.ana_labels.afni_def.atlas_unit=atlas.unit;
            datainfo.elec_info.bipolar.ana_labels.afni_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.bipolar.ana_labels.afni_def.elec_def=elecinfo;
        case 'brainweb'
             elecinfo=datainfo.elec_info.bipolar.elec_mni;
             cd(path_templates)
             atlas_file=strcat(path_templates,'brainweb\brainweb_discrete.mat');
             load(atlas_file)
             inputcoord='mni';
            datainfo.elec_info.bipolar.ana_labels.brainweb_def.file=atlas_file;
            datainfo.elec_info.bipolar.ana_labels.brainweb_def.atlas_resolution=[1 1 1];
            datainfo.elec_info.bipolar.ana_labels.brainweb_def.atlas_unit=atlas.unit;
            datainfo.elec_info.bipolar.ana_labels.brainweb_def.sphereradius_mm=(query_range-1).*0.5;
            datainfo.elec_info.bipolar.ana_labels.brainweb_def.elec_def=elecinfo;
        otherwise
    end
      
for r=1:numel(query_range)
sel_range=query_range(r);

%parfor possible
parfor e=1:numel(elecinfo.label)
sel_pos=elecinfo.elecpos(e,:);
cfg            = [];
cfg.roi        = sel_pos;
cfg.atlas      = atlas;
cfg.output     = 'label';
cfg.minqueryrange = sel_range;
cfg.maxqueryrange = sel_range;
cfg.inputcoord=inputcoord;
labels = ft_volumelookup(cfg, atlas);

if sum(labels.count)>0
tmp_label(e,r)={labels.name(labels.count>0)'};
else
tmp_label(e,r)={'no label found'};
end

end
end

eval(strcat('datainfo.elec_info.bipolar.ana_labels.',sel_atlas,'=tmp_label;'));
clear tmp_label
end

save(info_file,'datainfo')
end
