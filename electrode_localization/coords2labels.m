%% get labels using individualized freesurfer coordinates and normalized  coordinates

% use individual freesurfer atlas
% use mni based atlas
% 
% search with different query ranges: 1,3,5,7,9,11 (r=query/2-0.5),
% radius:0,1,2,3,4,5

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
%allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

allsubs = {'p_sub08'};

all_atlas={'freesurferDestrieux','afni','brainweb','freesurferDK','aal'}%,}
query_range=[1,3,5,7,9,11];

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
