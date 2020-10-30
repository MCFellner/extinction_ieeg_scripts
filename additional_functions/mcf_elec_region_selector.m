% select electrodes in specificied sub regions

% datainfo (datainfo structure (result of electrode_localization.m )
% cfg.atlas='aal' atlas on which selected roi is based {'freesurferDestrieux','afni','brainweb'};%'freesurferDK','aal'}
% cfg.region= {'Hippocampus_R'} or {'Hippocampus_R','Hippocampus_L'}
% cfg.distance_region= maximum distance to region alloed (in mm)
% cfg.ref= 'no_reref', 'bipolar'
% cfg.clean='preclean' 'clean' % currently only for bipolar

% output elec_selection
% elec_region_mat= elc*roi matrix (1 region match)
% roixroi_mat = roixroi matrix (number of electrodes in combination of regions)
% count = number of electrodes
% labels_per_roi = labels of electrodes in roi (use for data selection)
% elecpos_per_roi = position, use for veryfication and plotting

function elec_selection=mcf_elec_region_selector(cfg,datainfo)
sel_atlas=cfg.atlas;
sel_regions=cfg.region;
sel_dist=cfg.distance_region;


switch cfg.ref
    case 'no_reref'
sel_atlas_def=getfield(datainfo.elec_info.ana_labels, strcat(sel_atlas,'_def'));
    case 'bipolar'
sel_atlas_def=getfield(datainfo.elec_info.bipolar.ana_labels, strcat(sel_atlas,'_def'));
end
% check whether labels are part of the atlas
switch sel_atlas
    case 'freesurferDK'
        try
            load(sel_atlas_def.file)
            atlas=atlasDK;
            atlas_labels=atlas.aparclabel;
        catch
        end
    case 'freesurferDestrieux'
        try
            load(sel_atlas_def.file)
            atlas=atlasDest;
            atlas_labels=atlas.aparclabel;
        catch
        end
    case 'aal'
        try
            atlas=ft_read_atlas(sel_atlas_def.file);
            atlas_labels=atlas.tissuelabel;
        catch
        end
    case 'afni'
        try
            atlas=ft_read_atlas(sel_atlas_def.file);
            atlas_labels=[atlas.brick0label;atlas.brick1label];
        catch
        end
    case 'brainweb'
        try
            load(sel_atlas_def.file)
            atlas_labels=atlas.tissuelabel;
        catch
        end
    otherwise
end

if exist ('atlas_labels')
    for r=1:numel(sel_regions)
        check= strcmp(sel_regions{r},atlas_labels);
        if sum(check)==0
            error('selected region label not part of the selected atlas')
        end
    end
else
    warning('could not load atlas, cannot check whether labels are correct')
end
% query datainfo (make a subfunction out of this)
query_ind=nearest(sel_atlas_def.sphereradius_mm,sel_dist);
if sel_dist~=sel_atlas_def.sphereradius_mm(query_ind)
    warning(strcat('selected distance not possible, using ', num2str(sel_atlas_def.sphereradius_mm(query_ind))))
end

switch cfg.ref
    case 'no_reref'
    sel_labels=getfield(datainfo.elec_info.ana_labels, sel_atlas);
    case 'bipolar'
             sel_labels=getfield(datainfo.elec_info.bipolar.ana_labels, sel_atlas);

        switch cfg.clean
        case 'clean' 
    clean_labels=datainfo.elec_info.bipolar.montage_withoutartichan2.labelnew;
    [~,sel_ind,~]=intersect(datainfo.elec_info.bipolar.elec_mni.label,clean_labels,'stable');
    sel_labels=sel_labels(sel_ind,:);    
    datainfo.elec_info.bipolar.elec_mni.label= datainfo.elec_info.bipolar.elec_mni.label(sel_ind);
    datainfo.elec_info.bipolar.elec_mni.elecpos=datainfo.elec_info.bipolar.elec_mni.elecpos(sel_ind,:);
    end

        
    
end
        
sel_labels=sel_labels(:,query_ind);
for r=1:numel(sel_regions)
    sel_region=sel_regions{r};
    for e=1:numel(sel_labels)
        elec_region_mat(e,r)=sum(strcmp(sel_region,sel_labels{e}));
    end
end

% matrix of electrodes (also showing overlap)
elec_selection.mat=elec_region_mat;
elec_selection.roixroi_mat=elec_region_mat'*elec_region_mat;
% numbers of electrodes
elec_selection.count=elec_selection.roixroi_mat(find(diag(ones(numel(sel_regions),1))));
% select label of electrodes for each region
               
for r=1:numel(sel_regions)
    switch cfg.ref
    case 'no_reref'
    elec_selection.labels_per_roi{r}=datainfo.elec_info.elec_mni.label(find(elec_region_mat(:,r)));
    elec_selection.elecpos_per_roi{r}=datainfo.elec_info.elec_mni.elecpos(find(elec_region_mat(:,r)),:);
        case 'bipolar'         
    elec_selection.labels_per_roi{r}=datainfo.elec_info.bipolar.elec_mni.label(find(elec_region_mat(:,r)));
    elec_selection.elecpos_per_roi{r}=datainfo.elec_info.bipolar.elec_mni.elecpos(find(elec_region_mat(:,r)),:);
    end
end

elec_selection.cfg=cfg;