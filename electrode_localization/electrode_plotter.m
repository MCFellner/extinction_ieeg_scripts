% %% plot all electrodes on mni brain
%
%
% % load all datainfos and get normed electrode coordinates
% path_data='D:\Extinction\iEEG\';
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_figs='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\figure\';
% mkdir(path_figs)
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%          'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};
%
%
%      load('D:\Extinction\iEEG\scripts\additional_functions\sel_colorseries.mat')
% for n=1:numel(allsubs)
% sel_sub=allsubs{n};
%
% % electrodeinfo
% info_file=strcat(path_info,sel_sub,'_datainfo');
% load(info_file)
%
% all_pos{n}=datainfo.elec_info.elec_mni.elecpos;
% sub{n}=sel_sub;
% all_label{n}=datainfo.elec_info.elec_mni.label;
% end
%
%
%
% views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
% views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
% mesh.coordsys = 'mni';
% hemispheres={'left','right'};
% elec_def=[-1,1];
% for h=1:numel(hemispheres)
%     sel_hemi=hemispheres{h};
%     sel_elec_def=elec_def(h);
%     load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_',sel_hemi,'.mat')));
%     figure
%     ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
%     hold on
%     % elec to plot
%     elec_toplot.unit ='mm';
%     elec_toplot.coordsys ='mni';
%
%     for n=1:numel(sub)
%
%      sel_elec=sign(all_pos{n}(:,1))==sel_elec_def;
%
%     elec_toplot.label=all_label{n}(sel_elec,:);
%     elec_toplot.elecpos=all_pos{n}(sel_elec,:);
%     elec_toplot.chanpos=all_pos{n}(sel_elec,:);
%     sel_color=sel_col{n};
%
%
%     ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color);
%     end
%     view([-90 20]);
%     material dull;
%     view(squeeze(views(h,1,:))');
%     c1=camlight(0,0);
%     set(c1, 'style', 'infinite');
%
%     view(squeeze(views(h,2,:))');
%     c2=camlight(0, 0);
%     set(c2, 'style', 'infinite');
%
%     view(squeeze(views(h,3,:))');
%     print('-f1','-r600','-dtiff',fullfile(path_figs,strcat(sel_hemi,'_lat.tiff')))
%
%     view(squeeze(views(h,4,:))');
%     print('-f1','-r600','-dtiff',fullfile(path_figs,strcat(sel_hemi,'_med.tiff')))
%     clear c1 c2
% cd (path_figs)
%     % get movie
%     azis=-180:2:180;
%     for i=1:numel(azis)
%         view([azis(i) 0]);
%         F(i)=getframe(gcf);
%         drawnow
%     end
%     writerObj=VideoWriter(strcat('video_',sel_hemi,'.avi'));
%     writerObj.FrameRate=10;
%
%     open(writerObj);
%     for i=1:length(F)
%         frame=F(i);
%         writeVideo(writerObj,frame)
%     end
%     close(writerObj);
%
%     close all
% end
%
%% plot a specific region and its electrodes
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_templates='D:\matlab_tools\fieldtrip-20200130\template\atlas\';
path_figs='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\figure\';
sel_regions={'Amygdala_L','Hippocampus_L'};
%sel_regions={'Rectus_L','Frontal_Sup_Orb_L','Frontal_Med_Orb_L'};
%sel_regions={'Amygdala_R','Hippocampus_R','Rectus_R','Frontal_Sup_Orb_R','Frontal_Med_Orb_R'};
whole_brain='no';
distance_region=2; % distance of electrodes from region to be included in the plot
hemisphere='right';
sel_atlas='aal';
%sel_color= [1,0,0;0,1,0;0,0,1];%rgb value for each region
%sel_color= [0,0,1;0,0,1;0,0,1];%rgb value for each region
sel_color= [1,0,0;0,1,0;0,0,1;0,0,1;0,0,1]

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};


views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
f1=figure
hold on

switch sel_atlas
    case 'aal'
        atlas_def=strcat(path_templates,'aal\ROI_MNI_V4.nii');
    otherwise
end
atlas=ft_read_atlas(atlas_def);


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
    cfg.radbound=2;
    cfg.maxsurf=0;
    cfg.tissue='brain';
    cfg.numvertices=5000;
    cfg.smooth=3;
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
    elec_selection=mcf_elec_region_selector(cfg,datainfo);
    
    % sum elec in each region
    all_pos{n}=vertcat(elec_selection.elecpos_per_roi{:});
    sub{n}=sel_sub;
    all_label{n}=vertcat(elec_selection.labels_per_roi{:});
    
    clear elec_selection elec_region_mat
end

% elec to plot
elec_toplot.unit ='mm';
elec_toplot.coordsys ='mni';
load('D:\Extinction\iEEG\scripts\additional_functions\sel_colorseries.mat')

for n=1:numel(sub)   
    elec_toplot.label=all_label{n};
    elec_toplot.elecpos=all_pos{n};
    elec_toplot.chanpos=all_pos{n};
    sel_color=sel_col{n};
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color,'elecsize',2);
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
print(f1,'-r600','-dtiff',fullfile(path_figs,strcat(hemisphere,strcat(sel_regions{:}),'brain',whole_brain,'_lat.tiff')))

view(squeeze(views(h,4,:))');
print(f1,'-r600','-dtiff',fullfile(path_figs,strcat(hemisphere,strcat(sel_regions{:}),'brain',whole_brain,'_med.tiff')))
clear c1 c2