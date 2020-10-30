
%cfg
%cfg.elec_pos cell with elec for group col_def/trans_def in one cell
%cfg.elec_label cell with elec for group col_def/trans_def in one cell

%cfg.col_def groupx3 rgb definition
%cfg.trans_def group x 1 transperancy definition (1 opaque, 0 trans)

%cfg.views views to be plotted 
%cfg.hemispheres hemisphere to be plotted
%cfg.hemisphere_surf file with template

%cfg.path_fig path where figure will be saved

function mcf_ieegelectrodeplotter(cfg)

% elec to plot
elec_toplot.unit ='mm';
elec_toplot.coordsys ='mni';
mesh.coordsys = 'mni';

for h=1:numel(cfg.hemispheres)
    sel_hemi=cfg.hemispheres{h};
     % sign for left/right hemispere (assumes ras organization)
    switch sel_hemi
        case 'left'
        sel_elec_def=-1;
        case 'right'
        sel_elec_def=1;
    end
    
    load(cfg.hemisphere_surf{h});
    figure
    ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
    hold on
   
    for t=1:numel(cfg.elec_pos)
    % all elecs in selected hemisphere
    sel_elec=sign(cfg.elec_pos{t}(:,1))==sel_elec_def;
    
    sel_trans=cfg.trans_def(t);
    elec_toplot.label=cfg.elec_label{t}(sel_elec,:);
    elec_toplot.elecpos=cfg.elec_pos{t}(sel_elec,:);
    elec_toplot.chanpos=cfg.elec_pos{t}(sel_elec,:);
    sel_color=cfg.col_def(t,:);
    
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color,'facealpha',sel_trans);
    end
    view([-90 20]);
    material dull;
    view(squeeze(cfg.views(h,1,:))');
    c1=camlight(0,0);
    set(c1, 'style', 'infinite');

    view(squeeze(cfg.views(h,2,:))');
    c2=camlight(0, 0);
    set(c2, 'style', 'infinite');

    view(squeeze(cfg.views(h,3,:))');
    print('-f1','-r600','-dtiff',fullfile(cfg.path_fig,strcat(sel_hemi,'_lat.tiff'))) 

    view(squeeze(cfg.views(h,4,:))');
    print('-f1','-r600','-dtiff',fullfile(cfg.path_fig,strcat(sel_hemi,'_med.tiff'))) 
    clear c1 c2 
    close all
end
