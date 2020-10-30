%% get electrode positions in 3d file to overlay on mr/atlas
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

%allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%          'c_sub11','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17'};

% construct a colormap for specific freesurferDK subregions
       roi_def={'Left-Hippocampus','Right-Hippocampus','Left-Amygdala','Right-Amygdala',...
                'ctx-rh-medialorbitofrontal', 'ctx-lh-medialorbitofrontal',...
                'Left-Cerebral-White-Matter','Left-Cerebellum-White-Matter','Right-Cerebral-White-Matter','Right-Cerebellum-White-Matter'};
                
       color_def=[1 0 0; 1 0 0; 0 1 0; 0 1 0;...
                  0 0 1; 0 0 1;...
                  1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1];
       
    for i=1%6:numel(allsubs)
    sel_sub=allsubs{i};
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    path_anat=strcat(path_data,'data\freesurfer_anat\output\',sel_sub,'\fieldtrip\',sel_sub,'\');
    path_fs=strcat(path_data,'\data\freesurfer_anat\output\',sel_sub,'\freesurfer\mri\');
    path_out=strcat(path_data,'data\elec_masks\',sel_sub,'\');
    mkdir(path_out)
    
    
    cd(path_anat)
    mr_file=strcat('T1_fs.nii');
    mr = ft_read_mri(mr_file);
    mr.coordsys='acpc';
    load('atlasDK.mat')
    atlasDK.coordsys='acpc';
    
    tmp_colormap=0.5*ones(numel(atlasDK.aparclabel),3);
    for r=1:numel(roi_def)
    ind=strcmp(roi_def{r},atlasDK.aparclabel);
    tmp_colormap(ind,:)=color_def(r,:);
    end
    tmp_colormap=[0.5 0.5 0.5;tmp_colormap];
    
    keep_going=1;
    try 
    sel_elec=input('Which electrode/electrode selection do you want to plot? type like this {''A''1''} or {''A1'',''A2''}' );
    catch
    display('elec definition was wrong, try again')    
    sel_elec=input('Which electrode/electrode selection do you want to plot? type like this {''A''1''} or {''A1'',''A2''}' );

    end
    display(strcat('selected electrodes: ',sel_elec))
    
    
    
   while keep_going
    % plot one or a selected group of electrodes on a mr/ directly on atlas
   

        [~,~,ind]=intersect(sel_elec,datainfo.elec_info.elec_ct_mr.label);               
        cfg.roi =datainfo.elec_info.elec_ct_mr.elecpos(ind,:);
        cfg.sphere =repmat(2,numel(ind),1);%radius of each sphere in cm/mm dep on unit of input
        cfg.round2nearestvoxel ='yes';
        mask = ft_volumelookup(cfg, mr);
        mr.elecmask=(mask+((mask==0).*0.25)).*(atlasDK.aparc~=0);
        mr.atlas=atlasDK.aparc;
%        cfg=[];
%        cfg.method='ortho'; 
%        cfg.funparameter ='elec';
%        cfg.atlas=atlasDK;
%        ft_sourceplot(cfg,mr)
       

       
       
       cfg=[];
       cfg.method='ortho'; 
       cfg.funparameter ='atlas';
       cfg.maskparameter='elecmask';
       cfg.funcolormap=tmp_colormap;
      % cfg.opacitylim    = 'zeromax';
   %   cfg.maskstyle='colormix';
        cfg.opacitylim    = [0.0 1];
        cfg.opacitymap    = 'rampup';
       cfg.location      =datainfo.elec_info.elec_ct_mr.elecpos(ind(1),:);
       cfg.atlas=atlasDK;
       ft_sourceplot(cfg,mr)
              
%             % save mask as mr (as overlay for mricron and to load for plotting)    
%         elec_on_mr=mr;
%         elec_on_mr.anatomy=mask;
%         cfg=[];
%         cfg.parameter='anatomy';
%         cfg.filename  = strcat(path_out,sel_sub,'_electrode_',all_elec{e},'.nii')
%         ft_sourcewrite(cfg,elec_on_mr)
    
    try 
    sel_elec=input('Which electrode/electrode selection do you want to plot? type like this {''A''1''} or {''A1'',''A2''}' );
    catch
    display('elec definition was wrong, try again')    
    sel_elec=input('Which electrode/electrode selection do you want to plot? type like this {''A''1''} or {''A1'',''A2''}' );
    end
    display(strcat('selected electrodes: ',sel_elec))
    
    keep_going=~isempty(sel_elec);
    clear tmp_colormap
   end
   end