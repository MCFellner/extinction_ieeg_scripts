%% get electrode positions in 3d file to overlay on mr/atlas
path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

%allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%          'c_sub11','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17'};

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
    
%     keep_going=1;
%     
%     sel_elec=input('Which electrode/electrode selection do you want to plot? type like this {''A''1''} or {''A1'',''A2''}' );
%     display(strcat('selected electrodes: ',sel_elec))
%     
%    while keep_going
   
    %%%%%%%%%FIX
    % generate mask for every electrode ('A')
    all_elec=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));
    
    
    for e=1:numel(all_elec)
        sel_channel=datainfo.elec_info.reordered.sorted_ieeg(strcmp(datainfo.elec_info.reordered.cell_ieeg(:,1),all_elec{e}));
        [~,~,ind]=intersect(sel_channel,datainfo.elec_info.elec_ct_mr.label)
        
        
        cfg.roi                 = datainfo.elec_info.elec_ct_mr.elecpos(ind,:);
        cfg.sphere              = repmat(2,numel(ind),1);%radius of each sphere in cm/mm dep on unit of input
        cfg.round2nearestvoxel ='yes'
        mask = ft_volumelookup(cfg, mr);
        
%             % save mask as mr (as overlay for mricron and to load for plotting)    
        elec_on_mr=mr;
        elec_on_mr.anatomy=mask;
        cfg=[];
        cfg.parameter='anatomy';
        cfg.filename  = strcat(path_out,sel_sub,'_electrode_',all_elec{e},'.nii')
        ft_sourcewrite(cfg,elec_on_mr)
    end
    
    
   
end
