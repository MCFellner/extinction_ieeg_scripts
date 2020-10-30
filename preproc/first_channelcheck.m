path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20','p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% show bipolar reference below in same databrowser
% print labels of electrodes in command window for easy checking
% also select possible wm contact for referencing

for sub=5%:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    % first load readin data
    display(strcat('selected data of: ',sel_sub))
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % select iEEG channels
    cfg=[];
    cfg.channel= datainfo.elec_info.orderindata.label(strcmp(datainfo.elec_info.orderindata.channel_type,'iEEG'));
    %     cfg.reref         =  'yes';
    %     cfg.refchannel    =  cfg.channel;
    %     cfg.refmethod     = 'avg';
    data=ft_preprocessing(cfg,data);
    
    elecs=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));
    
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        sel_channel=datainfo.elec_info.reordered.sorted_ieeg(strcmp(datainfo.elec_info.reordered.cell_ieeg(:,1),sel_elec));
        
        cfg=[];
        cfg.channel=sel_channel;
        sel_data=ft_preprocessing(cfg,data);
        
        % fix channel order ( ft_preprocessing keeps old order)
        [~,~,ind]=intersect(sel_channel,sel_data.label,'stable');
        sel_data.label=sel_data.label(ind);
        sel_data.trial{1,1}=sel_data.trial{1,1}(ind,:);
        
        % bipolar referencing
        % select subset of montage
        [~,~,ind]=intersect(sel_channel,datainfo.elec_info.bipolar.montage.labelold,'stable');
        sel_montage.labelold=datainfo.elec_info.bipolar.montage.labelold(ind);
        sel_montage.tra=datainfo.elec_info.bipolar.montage.tra(:,ind);
        sel_montage.labelnew=datainfo.elec_info.bipolar.montage.labelnew(find(sum(sel_montage.tra~=0,2)))';
        sel_montage.tra=sel_montage.tra(find(sum(sel_montage.tra~=0,2)),:);
        
        sel_data_bip = ft_apply_montage(sel_data,sel_montage);
        clear sel_montage
        data_toplot = ft_appenddata(cfg, sel_data, sel_data_bip)
        
        [~,~,sel_ana]=intersect(sel_channel,datainfo.elec_info.ana_labels.labels,'stable');
        ana_labels=[sel_channel,datainfo.elec_info.ana_labels.freesurferDK(sel_ana,1),...
            datainfo.elec_info.ana_labels.freesurferDestrieux(sel_ana,1),...
            datainfo.elec_info.ana_labels.aal(sel_ana,1),...
            datainfo.elec_info.ana_labels.afni(sel_ana,1),...
            datainfo.elec_info.ana_labels.brainweb(sel_ana,1)]
        possible_WMref=datainfo.elec_info.whitematter.possible_WM_reference{e}
        
        cfg            = [];
        cfg.viewmode   = 'vertical';
        cfg.ylim   =  'maxmin';
        cfg.blocksize  =30;
        cfg.continuous              = 'yes';
        ft_databrowser(cfg, data_toplot);
        
        
        % put breakpoint here for easier datacheck
        electrode_check.elec{e}=sel_elec;
        electrode_check.labels{e}=sel_channel;
        electrode_check.all_ok{e}=input(strcat('Is electrode ',sel_elec,' generally ok? write [1] if ok, [0] if not'));
        electrode_check.flat_electrode{e}=input(strcat('Is there a flat channel on electrode ',sel_elec,'? write electrode like that: [10] or [] if none'));
        electrode_check.out_electrode{e}=input(strcat('Is there a out of brain channel on electrode ',sel_elec,'? write electrode like that: [10] or [] if none'));
        electrode_check.artifact_electrode{e}=input(strcat('Is there a artifact/noise channel on electrode ',sel_elec,'? write electrode like that: [10] or [] if none'));
        electrode_check.WM_ref{e}=input(strcat('Which WM electrode to select as possible WM ref ',[possible_WMref{:}],'? write electrode like that: [10] or [] if none'));
        close all
        delete(gcf)
    end
    datainfo.artifact_info.firstchannelcheck=electrode_check;
    save(info_file,'datainfo');
    keep path_info n allsubs path_preproc
end