%% first automatic artifact checks


path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
mkdir(path_artifactinfo)
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20','p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};



for sub=20:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    % first load readin data
    display(strcat('selected data of: ',sel_sub))
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % select ieeg channels & artifactfree channels
    firstchannelcheck=datainfo.artifact_info.firstchannelcheck;
    all_elec=[];
    all_out=[];
    for i=1:numel(firstchannelcheck.elec)
        sel_elec=firstchannelcheck.labels{i};
        out_elec_ind=[firstchannelcheck.out_electrode{i},firstchannelcheck.artifact_electrode{i}]
        out_elec=sel_elec(out_elec_ind);
        sel_elec(out_elec_ind)=[];
        all_elec=[all_elec;sel_elec];
        all_out=[all_out;out_elec];
    end
    
    cfg=[];
    cfg.channel= all_elec;
    data=ft_preprocessing(cfg,data);
    
    % bipolar referencing
    montage=datainfo.elec_info.bipolar.montage;
    for e=1:numel(all_out)
      sel_elec=all_out{e};
        ind_old=strcmp(montage.labelold,sel_elec);
        montage.tra(:,ind_old)=[];
        montage.labelold(ind_old)=[];
    end
    labelnew_ind=sum(montage.tra~=0,2)==2;
    montage.labelnew=montage.labelnew(labelnew_ind);
    montage.tra=montage.tra(labelnew_ind,:);
    
    data_bip = ft_apply_montage(data,montage);
    
    % run artifact check on "data reference" and bipolar reference
    cfg=[];
    cfg.segmentlength=0.1;
    cfg.segmentoverlap=0.05;
    cfg.zcutoff=6;
    cfg.metric={'range','kurtosis','zvalue','var_over_one'};
    cfg.hpfilter='yes';
    cfg.hpfreq=0.1;
    cfg.iteration=2;%'no_rejected';
    auto_artifact.dataref=mcf_continuous_autoarti(cfg, data)
    auto_artifact.databip=mcf_continuous_autoarti(cfg, data_bip)
   save (strcat(path_artifactinfo,sel_sub,'_autoarti_segment',num2str(cfg.segmentlength.*1000),'ms'),'auto_artifact')
keep allsubs path_info path_artifactinfo path_preproc
end

%