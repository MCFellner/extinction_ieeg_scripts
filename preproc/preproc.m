addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% iEEG preprocessing

% - montages for rereferencing
% - first channel check
% - filters
% - some automatic checks

%% build montages for rereferencing
%check data avg reference, bipolar, or nearest white matter
% average reference: use ft_preprocessing
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};
allsubs = {'c_sub20'};
% load datainfo for each subject
for n=1:numel(allsubs)
    sel_sub=allsubs{n};
    
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    % datainfo.elec_info.reordered.cell_ieeg: defines every electrode
    % get for every electrode a rereference scheme
    elecs=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));
    
    %%%%%bipolar
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        sel_chan=datainfo.elec_info.reordered.sorted_ieeg(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1)));
        sel_num=[datainfo.elec_info.reordered.cell_ieeg{(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1))),2}];
        %check for missing number
        if any(diff(sel_num)~=1)
            bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
            bip_elec{e}(diff(sel_num)~=1,:)=[];
        else
            bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
        end
        bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
    end
    
    bip_reref=vertcat(bip_elec{:});
    bip_labelold=datainfo.elec_info.reordered.sorted_ieeg;
    
    check=setdiff(bip_labelold,datainfo.elec_info.orderindata.label);
    if ~isempty(check)
        error('labels in elec_info.orderindata and reordered do not match')
    end
    
    for chan=1:size(bip_reref,1)
        % get new name
        bip_labelnew{chan}=strcat(bip_reref{chan,1},'-',bip_reref{chan,2});
        % define tra mat (1,-1)
        tra(chan,:)=strcmp(bip_reref{chan,1},bip_labelold)-strcmp(bip_reref{chan,2},bip_labelold);
        % recalculate elec pos
        elec1_ind=strcmp(bip_reref{chan,1},datainfo.elec_info.elec_ct_mr.label);
        elec2_ind=strcmp(bip_reref{chan,2},datainfo.elec_info.elec_ct_mr.label);
        bip_pos_ct_mr(chan,:)=(datainfo.elec_info.elec_ct_mr.elecpos(elec1_ind,:)+datainfo.elec_info.elec_ct_mr.elecpos(elec2_ind,:)).*0.5;
        
        elec1_ind=strcmp(bip_reref{chan,1},datainfo.elec_info.elec_mni.label);
        elec2_ind=strcmp(bip_reref{chan,2},datainfo.elec_info.elec_mni.label);
        bip_pos_mni(chan,:)=(datainfo.elec_info.elec_mni.elecpos(elec1_ind,:)+datainfo.elec_info.elec_mni.elecpos(elec2_ind,:)).*0.5;
    end
    
    % add new info to datainfo
    datainfo.elec_info.bipolar.elec_ct_mr.label=bip_labelnew;
    datainfo.elec_info.bipolar.elec_ct_mr.chanpos=bip_pos_ct_mr;
    datainfo.elec_info.bipolar.elec_ct_mr.elecpos=bip_pos_ct_mr;
    datainfo.elec_info.bipolar.elec_ct_mr.unit='mm';
    datainfo.elec_info.bipolar.elec_ct_mr.coordsys='acpc';
    datainfo.elec_info.bipolar.elec_ct_mr.tra=diag(ones(1,numel(bip_labelnew)));
    
    datainfo.elec_info.bipolar.elec_mni.label=bip_labelnew;
    datainfo.elec_info.bipolar.elec_mni.chanpos=bip_pos_mni;
    datainfo.elec_info.bipolar.elec_mni.elecpos=bip_pos_mni;
    datainfo.elec_info.bipolar.elec_mni.unit='mm';
    datainfo.elec_info.bipolar.elec_mni.coordsys='mni';
    datainfo.elec_info.bipolar.elec_mni.tra=diag(ones(1,numel(bip_labelnew)));
    
    % build montage for each way of rereferncing
    %    montage.tra      = MxN matrix
    %    montage.labelold = Nx1 cell-array
    %    montage.labelnew = Mx1 cell-array
    
    datainfo.elec_info.bipolar.montage.tra=tra;
    datainfo.elec_info.bipolar.montage.labelold=bip_labelold;
    datainfo.elec_info.bipolar.montage.labelnew=bip_labelnew;
    
    % check bipolar montage
    %   [data_bip]    = ft_apply_montage(data,datainfo.elec_info.bipolar.montage)
    
    %%%%%% white matter contact
    % load data
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % check for each electrode for white matter contacts (in ct_mr)
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        datainfo.elec_info.whitematter.electrodes=elecs;
        sel_chan=datainfo.elec_info.reordered.sorted_ieeg(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1)));
        datainfo.elec_info.whitematter.sel_chan_per_electrodes=sel_chan;
        
        for chan=1:numel(sel_chan)
            % check for white matter label 'White-Matter'
            all_labels=datainfo.elec_info.ana_labels.freesurferDK(strcmp(sel_chan{chan},datainfo.elec_info.ana_labels.labels),:);
            for l=1:numel(all_labels)
                sel_labels=all_labels{l};
                % get numbers of found regions for each
                num_labels(chan,l)=numel(sel_labels);
                check=(strfind(sel_labels, 'White-Matter'));
                wm_labels(chan,l)=~isempty([check{:}]);
            end
        end
        relation_WM2other_freesurferDK=wm_labels./num_labels;
        datainfo.elec_info.whitematter.relation_WM2other_freesurferDK{e}=relation_WM2other_freesurferDK;
        WM_score_freesurferDK=cumsum(relation_WM2other_freesurferDK==1,2,'reverse');
        datainfo.elec_info.whitematter.WM_score_freesurferDK{e}=WM_score_freesurferDK;
        datainfo.elec_info.whitematter.possible_WM_reference{e}=sel_chan(max(WM_score_freesurferDK(:,1))==WM_score_freesurferDK(:,1));
        clear num_labels wm_labels
        [~,ind]=intersect(sel_chan,data.label,'stable')
        %check for matching data
        if numel(ind)~=numel(sel_chan)
            error('datainfo channels and data channels do not match')
        end
        sel_data=data.trial{1,1}(ind,:);
        datainfo.elec_info.whitematter.corr_mat{e}=corr(sel_data');
        datainfo.elec_info.whitematter.cov_mat{e}=cov(sel_data');
        datainfo.elec_info.whitematter.var_mat{e}=var(sel_data');
        datainfo.elec_info.whitematter.mean_mat{e}=mean(sel_data');
    end
    save(info_file,'datainfo')
    keep n allsubs path_info path_preproc
end

%% check channels in raw data
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



%% browse & autoartifact
%% autoartifact check and browse

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
mkdir(path_artifactinfo)
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

sel_sub='p_sub08';

sel_ref='bip'; % or 'org'
all_metric={ 'range','kurtosis'};


% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)
auto_arti_file=strcat(path_artifactinfo,sel_sub,'_autoarti_segment100ms')
load(auto_arti_file)
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

switch sel_ref
    case 'org'
        sel_auto_artifact=(auto_artifact.dataref);
        ref='org';
    case 'bip'
        sel_auto_artifact=(auto_artifact.databip);
        ref='bip';
        data=data_bip;
end
% evaluate detected artifacts
for m=1:numel(all_metric)
    sel_metric=all_metric{m};
    detected_artis=getfield( sel_auto_artifact,sel_metric);
    artis=detected_artis.arti_pertrial_chan;
    % relative number of artis
    rel_chan= sum(artis,2)./size(artis,2);
    rel_trial=sum(artis,1)./size(artis,1);
    figure
    subplot(3,2,3)
    imagesc(artis)
    xlabel('segments')
    ylabel('channel')
    title(strcat(ref,': ',sel_metric))
    subplot(3,2,1)
    plot(rel_trial)
    subplot(3,2,4)
    plot(rel_chan,1:size(artis,1))
    subplot(3,2,5)
    hist(rel_trial)
    title('distribution trial')
    subplot(3,2,6)
    hist(rel_chan)
    title('distribution chan')
    % relativ for electrodes
    % get measures for each electrode
    artis_all(m,:,:)=artis;
    arti_chan_all(m,:,:)=detected_artis.arti_z_across_chan;
    
end
artis=squeeze(sum(artis_all))>0;
artis_chan=squeeze(sum(arti_chan_all))>0;
% relative number of artis
rel_chan= sum(artis,2)./size(artis,2);
rel_trial=sum(artis,1)./size(artis,1);
figure
subplot(3,2,3)
imagesc(artis)
xlabel('segments')
ylabel('channel')
title(strcat(ref,': combined'))
subplot(3,2,1)
plot(rel_trial)
subplot(3,2,4)
plot(rel_chan,1:size(artis,1))
subplot(3,2,5)
hist(rel_trial)
title('distribution trial')
subplot(3,2,6)
hist(rel_chan)
title('distribution chan')

%sugest channels to remove
z_chan=zscore(sum(artis_chan,2));
rel_arti_chan=sum(artis_chan,2)./size(artis_chan,2);
reject_channel=[sel_auto_artifact.labels(z_chan>3)',num2cell(rel_arti_chan(z_chan>3))]

% check for dead channels/low variance
z_chan_novar=nanzscore(sum(sel_auto_artifact.var_over_one.arti_z_across_chan,2));
rel_novar_chan=sum(sel_auto_artifact.var_over_one.arti_z_across_chan,2)./size(sel_auto_artifact.var_over_one.arti_z_across_chan,2);
dead_channel=[sel_auto_artifact.labels(z_chan_novar>3)',num2cell(rel_novar_chan(z_chan_novar>3))]

% select threshold to mark segment as artifact
relchan_thres=0.1;
rel_trial= sum(artis)./size(artis,1);
sel_artis=rel_trial>relchan_thres;
arti_sp=sel_auto_artifact.samples(sel_artis,1:2);

% show info regarding channels in a figure
figure
subplot(1,2,1)
for e=1:numel(datainfo.elec_info.ana_labels.labels)
    text(0,numel(datainfo.elec_info.ana_labels.labels)-e,strcat(datainfo.elec_info.ana_labels.labels{e},': ',datainfo.elec_info.ana_labels.freesurferDK{e,1}))
    ylim([-2 numel(datainfo.elec_info.ana_labels.labels)])
end
axis off
subplot(1,2,2)
axis off
text(0,10,strcat('artifact elec:',[reject_channel{:}]))
hold on
text(0,5,strcat('no var elec:',[dead_channel{:}]))
ylim([0 12])

% check for channels to exclude & mark undetected artifact & remove wrong
% artifacts
cfg=[];
cfg.viewmode = 'vertical';
cfg.ylim =  'maxmin';
cfg.blocksize =30;
cfg.continuous  = 'yes';
cfg.artfctdef.xxx.artifact  =arti_sp;
artif_check=ft_databrowser(cfg,data_bip)

% add breakpoint for easier channel selection
sel_spikechan=input('spike channels? write as [{''}]')
sel_artichan=input('arti channels? write as [{''}]')
sel_arti=artif_check.artfctdef.xxx;

browsercheck.spikechan=sel_spikechan;
browsercheck.artichan=sel_artichan;
browsercheck.arti_sp=sel_arti;

datainfo.artifact_info=setfield(datainfo.artifact_info,strcat('browsercheck_',sel_ref),browsercheck);
datainfo.artifact_info.autoartifact=sel_auto_artifact;

save(info_file,'datainfo')
clear
close all
delete(gcf)

%% check trial segments

%% add additional info to datainfo: which channels to keep/reject and artifactvector

% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
% mkdir(path_artifactinfo)
% path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% for sub=1:length(allsubs)
%     sel_sub=allsubs{sub};
%     % electrodeinfo
%     info_file=strcat(path_info,sel_sub,'_datainfo');
%     load(info_file)
%     load(strcat(path_preproc,sel_sub,'_data.mat'))
%
%       % select ieeg channels & artifactfree channels
%       % remove electrode first channel check
%     firstchannelcheck=datainfo.artifact_info.firstchannelcheck;
%     allfirst_elec=[];
%     allfirst_out=[];
%     for i=1:numel(firstchannelcheck.elec)
%         sel_elec=firstchannelcheck.labels{i};
%         out_elec_ind=[firstchannelcheck.out_electrode{i},firstchannelcheck.artifact_electrode{i}]
%         out_elec=sel_elec(out_elec_ind);
%         sel_elec(out_elec_ind)=[];
%         allfirst_elec=[allfirst_elec;sel_elec];
%        allfirst_out=[allfirst_out;out_elec];
%     end
%
%     datainfo.artifact_info.firstchannelcheck.eleclabel_in=allfirst_elec;
%     datainfo.artifact_info.firstchannelcheck.eleclabel_out=allfirst_out;
%
%     % remove electrode identified as artifact in second round
%     second_round_arti=datainfo.artifact_info.browsercheck_bip;
%
%     % check for rejected electrodes
%     allsecond_out=[second_round_arti.spikechan,second_round_arti.artichan];
%     % check for correct spelling
%         % channels needs to be either part of bip or normal cha
%         labelold=datainfo.elec_info.bipolar.montage.labelold;
%         labelnew=datainfo.elec_info.bipolar.montage.labelnew;
%         no_bip=setdiff(allsecond_out,labelnew);
%         no_old=setdiff(no_bip,labelold)
%         if ~isempty(no_old)
%            input('check names of artifact electrodes again')
%            % save(info_file,'datainfo')
%         end
%     % search for matches with channel in the bip labels
%    all_out=[allfirst_out',allsecond_out]
%
%     out_elec=[];
%     allout_elec=[];
%     for e=1:numel(all_out)
%     sel_elec=all_out{e};
%     if str2double(sel_elec(end))==1 & isnan(str2double(sel_elec(end-1))) % exception for first electrode
%     out_elec=labelnew(cellfun(@isempty,(strfind(labelnew,(strcat(sel_elec,'-')))))==0);
%     else
%     out_elec=labelnew(cellfun(@isempty,(strfind(labelnew,sel_elec)))==0);
%     end
%     allout_elec=[allout_elec,out_elec];
%     end
%     allin_elec=setdiff(labelnew,allout_elec);
%
%     datainfo.artifact_info.browsercheck_bip.eleclabel_in=allin_elec;
%     datainfo.artifact_info.browsercheck_bip.eleclabel_out=allout_elec;
%
%
%     % remove rejected electrodes from bipolar referencing
%     montage=datainfo.elec_info.bipolar.montage;
%     for e=1:numel(allout_elec)
%       sel_elec=allout_elec{e};
%         ind_new=strcmp(montage.labelnew,sel_elec);
%         montage.tra(ind_new,:)=[];
%         montage.labelnew(ind_new)=[];
%     end
%     % sanity check
%     setdiff(montage.labelnew,allin_elec)
%
%     labelold_ind=sum(montage.tra~=0,1)~=0;
%     montage.labelold=montage.labelold(labelold_ind);
%     montage.tra=montage.tra(:,labelold_ind);
%
%     datainfo.elec_info.bipolar.montage_withoutartichan=montage;
%
%     % construct artif vec for easier trial selection
%     arti_vec=zeros(1,numel(data.time{1}));
%    % set arti_vec 1 for all epochs with marked artifact
%    marked_arti=datainfo.artifact_info.browsercheck_bip.arti_sp.artifact;
%    for a=1:size(marked_arti,1)
%        arti_vec(marked_arti(a,1):marked_arti(a,2))=1;
%    end
% datainfo.artifact_info.browsercheck_bip.artifact_vec=arti_vec;
% save(info_file,'datainfo')
% keep path_info  path_preproc allsubs n
% end

%% check trial segments

% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
% mkdir(path_artifactinfo)
% path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% for sub=27:length(allsubs)
%     sel_sub=allsubs{sub};
%     % electrodeinfo
%     info_file=strcat(path_info,sel_sub,'_datainfo');
%     load(info_file)
%     load(strcat(path_preproc,sel_sub,'_data.mat'))
%
%     % apply bipolar montage
%     montage=datainfo.elec_info.bipolar.montage_withoutartichan;
%     data = ft_apply_montage(data,montage);
%
%     % add arti_vec as additional channel (for easy artifact removal)
%     data.trial{1}(end+1,:)=datainfo.artifact_info.browsercheck_bip.artifact_vec;
%     data.label(end+1)={'artifact'};
%
%     trlinfo=datainfo.trialinfo;
% % segment data in different trial parts
% % item window: -1 to 4
% pre_item=1;
% post_item=4;
% trl_item(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_item);
% trl_item(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
% trl_item(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_item.*data.fsample);
% datainfo.trigger.trl_item=trl_item;
%
% pre_us=3;
% post_us=2;
% us_onset=datainfo.trigger.trigger_sp'+round(((trlinfo(:,13)-trlinfo(:,11))./10000).*data.fsample);
% trl_us(:,1)=us_onset-(data.fsample.*pre_us);
% trl_us(:,2)=us_onset+(data.fsample.*post_us);
% trl_us(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_us.*data.fsample);
% datainfo.trigger.trl_us=trl_us;
%
% cfg=[];
% cfg.trl=trl_item;
% data_item=ft_redefinetrial(cfg,data);
% cfg=[];
% cfg.keeptrials='yes';
% erp_item=ft_timelockanalysis(cfg,data_item)
% artifree_ind=squeeze(sum(erp_item.trial(:,end,:),3))==0;
% clear erp_item
% cfg=[];
% cfg.trials=find(artifree_ind);
% data_item=ft_preprocessing(cfg,data_item)
%
%
% cfg=[];
% cfg.trl=trl_us;
% data_us=ft_redefinetrial(cfg,data);
% cfg=[];
% cfg.keeptrials='yes';
% erp_us=ft_timelockanalysis(cfg,data_us)
% artifree_ind=squeeze(sum(erp_us.trial(:,end,:),3))==0;
% clear erp_us
% cfg=[];
% cfg.trials=find(artifree_ind);
% data_us=ft_preprocessing(cfg,data_us)
%
%
% % run reject visual-summary as added check
% cfg=[];
% cfg.method='summary';
% cfg.keeptrial ='nan';
% cfg.keepchannel='no';
% cfg.channel={'all','-artifact'};
% summarycheck=ft_rejectvisual(cfg,data_item)
%
% % add new rejected trials to artifacts
% arti_vec=datainfo.artifact_info.browsercheck_bip.artifact_vec;
% for a=1:size(summarycheck.cfg.artfctdef.summary.artifact,1)
%     sel_arti=summarycheck.cfg.artfctdef.summary.artifact(a,:);
% arti_vec(sel_arti(1):sel_arti(2))=1;
% end
% % add rejected channels to datainfo
% allout_elec=setdiff(data_item.label(1:end-1),summarycheck.label);
% datainfo.artifact_info.rejectvisual_bip.elecsin=data_item.label(1:end-1);
% datainfo.artifact_info.rejectvisual_bip.elecsout=allout_elec;
%
% % update montage
% montage=datainfo.elec_info.bipolar.montage_withoutartichan;
%     for e=1:numel(allout_elec)
%        sel_elec=allout_elec{e};
%          ind_new=strcmp(montage.labelnew,sel_elec);
%          montage.tra(ind_new,:)=[];
%          montage.labelnew(ind_new)=[];
%     end
%
%     labelold_ind=sum(montage.tra~=0,1)~=0;
%     montage.labelold=montage.labelold(labelold_ind);
%     montage.tra=montage.tra(:,labelold_ind);
%     datainfo.elec_info.bipolar.montage_withoutartichan2=montage;
%
%
% % run reject visual-summary as added check
% cfg=[];
% cfg.method='summary';
% cfg.keeptrial ='nan';
% cfg.keepchannel='nan';
% cfg.channel={'all','-artifact'};
% summarycheck=ft_rejectvisual(cfg,data_us)
%
% % add new rejected trials to artifacts
% for a=1:size(summarycheck.cfg.artfctdef.summary.artifact,1)
%     sel_arti=summarycheck.cfg.artfctdef.summary.artifact(a,:);
% arti_vec(sel_arti(1):sel_arti(2))=1;
% end
%
% datainfo.artifact_info.browsercheck_bip.artifact_vec_rejectvisual=arti_vec;
%
% save(info_file,'datainfo')
% end

%% data checker
% path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
% path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
%
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
%
% for sub=27:length(allsubs)
%     sel_sub=allsubs{sub};
%     % electrodeinfo
%     info_file=strcat(path_info,sel_sub,'_datainfo');
%     load(info_file)
%     load(strcat(path_preproc,sel_sub,'_data.mat'))
%
%     % apply bipolar montage
%     montage=datainfo.elec_info.bipolar.montage_withoutartichan;
%     data = ft_apply_montage(data,montage);
%
%     % add arti_vec as additional channel (for easy artifact removal)
%     data.trial{1}(end+1,:)=datainfo.artifact_info.browsercheck_bip.artifact_vec;
%     data.label(end+1)={'artifact'};
%
%     trlinfo=datainfo.trialinfo;
% % segment data in different trial parts
% % item window: -1 to 4
% pre_item=1;
% post_item=4;
% trl_item(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_item);
% trl_item(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
% trl_item(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_item.*data.fsample);
% datainfo.trigger.trl_item=trl_item;
%
% pre_us=1;
% post_us=4;
% us_onset=datainfo.trigger.trigger_sp'+round(((trlinfo(:,13)-trlinfo(:,11))./10000).*data.fsample);
% trl_us(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_us);
% trl_us(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_us);
% trl_us(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_us.*data.fsample);
% datainfo.trigger.trl_us=trl_us;
% cfg=[];
% cfg.trl=trl_item;
% data_item=ft_redefinetrial(cfg,data);
% cfg=[];
% cfg.keeptrials='yes';
% erp_item=ft_timelockanalysis(cfg,data_item)
% artifree_ind=squeeze(sum(erp_item.trial(:,end,:),3))==0;
% clear erp_item
% cfg=[];
% cfg.trials=find(artifree_ind);
% data_item=ft_preprocessing(cfg,data_item)
%
%
% cfg=[];
% cfg.viewmode='vertical';
% ft_databrowser(cfg,data_item)
% close all
% delete(gcf)
% clear
% end
%% count trial numbers


conditions={'A_csplusplus','B_csplusplus','C_csplusplus',...
    'A_csplusminus','B_csplusminus','C_csplusminus',...
    'A_csminusminus','B_csminusminus','C_csminusminus'};
cond_def=[1,1;2,1;3,1;...
    1,2;2,2;3,2;...
    1,3;2,3;3,3]

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);
    
    % add arti_vec as additional channel (for easy artifact removal)
    data.trial{1}(end+1,:)=datainfo.artifact_info.browsercheck_bip.artifact_vec_rejectvisual;
    data.label(end+1)={'artifact'};
    
    trlinfo=datainfo.trialinfo;
    % segment data in different trial parts
    % item window: -1 to 4
    pre_item=1;
    post_item=4;
    trl_item(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_item);
    trl_item(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl_item(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_item.*data.fsample);
    datainfo.trigger.trl_item=trl_item;
    
    pre_us=1;
    post_us=3;
    us_onset=datainfo.trigger.trigger_sp'+round(((trlinfo(:,13)-trlinfo(:,11))./10000).*data.fsample);
    trl_us(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_us);
    trl_us(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_us);
    trl_us(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_us.*data.fsample);
    datainfo.trigger.trl_us=trl_us;
    
    cfg=[];
    cfg.trl=trl_item;
    data_item=ft_redefinetrial(cfg,data);
    cfg=[];
    cfg.keeptrials='yes';
    erp_item=ft_timelockanalysis(cfg,data_item)
    artifreeitem_ind=squeeze(sum(erp_item.trial(:,end,:),3))==0;
    clear erp_item
    
    cfg=[];
    cfg.trl=trl_us;
    data_us=ft_redefinetrial(cfg,data);
    cfg=[];
    cfg.keeptrials='yes';
    erp_us=ft_timelockanalysis(cfg,data_us)
    artifreeus_ind=squeeze(sum(erp_us.trial(:,end,:),3))==0;
    clear erp_us
    
    item_trlinfo=datainfo.trialinfo(artifreeitem_ind,:);
    us_trlinfo=datainfo.trialinfo(artifreeus_ind,:);
    
    datainfo.artifact_info.clean_trials.item.trl=trl_item;
    datainfo.artifact_info.clean_trials.item.artifactfree=artifreeitem_ind;
    
    datainfo.artifact_info.clean_trials.us.trl=trl_us;
    datainfo.artifact_info.clean_trials.us.artifactfree=artifreeus_ind;
   save(info_file,'datainfo')
    
    for c=1:numel(conditions)
        % item trials
        trialcount_item(sub,c)=sum(item_trlinfo(:,2)==cond_def(c,1)&item_trlinfo(:,5)==cond_def(c,2));
        % us trials
        trialcount_us(sub,c)=sum(us_trlinfo(:,2)==cond_def(c,1)&us_trlinfo(:,5)==cond_def(c,2));
        
    end
    keep trialcount_item trialcount_us path_preproc path_info allsubs sub cond_def conditions
end

%% check location of remaining electrodes: bip processed data, checking notreref labels
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

all_atlas={'freesurferDK','aal','freesurferDestrieux'}%,}
for at=1:numel(all_atlas)
   sel_atlas=all_atlas{at};
   
    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        % electrodeinfo
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
        
   switch sel_atlas
       case 'freesurferDK'
           all_labels=datainfo.elec_info.ana_labels.freesurferDK; 
       case 'aal'
            all_labels=datainfo.elec_info.ana_labels.aal; 
       case 'freesurferDestrieux'
            all_labels=datainfo.elec_info.ana_labels.freesurferDestrieux; 
   end
        
        
        % channels surviving preprocessing
        clean_chan_bip=datainfo.artifact_info.rejectvisual_bip.elecsin;
        
        % clean chan bip, seperate channel names to get locations
        clean_chan=[];
        for i=1:numel(clean_chan_bip)
            sel_chan=clean_chan_bip{i};
            ind=strfind(sel_chan,'-');
            chan1={sel_chan(1:ind-1)};
            chan2={sel_chan(ind+1:end)};
            clean_chan=[clean_chan;chan1;chan2];
        end
        clean_chan=unique(clean_chan,'stable');
        
        % get indices of channel in analabel
        [chan,~,ind]=intersect(clean_chan,datainfo.elec_info.ana_labels.labels,'stable')
        % check whether each channel is matched
        if numel(chan)~=numel(clean_chan)
            error('channel labels do not match')
        end
        datainfo.artifact_info.rejectvisual_bip.clean_chan_nobip=clean_chan;
        save(info_file,'datainfo')
        
   
        % get analabels for the clean chan
        sel_labels=  all_labels(ind,:);       
        % goal: count of how many subjects have at least one contact in an area
        for r=1:size(sel_labels,2)
            all_labels_sub{sub,r}=unique([sel_labels{:,r}])';
        end
 end
 
% get all regions with an electrodes in at least one patient
all_areas=unique(vertcat(all_labels_sub{:,end}));

% count for every region how many patients have at least one electrode in
% the area
for r=1:size(sel_labels,2)
    for a=1:numel(all_areas)
        sel_area=all_areas{a};
        tmp_count=0;
        for sub=1:numel(allsubs)
            tmp_count=tmp_count+any(strcmp(all_labels_sub{sub,r},sel_area));
        end
        area_count{a,r}=tmp_count;
    end
end

all_areas=[all_areas,area_count];
table_areas=table(all_areas(:,1),all_areas(:,2),all_areas(:,3),all_areas(:,4),all_areas(:,5),all_areas(:,6),all_areas(:,7),....
    'VariableNames',{'label_freesurferDK','radius0','radius1','radius2','radius3','radius4','radius5'})

table_areas=sortrows(table_areas,2,'descend')
file_out=fullfile(path_info,strcat('patcountlabel_',sel_atlas));
save(file_out,'table_areas')
keep path_info path_artifactinfo path_preproc allsubs all_atlas at
end


%% add sampling rate to datainfo and check for lowest sampling rate
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
   % load(strcat(path_preproc,sel_sub,'_data.mat'))
   % datainfo.header=data.hdr;
  %save(info_file,'datainfo')

     sr(sub)=datainfo.header.Fs;
    keep path_info path_artifact path_preproc allsubs sub sr
end