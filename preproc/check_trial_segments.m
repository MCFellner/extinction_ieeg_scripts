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
post_us=4;
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

for c=1:numel(conditions)
% item trials
trialcount_item(sub,c)=sum(item_trlinfo(:,2)==cond_def(c,1)&item_trlinfo(:,5)==cond_def(c,2));
% us trials
trialcount_us(sub,c)=sum(us_trlinfo(:,2)==cond_def(c,1)&us_trlinfo(:,5)==cond_def(c,2));

end
keep trialcount_item trialcount_us path_preproc path_info allsubs sub cond_def conditions
end


%% check elecs again
% elecs counts after preproc

% bipolar electrode contact positions