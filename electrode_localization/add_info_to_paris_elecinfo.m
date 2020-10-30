%% add info to paris electrode info

% info about electrodes

% info about the channeltype
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};

for n=1:numel(allsubs)
sel_sub=allsubs{n};  

% electrodeinfo
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

% there are some electrode inconsistencies (electrodes in epiloc not in
% data and vice versa), these issues are fixed here

iEEG_labels=datainfo.elec_info.elec_ct_mr.label;
inctmr_not_indata=setdiff(datainfo.elec_info.elec_ct_mr.label, datainfo.elec_info.orderindata.label)
% delete channels  inctmr_not_indata (no need for the elecpos then)
if ~isempty(inctmr_not_indata)
     %delete in elec_ct_mr
     [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.elec_ct_mr.label)
     datainfo.elec_info.elec_ct_mr.label(ind)=[];
     datainfo.elec_info.elec_ct_mr.chanpos(ind,:)=[];
     datainfo.elec_info.elec_ct_mr.elecpos(ind,:)=[];
     datainfo.elec_info.elec_ct_mr.tra(ind,:)=[];
     datainfo.elec_info.elec_ct_mr.tra(:,ind)=[];
     
    %delete in elec_mni
     [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.elec_mni.label)
     datainfo.elec_info.elec_mni.label(ind)=[];
     datainfo.elec_info.elec_mni.chanpos(ind,:)=[];
     datainfo.elec_info.elec_mni.elecpos(ind,:)=[];
     datainfo.elec_info.elec_mni.tra(ind,:)=[];
     datainfo.elec_info.elec_mni.tra(:,ind)=[];
    %delete in ana_labels
     [~,~,ind]= intersect(inctmr_not_indata,datainfo.elec_info.ana_labels.labels)
     datainfo.elec_info.ana_labels.freesurferDK(ind,:)=[];
     datainfo.elec_info.ana_labels.freesurferDestrieux(ind,:)=[];
     datainfo.elec_info.ana_labels.aal(ind,:)=[];
     datainfo.elec_info.ana_labels.afni(ind,:)=[];
     datainfo.elec_info.ana_labels.brainweb(ind,:)=[];
     
      datainfo.elec_info.ana_labels.freesurferDK_def.elec_def=datainfo.elec_info.elec_ct_mr;
      datainfo.elec_info.ana_labels.freesurferDestrieux_def.elec_def=datainfo.elec_info.elec_ct_mr;
      datainfo.elec_info.ana_labels.aal_def.elec_def=datainfo.elec_info.elec_mni;
      datainfo.elec_info.ana_labels.afni_def.elec_def=datainfo.elec_info.elec_mni;
      datainfo.elec_info.ana_labels.brainweb_def.elec_def=datainfo.elec_info.elec_mni;
    end
  datainfo.elec_info.ana_labels.labels=datainfo.elec_info.elec_ct_mr.label;

indata_not_inctmr=setdiff(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label)
% define channels indata_not_inctmr specifically in channeltype
% 'nopos_iEEG'
    % define datainfo.elec_info.orderindata.channel_type
% get channel_type for indata_not_inctmr
    [~,~,ind]=intersect(indata_not_inctmr,datainfo.elec_info.orderindata.label)
    chantype_indata_not_inctmr=datainfo.elec_info.orderindata.channel_type(ind);
    [num2cell(ind),indata_not_inctmr,chantype_indata_not_inctmr]
    nopos_iEEG_ind=input('Which channels are iEEG?, type indices[]');
    if ~isempty(nopos_iEEG_ind)
       datainfo.elec_info.orderindata.channel_type(nopos_iEEG_ind)={'nopos_iEEG'};
    end

% now channels are consistent, fix datainfo.elec_info.orderindata.channel_type
    % set for all labels with elecpos (aka iEEG) the channel type label 'iEEG'
  [~,ind,~]=intersect(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label);
  datainfo.elec_info.orderindata.channel_type(ind)={'iEEG'};
     [datainfo.elec_info.orderindata.label,datainfo.elec_info.orderindata.channel_type]
    check=input('channel_types generally ok?')
    
% add datainfo.elec_info.reordered (only for iEEG channels)
    ieeg_label=datainfo.elec_info.orderindata.label(strcmp('iEEG',datainfo.elec_info.orderindata.channel_type));
    %add datainfo.elec_info.reordered.cell: first column name of electrode, 2 name of contact  
    for e=1:numel(ieeg_label)
    sel_label=ieeg_label{e};   
    ind_=strfind(sel_label,'_');
    datainfo.elec_info.reordered.cell_ieeg{e,1}=sel_label(1:ind_-1);
    datainfo.elec_info.reordered.cell_ieeg{e,2}=str2num(sel_label(ind_+1:end));
    end
  [datainfo.elec_info.reordered.cell_ieeg,ind]=  sortrows( datainfo.elec_info.reordered.cell_ieeg);
    
    datainfo.elec_info.reordered.sorted_ieeg=ieeg_label(ind);
    [datainfo.elec_info.reordered.sorted_ieeg,datainfo.elec_info.reordered.cell_ieeg]
    check=input('are electrodes and numbers well seperated?')
    
save(info_file,'datainfo')
clear datainfo ind sel_label ind_ chantype_indata_not inctmr check e ieeg_label iEEG_labels inctmr_not_indata indata_not_inctmr nopo_iEEG_ind
end