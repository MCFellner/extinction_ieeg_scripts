%% copy files back from Q  (and fix electrode error)

% %datainfo, ct, mr
% path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
% path_q='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\electrode_localization\';
% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
% 
% for i=12%:numel(allsubs)
% sel_sub=allsubs{i};  
% 
% path_info=strcat(path_q,sel_sub,'\');
% info_file=strcat(path_info,sel_sub,'_datainfo_tmp.mat');
% load(info_file)
% 
% % check for wrong postions in datainfo.elec_info.reordered
% extra_labels=setdiff(datainfo.elec_info.reordered.sorted_ieeg,datainfo.elec_info.labels_orderindata);
% 
% if ~isempty(extra_labels)
%     display('fix electrodes')
%     data.label= datainfo.elec_info.labels_orderindata;
% % % resort labels for easier labeling
% % % define electrodes
%  ekg_count=0;
% for i=1:numel(data.label)
%      sel_lab=data.label{i};
%     if strncmp(sel_lab, 'EKG',3)
%      channel_type{i}='EKG';
%      electrode{i,1}='EKG';
%      ekg_count=ekg_count+1;
%     electrode{i,2}=ekg_count;
%      else
%      channel_type{i}='iEEG';
%      for j=1:numel(sel_lab)
%      num_ind(j)=str2double(sel_lab(j));
%      end
%      electrode{i,1}=sel_lab(isnan(num_ind));
%      electrode{i,2}=str2double(sel_lab(~isnan(num_ind)));
%      clear num_ind
%      end
%  end
%  
%  cell_ieeg=sortrows(electrode(strcmp(channel_type,'iEEG'),:));
%  
%  for i=1:size(cell_ieeg,1)
%     sorted_ieeg{i,1}=strcat(cell_ieeg{i,1},num2str(cell_ieeg{i,2}));
%  end
% datainfo.elec_info.orderindata.channel_type=channel_type;
% datainfo.elec_info.reordered.cell_ieeg=cell_ieeg;
% datainfo.elec_info.reordered.sorted_ieeg=sorted_ieeg;   
% clear sorted_ieeg cell_ieeg num_ind ekg_count channel_type electrode
% end
% 
% % check for not localizable electrodes:
% iEEG_ind=strcmp(datainfo.elec_info.orderindata.channel_type,'iEEG');
% check_reordered=sum(iEEG_ind)==numel(datainfo.elec_info.reordered.sorted_ieeg);
% if check_reordered==0
% error('elec reordered and iEEG do not match')
% end
% 
% non_found=setdiff(datainfo.elec_info.orderindata.label(iEEG_ind),datainfo.elec_info.elec_ct_mr.label)
% datainfo.elec_info.not_found_in_ct=non_found;
% save(strcat(path_datainfo,sel_sub,'_datainfo'), 'datainfo')
% end