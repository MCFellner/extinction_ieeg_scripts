path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\electrode_localization\';
path_rawdata='D:\Extinction\iEEG\rawdata\extinction_ieeg\';


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

for i=4:numel(allsubs)
sel_sub=allsubs{i};  
% electrode names from ieeg files
 eeg_file=strcat(path_data,'data\preproc\ieeg\readin\', sel_sub,'_data.mat');
load(eeg_file)

info_file=strcat(path_info,sel_sub,'\',sel_sub,'_datainfo_tmp');
load(info_file)

elec_info.labels_orderindata=data.label;
% resort labels for easier labeling
% define electrodes
ekg_count=0;
for i=1:numel(data.label)
    sel_lab=data.label{i};
    if strncmp(sel_lab, 'EKG',3)
    channel_type{i}='EKG';
    electrode{i,1}='EKG';
    ekg_count=ekg_count+1;
    electrode{i,2}=ekg_count;
    else
    channel_type{i}='iEEG';
    for j=1:numel(sel_lab)
    num_ind(j)=str2double(sel_lab(j));
    end
    electrode{i,1}=sel_lab(isnan(num_ind));
    electrode{i,2}=str2double(sel_lab(~isnan(num_ind)));
    clear num_ind
    end
end

cell_ieeg=sortrows(electrode(strcmp(channel_type,'iEEG'),:));

for i=1:size(cell_ieeg,1)
   sorted_ieeg{i,1}=strcat(cell_ieeg{i,1},num2str(cell_ieeg{i,2}));
end
datainfo.elec_info.labels_orderindata=data.label;
datainfo.elec_info.orderindata.label=data.label;
datainfo.elec_info.orderindata.channel_type=channel_type;
datainfo.elec_info.reordered.cell_ieeg=cell_ieeg;
datainfo.elec_info.reordered.sorted_ieeg=sorted_ieeg;


% if strcmp(sel_sub,'c_sub20')
%     clean_chan=[];
% else
% sel_folder=strcat(path_rawdata,sel_sub,'\ieeg\');
% cd(sel_folder)
% file_chan=dir('*_chaninfo.mat');
% load(file_chan.name)
% end
% data_label=data.label;
% 
% % check for name doubles
% 
% if numel(unique(data_label))~=numel(data_label)
% input('data_label ok?')
% end
% 
% if numel(unique(clean_chan))~=numel(clean_chan)
% input('clean_chan ok?')
% end
% 
% if numel(unique(datainfo.elec_info.labels_orderindata))~=numel(datainfo.elec_info.labels_orderindata)
% input('orderindata ok?')
% end
% 
% if numel(unique(datainfo.elec_info.orderindata.label))~=numel(datainfo.elec_info.orderindata.label)
% input('orderindata ok?')
% end
% 
% if numel(unique(datainfo.elec_info.reordered.sorted_ieeg))~=numel(datainfo.elec_info.reordered.sorted_ieeg)
% input('reordered ok?')
% [~,ind]=unique(datainfo.elec_info.reordered.sorted_ieeg);
% datainfo.elec_info.reordered.sorted_ieeg=datainfo.elec_info.reordered.sorted_ieeg(sort(ind));
% datainfo.elec_info.reordered.cell_ieeg=datainfo.elec_info.reordered.cell_ieeg(sort(ind),:);
% end
% 
% if numel(unique(datainfo.elec_info.elec_ct_mr.label))~=numel(datainfo.elec_info.elec_ct_mr.label)
% input('elec_ct_mr.label ok?')
% end
% 
% 
% % check clean_chan and data label
% indata_notinchaninfo=setdiff(data_label,clean_chan)
% inchaninfo_notindata=setdiff(clean_chan,data_label)
% 
% input('channels ok?')
% 
% % check' all labels, remove ghost electrodes
% indata_notinorderindata=setdiff(data_label,datainfo.elec_info.labels_orderindata)
% indata_notinorderindata=setdiff(data_label,datainfo.elec_info.orderindata.label)
% indata_notreordered=setdiff(data_label,datainfo.elec_info.reordered.sorted_ieeg)
% input('channels ok?')
% 
% inorderindata_notindata=setdiff(datainfo.elec_info.labels_orderindata,data_label)
% inorderindata_notindata=setdiff(datainfo.elec_info.orderindata.label,data_label)
% reordered_notindata=setdiff(datainfo.elec_info.reordered.sorted_ieeg,data_label)
% not_ok=input('channels ok?')
% if strcmp(not_ok,'ok')
% elseif strcmp(not_ok,'reordered')
% [~,ind]=setdiff(datainfo.elec_info.reordered.sorted_ieeg,data_label)
% datainfo.elec_info.reordered.sorted_ieeg(ind)=[];
% if length(datainfo.elec_info.reordered.cell_ieeg)==(numel(ind)+numel(datainfo.elec_info.reordered.sorted_ieeg))
% datainfo.elec_info.reordered.cell_ieeg(ind,:)=[];
% end
% for j=1:length(datainfo.elec_info.reordered.cell_ieeg)
% tmp{j}=strcat(datainfo.elec_info.reordered.cell_ieeg{j,1},num2str(datainfo.elec_info.reordered.cell_ieeg{j,2}))
% end
% 
% setdiff(tmp,datainfo.elec_info.reordered.sorted_ieeg)
% setdiff(datainfo.elec_info.reordered.sorted_ieeg,tmp)
% celllabel=input('reordered ok?')
% if strcmp(celllabel,'no')
% [~,ind]=setdiff(tmp,datainfo.elec_info.reordered.sorted_ieeg)
% datainfo.elec_info.reordered.cell_ieeg(ind,:)=[];
% end
% end


indata_notinelec_ct_mr=setdiff(data.label,datainfo.elec_info.elec_ct_mr.label)
inelec_ct_mr_notindata=setdiff(datainfo.elec_info.elec_ct_mr.label,data.label)
ct_mr=input('channels mr ct ok?')
if strcmp(ct_mr,'no')
 [~,ind] =setdiff(datainfo.elec_info.elec_ct_mr.label,data.label)
datainfo.elec_info.elec_ct_mr.label(ind)=[];
datainfo.elec_info.elec_ct_mr.elecpos(ind,:)=[];
datainfo.elec_info.elec_ct_mr.chanpos(ind,:)=[];
datainfo.elec_info.elec_ct_mr.tra(ind,:)=[];
datainfo.elec_info.elec_ct_mr.tra(:,ind)=[];
datainfo.elec_info.elec_ct_mr.cfg.channel=datainfo.elec_info.elec_ct_mr.label;
end


save(info_file,'datainfo')
%save(eeg_file,'data','-v7.3')
keep path_data path_info path_rawdata allsubs

end