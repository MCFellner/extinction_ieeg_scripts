% check for correct channeltype definition

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
        'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

 for i=1:numel(allsubs)
% save labels in file & as excel sheet
sel_sub=allsubs{i};
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)
% check for correct datainfo
check=numel(datainfo.elec_info.orderindata.label)==numel(datainfo.elec_info.orderindata.channel_type)
if check==0
[~,iEEGind,~]=intersect(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label,'stable')
channel_type(iEEGind)={'iEEG'};

[~,noiEEGind]=setdiff(datainfo.elec_info.orderindata.label,datainfo.elec_info.elec_ct_mr.label,'stable')
noiEEG_label=[datainfo.elec_info.orderindata.label(noiEEGind)]
ekg=input('are the channels EKG?');
if ekg
channel_type(noiEEGind)={'ecg'};
end
datainfo.elec_info.orderindata.channel_type=channel_type;
save(info_file,'datainfo')
end
      
end