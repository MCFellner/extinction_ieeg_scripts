%% fix trialinfo in datainfo (get nous times for c phase)
path_trlinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

for sub=19%:length(allsubs)
    sel_sub=allsubs{sub};
    info_file=strcat(path_datainfo,sel_sub,'_datainfo');
    load(info_file) 
    load(strcat(path_trlinfo,sel_sub,'_trlinfo'))

% check datainfo trlinfo and trlinfo for mismatch,
tmp_trlinfo=trlinfo;
tmp_datainfo=datainfo.trialinfo;
tmp_trlinfo(isnan(tmp_trlinfo))=99;
tmp_datainfo(isnan(tmp_datainfo))=99;

[delete_row,delete_trlind]=setdiff(tmp_trlinfo(:,1:11),tmp_datainfo(:,1:11),'rows','stable')
if strcmp(sel_sub,'p_sub08')
trlinfo(133:134,:)=[];
end
datainfo.trialinfo=trlinfo;
save(info_file,'datainfo')
end
