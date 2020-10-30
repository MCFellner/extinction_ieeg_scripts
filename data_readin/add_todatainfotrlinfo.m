%% add to trialinfo


path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08',...
            'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};


for i=1:numel(allsubs)
% save labels in file & as excel sheet
sel_sub=allsubs{i};
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

trlinfo=datainfo.trialinfo;
% info to add
% count presentation number of cond in each block
trlinfo(:,15)=zeros(length(trlinfo),1);
for b=1:3
for c=1:3
    sel_ind=find(trlinfo(:,5)==c &trlinfo(:,2)==b);
    trlinfo(sel_ind,15)=1:numel(sel_ind);
end
end
% count number of us for each cond
for c=1:3
    sel_ind=find(trlinfo(:,5)==c &trlinfo(:,9)==1);
    trlinfo(sel_ind,16)=1:numel(sel_ind);
end
datainfo.trialinfo=trlinfo;
save(info_file,'datainfo')
end