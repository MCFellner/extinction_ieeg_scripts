path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
path_trialinfo='D:\Extinction\iEEG\data\preproc\trialinfo\';
path_out='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_datainfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
mkdir(path_datainfo)
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18'};

for sub=2:length(allsubs)
sel_sub=allsubs{sub};

sel_folder=strcat(path_data,sel_sub,'\ieeg\');
cd (sel_folder)
filename=dir('*.edf');

% define trials
cfg            = [];
cfg.dataset    = filename.name;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data_all           = ft_preprocessing(cfg);



% get rid of 'POL' in label names
for i=1:numel(data_all.label)
   data_all.label{i}=data_all.label{i}(5:end); 
   if isempty(data_all.label{i})
       data_all.label{i}='POL';
   end
end

if strcmp(sel_sub, 'c_sub03')|| strcmp(sel_sub, 'c_sub04') || strcmp(sel_sub, 'c_sub05') || strcmp(sel_sub, 'c_sub06')  || strcmp(sel_sub, 'c_sub08')...
        || strcmp(sel_sub, 'c_sub11') || strcmp(sel_sub, 'c_sub12')|| strcmp(sel_sub, 'c_sub14') || strcmp(sel_sub, 'c_sub15') || strcmp(sel_sub, 'c_sub18')
 for i=1:numel(data_all.label)
    refmatch=strfind(data_all.label{i},'-Ref');
   if    ~isempty(refmatch)
       data_all.label{i}=data_all.label{i}(1:refmatch-1);
   end
end
end



file_chan=dir('*_chaninfo.mat');
load(file_chan.name)
save(strcat(file_chan.name,'_old'),'clean_chan')
clear clean_chan

clean_chan=data_all.label;
save(file_chan.name,'clean_chan')
end