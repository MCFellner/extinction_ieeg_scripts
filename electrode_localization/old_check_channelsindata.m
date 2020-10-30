addpath('T:\hanna\fieldtrip-20190611')
ft_defaults


%% check channels in raw data
path_data='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\rawdata\extinction_ieeg\';
path_preproc='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\preproc\ieeg\readin\';
path_info='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\electrode_localization\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

for sub=1%:length(allsubs)
sel_sub=allsubs{sub};
cd(strcat(path_info,sel_sub))
% electrodeinfo
info_file=strcat(sel_sub,'_datainfo');
load(info_file)

% first load readin data
display(strcat('selected data of: ',sel_sub))
load(strcat(path_preproc,sel_sub,'_data.mat'))

all_elec=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));

for e=1:numel(all_elec)
sel_channel=datainfo.elec_info.reordered.sorted_ieeg(strcmp(datainfo.elec_info.reordered.cell_ieeg(:,1),all_elec{e}));
cfg            = [];
cfg.viewmode   = 'vertical';
cfg.ylim   =  'maxmin';
cfg.channel    = sel_channel;
cfg.blocksize  =30;
cfg.continuous              = 'yes';
ft_databrowser(cfg, data);

electrode_check.label=all_elec;
electrode_check.all_ok{e}=input(strcat('Is electrode ',all_elec{e},' generally ok? write [1] if ok, [0] if not'))
electrode_check.flat_electrode{e}=input(strcat('Is there a flat channel on electrode ',all_elec{e},'? write electrode like that: [10] or [] if none'))
electrode_check.dampened_electrode{e}=input(strcat('Is there a dampened channel on electrode ',all_elec{e},'? write electrode like that: [10] or [] if none'))
electrode_check.artifact_electrode{e}=input(strcat('Is there a artifact channel on electrode ',all_elec{e},'? write electrode like that: [10] or [] if none'))

end
save(strcat(path_info,sel_sub,'_firstelectrodecheck'));

end

%% for datasets with missing electrodes check original data, for ok electrodes

path_data='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\rawdata\extinction_ieeg\';
path_preproc='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\preproc\ieeg\readin\';
path_info='Q:\Ongoing_projects_Fellner_MC\Extinction\iEEG\data\electrode_localization\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};

for sub=1%:length(allsubs)
sel_sub=allsubs{sub};
cd(strcat(path_info,sel_sub))
sel_folder=strcat(path_data,sel_sub,'\ieeg\');
cd (sel_folder)
filename=dir('*.edf');

% define trials
cfg            = [];
cfg.dataset    = filename.name;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data_all           = ft_preprocessing(cfg);

cfg            = [];
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data_all);
end