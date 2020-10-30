addpath('E:\matlab_tools\fieldtrip-20190108')
ft_defaults

% edf readin

filename='E:\Extinction\iEEG\rawdata\China\QichenGroup\1-DBX\1-DBX\DBX-Learning test.edf';

% define trials
cfg            = [];
cfg.dataset    = filename;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);


% visually inspect the data
cfg            = [];
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);


% trigger channels: POL DC09-POL DC12