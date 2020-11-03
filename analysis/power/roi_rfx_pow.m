addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')


%% roi specific rfx power cluster analysis

% get vp with elec in roi
% 
% roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
% roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};
% roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
% roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
% roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
% roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
% roi.amy_r={'Right-Amygdala'};
% roi.amy_l={'Left-Amygdala'};
% roi.hip_l={'Left-Hippocampus'};
% roi.hip_r={'Right-Hippocampus'};

% roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
% roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.ifg={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis','ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
%roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.dm_pfc ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal','ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
%roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};
roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};
roi.ventraltempocci={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole','ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};
%roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};

rois=fieldnames(roi);


% contrasts of interest
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\pow\rfx\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
stat_windows=[0 4.5;2 3;3 4;4 5]; % add rows for more windows
%stat_windows=[4 4.5];
% downsample to smallest sr
sr=1000;

% define all conditions to run

% conditions={'B_switch','B_plus'};
% cond1_def=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% cond2_def=[{'2'},{'==2'};{'6'},{'==1'}];
% conditions={'A_cs_2half','A_nocs_2half'};
% cond1_def=[{'2'},{'==1'};{'8'},{'==1'};{'15'},{'>=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% cond2_def=[{'2'},{'==1'};{'8'},{'==0'};{'15'},{'>=12'}];


% all_contrasts.Aplusplus_halfs.conditions={'Aplusplus_1half','Aplusplus_2half'};
% all_contrasts.Aplusplus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aplusplus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.Aplusplus_halfs.type='simple_con';
% 
% all_contrasts.Aplusminus_halfs.conditions={'Aplusminus_1half','Aplusminus_2half'};
% all_contrasts.Aplusminus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aplusminus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.Aplusminus_halfs.type='simple_con';
% 
% all_contrasts.Aminusminus_halfs.conditions={'Aminusminus_1half','Aminusminus_2half'};
% all_contrasts.Aminusminus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aminusminus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.Aminusminus_halfs.type='simple_con';
% 
% all_contrasts.Bplusplus_halfs.conditions={'Bplusplus_1half','Bplusplus_2half'};
% all_contrasts.Bplusplus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bplusplus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.Bplusplus_halfs.type='simple_con';
% 
% all_contrasts.Bplusminus_halfs.conditions={'Bplusminus_1half','Bplusminus_2half'};
% all_contrasts.Bplusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bplusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.Bplusminus_halfs.type='simple_con';
% 
% all_contrasts.Bminusminus_halfs.conditions={'Bminusminus_1half','Bminusminus_2half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% 
% all_contrasts.int_Aplusminus_Aminusminus_halfs.conditions={'Aplusminus_1half','Aplussminus_2half';'Aminusminus_1half','Aminusminus_2half'};
% all_contrasts.int_Aplusminus_Aminusminus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aplusminus_Aminusminus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.int_Aplusminus_Aminusminus_halfs.cond_def{3}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aplusminus_Aminusminus_halfs.cond_def{4}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.int_Aplusminus_Aminusminus_halfs.type='interaction';
% 
% 
% all_contrasts.int_Bplusminus_Bminusminus_halfs.conditions={'Bplusminus_1half','Bplussminus_2half';'Bminusminus_1half','Bminusminus_2half'};
% all_contrasts.int_Bplusminus_Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bplusminus_Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.int_Bplusminus_Bminusminus_halfs.cond_def{3}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bplusminus_Bminusminus_halfs.cond_def{4}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.int_Bplusminus_Bminusminus_halfs.type='interaction';
% 
% 
% all_contrasts.int_Aplusminus_Aplusplus_halfs.conditions={'Aplusminus_1half','Aplusplus_2half';'Aplusplus_1half','Aminusminus_2half'};
% all_contrasts.int_Aplusminus_Aplusplus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aplusminus_Aplusplus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.int_Aplusminus_Aplusplus_halfs.cond_def{3}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aplusminus_Aplusplus_halfs.cond_def{4}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.int_Aplusminus_Aplusplus_halfs.type='interaction';
% 
% 
% all_contrasts.int_Bplusminus_Bplusplus_halfs.conditions={'Bplusminus_1half','Bplusplus_2half';'Bplusplus_1half','Bminusminus_2half'};
% all_contrasts.int_Bplusminus_Bplusplus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bplusminus_Bplusplus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];
% all_contrasts.int_Bplusminus_Bplusplus_halfs.cond_def{3}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bplusminus_Bplusplus_halfs.cond_def{4}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.int_Bplusminus_Bplusplus_halfs.type='interaction';
% 
% 
% all_contrasts.int_Aminusminus_Aplusplus_halfs.conditions={'Aminusminus_1half','Aminusminus_2half';'Aplusplus_1half','Aminusminus_2half'};
% all_contrasts.int_Aminusminus_Aplusplus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aminusminus_Aplusplus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.int_Aminusminus_Aplusplus_halfs.cond_def{3}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Aminusminus_Aplusplus_halfs.cond_def{4}=[{'2'},{'==1'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.int_Aminusminus_Aplusplus_halfs.type='interaction';
% 
% 
% all_contrasts.int_Bminusminus_Bplusplus_halfs.conditions={'Bminusminus_1half','Bminusminus_2half';'Bplusplus_1half','Bminusminus_2half'};
% all_contrasts.int_Bminusminus_Bplusplus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bminusminus_Bplusplus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.int_Bminusminus_Bplusplus_halfs.cond_def{3}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_Bminusminus_Bplusplus_halfs.cond_def{4}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.int_Bminusminus_Bplusplus_halfs.type='interaction';

% 
% all_contrasts.Bminusminus_halfs.conditions={'Aplus_2half','Aminus_2half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% 
% all_contrasts.Bminusminus_halfs.conditions={'Aplus_1half','Aminus_1half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% all_contrasts.Bminusminus_halfs.conditions={'Bswitch_2half','Bminusminus_2half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% 
% all_contrasts.Bminusminus_halfs.conditions={'Bswitch_1half','Bminusminus_1half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% all_contrasts.Bminusminus_halfs.conditions={'Bswitch_2half','Bplusplus_2half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% 
% all_contrasts.Bminusminus_halfs.conditions={'Bswitch_1half','Bplusplus_1half'};
% all_contrasts.Bminusminus_halfs.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Bminusminus_halfs.cond_def{2}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];
% all_contrasts.Bminusminus_halfs.type='simple_con';
% 
% all_contrasts.Aus_nous.conditions={'Aus','Ano_us'};
% all_contrasts.Aus_nous.cond_def{1}=[{'2'},{'==1'};{'9'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aus_nous.cond_def{2}=[{'2'},{'==1'};{'9'},{'==0'}];
% all_contrasts.Aus_nous.type='simple_con';

% all_contrasts.Aplus_halfs.conditions={'Aplus_1half','Aplus_2half'};
% all_contrasts.Aplus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.Aplus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'>12'}];
% all_contrasts.Aplus_halfs.type='simple_con';
% 
% all_contrasts.int_A_plus_minus_halfs.conditions={'Aplus_1half','Aplus_2half','Bminusminus_1half','Bminusminus_2half'};
% all_contrasts.int_A_plus_minus_halfs.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_A_plus_minus_halfs.cond_def{2}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'>12'}];
% all_contrasts.int_A_plus_minus_halfs.cond_def{3}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
% all_contrasts.int_A_plus_minus_halfs.cond_def{4}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
% all_contrasts.int_A_plus_minus_halfs.type='interaction';


all_contrasts.Aplus_minus.conditions={'Aplus','Aminus'};
all_contrasts.Aplus_minus.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Aplus_minus.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'}];
all_contrasts.Aplus_minus.type='simple_con';

all_contrasts.Bplusplus_Bplusminus.conditions={'Bplusplus','Bplusminus'};
all_contrasts.Bplusplus_Bplusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bplusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'}];
all_contrasts.Bplusplus_Bplusminus.type='simple_con';

all_contrasts.Bplusminus_Bminusminus.conditions={'Bplusminus','Bminusminus'};
all_contrasts.Bplusminus_Bminusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusminus_Bminusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'}];
all_contrasts.Bplusminus_Bminusminus.type='simple_con';


all_contrasts.Bplusplus_Bminusminus.conditions={'Bplusplus','Bminusminus'};
all_contrasts.Bplusplus_Bminusminus.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bminusminus.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'}];
all_contrasts.Bplusplus_Bminusminus.type='simple_con';


all_contrasts.Aplus_minus_2half.conditions={'Aplus_2half','Aminus_2half'};
all_contrasts.Aplus_minus_2half.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Aplus_minus_2half.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'>12'}];
all_contrasts.Aplus_minus_2half.type='simple_con';

all_contrasts.Bplusplus_Bplusminus_2half.conditions={'Bplusplus_2half','Bplusminus_2half'};
all_contrasts.Bplusplus_Bplusminus_2half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bplusminus_2half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];
all_contrasts.Bplusplus_Bplusminus_2half.type='simple_con';

all_contrasts.Bplusminus_Bminusminus_2half.conditions={'Bplusminus_2half','Bminusminus_2half'};
all_contrasts.Bplusminus_Bminusminus_2half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusminus_Bminusminus_2half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
all_contrasts.Bplusminus_Bminusminus_2half.type='simple_con';


all_contrasts.Bplusplus_Bminusminus_2half.conditions={'Bplusplus_2half','Bminusminus_2half'};
all_contrasts.Bplusplus_Bminusminus_2half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'>12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bminusminus_2half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'>12'}];
all_contrasts.Bplusplus_Bminusminus_2half.type='simple_con';

all_contrasts.Aplus_minus_1half.conditions={'Aplus_1half','Aminus_1half'};
all_contrasts.Aplus_minus_1half.cond_def{1}=[{'2'},{'==1'};{'6'},{'<=2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Aplus_minus_1half.cond_def{2}=[{'2'},{'==1'};{'6'},{'==3'};{'15'},{'<=12'}];
all_contrasts.Aplus_minus_1half.type='simple_con';

all_contrasts.Bplusplus_Bplusminus_1half.conditions={'Bplusplus_1half','Bplusminus_1half'};
all_contrasts.Bplusplus_Bplusminus_1half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bplusminus_1half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];
all_contrasts.Bplusplus_Bplusminus_1half.type='simple_con';

all_contrasts.Bplusminus_Bminusminus_1half.conditions={'Bplusminus_1half','Bminusminus_1half'};
all_contrasts.Bplusminus_Bminusminus_1half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==2'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusminus_Bminusminus_1half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];
all_contrasts.Bplusminus_Bminusminus_1half.type='simple_con';


all_contrasts.Bplusplus_Bminusminus_1half.conditions={'Bplusplus_1half','Bminusminus_1half'};
all_contrasts.Bplusplus_Bminusminus_1half.cond_def{1}=[{'2'},{'==2'};{'6'},{'==1'};{'15'},{'<=12'}];% column, value (through eval also <= or ~=), definition across columns combined with &
all_contrasts.Bplusplus_Bminusminus_1half.cond_def{2}=[{'2'},{'==2'};{'6'},{'==3'};{'15'},{'<=12'}];
all_contrasts.Bplusplus_Bminusminus_1half.type='simple_con';




all_cons=fieldnames(all_contrasts);



freqs={'lf','hf'}
for f=1:numel(freqs)
    sel_freq=freqs{f};
    
    % get electrodes for each sub/roi
    for r=1:numel(rois)
        sel_roi=rois{r};
        roi_def=getfield(roi,sel_roi);
        sub_num=0;
        for sub=1:length(allsubs)
            sel_sub=allsubs{sub};
            % electrodeinfo
            info_file=strcat(path_info,sel_sub,'_datainfo');
            load(info_file)
            sel_elec_tmp=datainfo.elec_info.bipolar.elec_ct_mr.label(ismember([datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK{:}],[roi_def]));
            % only select after preproc elecs
            sel_elec=intersect(sel_elec_tmp,datainfo.artifact_info.rejectvisual_bip.elecsin);
            
            count_elec(sub,r)=numel(sel_elec);
            
            if numel(sel_elec)>0
                sub_num=sub_num+1;
                
                load(strcat(path_preproc,sel_sub,'_data.mat'))
                
                
                % apply bipolar montage
                montage=datainfo.elec_info.bipolar.montage_withoutartichan;
                [~,sel_ind]=intersect(montage.labelnew,sel_elec,'stable');
                montage.labelnew=montage.labelnew(sel_ind);
                montage.tra=montage.tra(sel_ind,:);
                data = ft_apply_montage(data,montage);
                
                
                % cut trials
                trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
                trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
                trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
                cfg=[];
                cfg.trl=trl;
                data=ft_redefinetrial(cfg,data);
                clear trl
                % downsample to common sampling rate
                cfg=[];
                cfg.resamplefs      = sr;
                cfg.detrend='yes';
                data=ft_resampledata(cfg,data);
                
                % only select artifree trials in trialinfo
                trlinfo=datainfo.trialinfo;
                trlinfo=trlinfo(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree,:);
                
                cfg=[]
                cfg.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
                data=ft_selectdata(cfg,data);
                
                
                cfg=[];
                cfg.output='pow';
                cfg.toi =pre_item:0.05:post_item;
                cfg.keeptrials  = 'yes';
                cfg.pad='nextpow2';
                switch sel_freq
                    case 'lf'
                        cfg.foi     = 1:1:30;
                        cfg.method='wavelet';
                        cfg.width = 7;
                    case 'hf'
                        cfg.method      = 'mtmconvol';
                        cfg.taper       = 'dpss';
                        cfg.foi         = 30:5:200;
                        cfg.t_ftimwin   = ones(1,length(cfg.foi))*0.3;
                        cfg.tapsmofrq   = ones(1,length(cfg.foi))*10;
                end
                freq=ft_freqanalysis(cfg,data);
                % z trans
                cfg=[];
                cfg.time=[pre_item+0.5,post_item-0.5];
                freq=z_trans_TF_seltime(cfg, freq);
                
                
                
               for con=1:numel(all_cons)
                    % seperate in condition 1 & 2 & average over freqs
                    sel_con=all_cons{con};
                    sel_conditions=getfield(all_contrasts,sel_con);
                    conditions=getfield(sel_conditions,'conditions');
                    cond_def=getfield(sel_conditions,'cond_def');
                    
                    % define two conditions
                    for c=1:numel(cond_def)
                        sel_cond_def=cond_def{c};
                        for d=1:size(sel_cond_def,1)
                            eval(strcat('tmp(:,d)=trlinfo(:,',sel_cond_def{d,1},')',sel_cond_def{d,2}));
                        end
                        sel_trials=find(sum(tmp,2)==size(sel_cond_def,1));
                        clear tmp
                        
                        
                        freq_all{con,c,sub_num}=freq;
                        freq_all{con,c,sub_num}.powspctrm=squeeze(nanmean(freq.powspctrm(sel_trials,:,:,:),1));
                        freq_all{con,c,sub_num}.powspctrm=nanmean(freq_all{con,c,sub_num}.powspctrm(:,:,:),1);
                        freq_all{con,c,sub_num}.cumtapcnt=nanmean(freq.cumtapcnt(sel_trials,:),1);
                        freq_all{con,c,sub_num}.dimord='chan_freq_time';
                        freq_all{con,c,sub_num}.label={sel_roi};
                        
                        
                    end
                end
            else
            end
            clear freq data datainfo montage sel_elec sel_elec_tmp sel_ind sel_trials  trl trlinfo
            
        end
 
        for con=1:numel(all_cons)
            % seperate in condition 1 & 2 & average over freqs
            sel_con=all_cons{con};
            sel_conditions=getfield(all_contrasts,sel_con);
            conditions=getfield(sel_conditions,'conditions');
            type=getfield(sel_conditions,'type');
            
            for c=1:numel(conditions)
            cfg=[];
            cfg.keepindividual='yes';
            ga{c}=ft_freqgrandaverage(cfg,freq_all{con,c,:});
            end
            
            switch type
                case 'interaction'
                    ga1=ga{1};
                    ga2=ga{1};
                    ga1.powspctrm=ga{1}.powspctrm-ga{2}.powspctrm;
                    ga2.powspctrm=ga{3}.powspctrm-ga{4}.powspctrm;
                    
                case 'simple_con'
                    ga1=ga{1};
                    ga2=ga{2};
            end
            clear ga
            
            for w=1:size(stat_windows)
                sel_window=stat_windows(w,:);
                
                cfg=[];
                cfg.latency     = sel_window;
                cfg.avgoverchan =  'yes';
                cfg.avgovertime =  'no';
                cfg.avgoverfreq =  'no';
                
                % first level
                cfg.method           = 'montecarlo';
                cfg.numrandomization = 1000;
                cfg.correctm         =  'cluster';
                cfg.correcttail      = 'prob';
                cfg.statistic ='depsamplesT'
                % for within-subjects (depsamplesT)
                Nsub = size(ga1.powspctrm,1);                                       %# of subjects?
                design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
                design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
                
                cfg.uvar     = 2;
                cfg.ivar     = 1;
                cfg.design = design;
                stat=ft_freqstatistics (cfg,ga1,ga2)
                stat.sel_conditions=sel_conditions;
                mkdir(fullfile(path_out,strcat(sel_con),sel_roi,sel_freq))
                save(fullfile(path_out,strcat(sel_con),sel_roi,sel_freq,strcat('clusterstat_',num2str(sel_window(1).*1000),'to', num2str(sel_window(2).*1000))),'stat')
                clear design
                
            end
        end
        clear freq_all
    end
end




