addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')


%% learn curve analysis
%
% how to
% get power (separete lf/hf)
% do regression for every electrode/condition: yield beta value map

% save beta value maps for each subject/conditions

% conditions: linear regressor type 1-3/block A-C
% response: response on trial? change on trial?

% extract power for learncurvees

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\pow\regression\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
stat_windows=[2 3; 3  4]; % add rows for more windows

% downsample to smallest sr
sr=1000;

% condition definitions: conditions struct with fields defining reg/trl for
% each condition

conditions.type1_A.trl_def=[{'2'},{'==1'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_A.reg_def=['15'];

conditions.type2_A.trl_def=[{'2'},{'==1'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_A.reg_def=['15'];

conditions.type3_A.trl_def=[{'2'},{'==1'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_A.reg_def=['15'];

conditions.type1_B.trl_def=[{'2'},{'==2'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_B.reg_def=['15'];

conditions.type2_B.trl_def=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_B.reg_def=['15'];

conditions.type3_B.trl_def=[{'2'},{'==2'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_B.reg_def=['15'];

conditions.type1_C.trl_def=[{'2'},{'==3'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_C.reg_def=['15'];

conditions.type2_C.trl_def=[{'2'},{'==3'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_C.reg_def=['15'];

conditions.type3_C.trl_def=[{'2'},{'==3'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_C.reg_def=['15'];


freqs={'hf'};%'lf',
for f=1:numel(freqs)
    freq_def=freqs{f};

    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        % electrodeinfo
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
        load(strcat(path_preproc,sel_sub,'_data.mat'))

        % apply bipolar montage
        montage=datainfo.elec_info.bipolar.montage_withoutartichan;
        data = ft_apply_montage(data,montage);


        % cut trials
        trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
        trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
        trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
        cfg=[];
        cfg.trl=trl;
        data=ft_redefinetrial(cfg,data);

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
        switch freq_def
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

        conds=fieldnames(conditions);
        % select trials for cond
        for c=1:numel(conds)
            sel_cond=conds{c};
            cond_def=getfield(conditions,sel_cond);
            %sel_trials
            % sel_reg
            for d=1:size(cond_def.trl_def,1)
                eval(strcat('tmp(:,d)=trlinfo(:,',cond_def.trl_def{d,1},')',cond_def.trl_def{d,2}));
            end
            sel_trials=find(sum(tmp,2)==size(cond_def.trl_def,1));
            clear tmp
           eval(strcat('sel_reg=trlinfo(sel_trials,',cond_def.reg_def,');'))

            t_map=zeros(numel(freq.label),numel(freq.freq),numel(freq.time));
            p_pos=zeros(numel(freq.label),1);
            p_neg=zeros(numel(freq.label),1);

            for e=1:numel(freq.label)

                sel_freq=freq;
                sel_freq.powspctrm=zscore(freq.powspctrm(sel_trials,e,:,:));
                sel_freq.label=freq.label(e);

                cfg=[];
                cfg.avgoverchan =  'yes';
                cfg.avgovertime =  'no';
                cfg.avgoverfreq =  'no';
                cfg.method           = 'montecarlo';
                cfg.statistic ='indepsamplesregrT';
                cfg.numrandomization = 1000;
                cfg.correctm         =  'cluster';
                cfg.correcttail      = 'prob';
                cfg.design = sel_reg;
                cfg.ivar=1;
                cfg.cvar=[];
                stat=ft_freqstatistics (cfg,sel_freq)

                % get p values for each electrode
                if isfield(stat,'posclusters')
                    if ~isempty(stat.posclusters)
                        p_pos(e)=stat.posclusters(1).prob;
                    else
                        p_pos(e)=1;
                    end
                else
                    p_pos(e)=1
                end

                if isfield(stat,'negclusters')
                    if ~isempty(stat.negclusters)
                        p_neg(e)=stat.negclusters(1).prob;
                    else
                        p_neg(e)=1;
                    end
                else
                    p_neg(e)=1;
                end

                t_map(e,:,:)=stat.stat;


            end
           clear sel_reg sel_trials
           freq_reg=freq;
           freq_reg.powspctrm=t_map;
           freq_reg.dimord='chan_freq_time';
           freq_reg.cfg_reg=stat.cfg;
           freg_reg.cond_def=cond_def;
           sel_path=fullfile(path_out,sel_cond);
           mkdir(sel_path)
           save(fullfile(sel_path,[sel_sub,'_',freq_def,'_reg_tmap']),'freq_reg')
        end
    end
end
%% elec wise test

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\pow\regression\elecwise\';
mkdir(path_out)


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=-1;
post_item=5.5;
stat_windows=[2 2.5; 2.5 3; 3 3.5; 3.5 4]; % add rows for more windows
fois=[2 4;5 7;8 12;13 20;20 30; 35  60; 70 100];
foi_names={'theta1','theta2','alpha','beta1','beta2','gamma1','gamma2'};
sr=1000;
% condition definitions: conditions struct with fields defining reg/trl for
% each condition

conditions.type1_A.trl_def=[{'2'},{'==1'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_A.reg_def=['15'];

conditions.type2_A.trl_def=[{'2'},{'==1'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_A.reg_def=['15'];

conditions.type3_A.trl_def=[{'2'},{'==1'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_A.reg_def=['15'];

conditions.type1_B.trl_def=[{'2'},{'==2'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_B.reg_def=['15'];

conditions.type2_B.trl_def=[{'2'},{'==2'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_B.reg_def=['15'];

conditions.type3_B.trl_def=[{'2'},{'==2'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_B.reg_def=['15'];

conditions.type1_C.trl_def=[{'2'},{'==3'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type1_C.reg_def=['15'];

conditions.type2_C.trl_def=[{'2'},{'==3'};{'6'},{'==2'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type2_C.reg_def=['15'];

conditions.type3_C.trl_def=[{'2'},{'==3'};{'6'},{'==3'}];% column, value (through eval also <= or ~=), definition across columns combined with &
conditions.type3_C.reg_def=['15'];




for sub=27%1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);
    
    
    % cut trials
    trl(:,1)=datainfo.trigger.trigger_sp+(data.fsample.*pre_item);
    trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(1.*pre_item.*data.fsample);
    cfg=[];
    cfg.trl=trl;
    data=ft_redefinetrial(cfg,data);
    
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
    
    
    for f=1:length(fois)
        foi=fois(f,:);
        f_name=foi_names{f};
        
        cfg=[];
        cfg.output='pow';
        cfg.toi =pre_item:0.05:post_item;
        cfg.keeptrials  = 'yes';
        cfg.pad='nextpow2';
        switch f_name
            case {'theta1','theta2','alpha','beta1','beta2'}
                cfg.foilim     = foi;
                cfg.method='wavelet';
                cfg.width = 5;
            case {'gamma1','gamma2'}
                cfg.method      = 'mtmconvol';
                cfg.taper       = 'dpss';
                cfg.foi         = foi(1):5:foi(2);
                cfg.t_ftimwin   = ones(1,length(cfg.foi))*0.3;
                cfg.tapsmofrq   = ones(1,length(cfg.foi))*10;
        end
        freq=ft_freqanalysis(cfg,data);
        
        % z trans
        cfg=[];
        cfg.time=[pre_item+0.5,post_item-0.5];
        freq=z_trans_TF_seltime(cfg, freq);
        
        conds=fieldnames(conditions);
        % select trials for cond
        for c=1:numel(conds)
            sel_cond=conds{c};
            cond_def=getfield(conditions,sel_cond);
            %sel_trials
            % sel_reg
            for d=1:size(cond_def.trl_def,1)
                eval(strcat('tmp(:,d)=trlinfo(:,',cond_def.trl_def{d,1},')',cond_def.trl_def{d,2}));
            end
            sel_trials=find(sum(tmp,2)==size(cond_def.trl_def,1));
            clear tmp
            eval(strcat('sel_reg=trlinfo(sel_trials,',cond_def.reg_def,');'))
            
            
            for t=1:size(stat_windows,1)
                sel_window=   stat_windows(t,:);
                p_pos=zeros(numel(freq.label),1);
                p_neg=zeros(numel(freq.label),1);
                
                parfor e=1:numel(freq.label)
                    
                    sel_freq=freq;
                    sel_freq.powspctrm=zscore(freq.powspctrm(sel_trials,e,:,:));
                    sel_freq.label=freq.label(e);
                    
                    cfg=[];
                    cfg.avgoverchan =  'yes';
                    cfg.avgovertime =  'no';
                    cfg.avgoverfreq =  'yes';
                    cfg.latency     = sel_window;
                    cfg.method           = 'montecarlo';
                    cfg.statistic ='indepsamplesregrT';
                    cfg.numrandomization = 1000;
                    cfg.correctm         =  'cluster';
                    cfg.correcttail      = 'prob';
                    cfg.design = sel_reg;
                    cfg.ivar=1;
                    cfg.cvar=[];
                    stat=ft_freqstatistics (cfg,sel_freq)
                    
                    if isfield(stat,'posclusters')
                        if ~isempty(stat.posclusters)
                            pos_prob(e)=stat.posclusters(1).prob;
                        else
                            pos_prob(e)=1;
                        end
                    else
                        pos_prob(e)=1;
                    end
                    
                    if isfield(stat, 'negclusters')
                        if ~isempty(stat.negclusters)
                            neg_prob(e)=stat.negclusters(1).prob;
                        else
                            neg_prob(e)=1;
                        end
                    else
                        neg_prob(e)=1;
                    end
                end
                
                freqstat.label=freq.label;
                freqstat.cfg=cfg;
                freqstat.prob_pos=pos_prob;
                freqstat.prob_neg=neg_prob;
                freqstat.time_window= sel_window;
                freqstat.conditions=conditions;
                freqstat.sel_reg=sel_reg;
                freqstat.sel_cond=sel_cond;
                freqstat.sel_trials=sel_trials;
                freqstat.foi=foi;
                freqstat.f_name=f_name;
                
                % save in specific folder
                sel_path=fullfile(path_out,sel_cond,[sel_cond,f_name,num2str(foi(1)),'to',num2str(foi(2)),'Hz','_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec']);
                mkdir(sel_path)
                save(fullfile(sel_path,strcat(sel_sub,'freqstat')),'freqstat')
                clear freqstat pos_prob neg_prob stat
                
            end
        end
    end
end
%% get roi group stats

