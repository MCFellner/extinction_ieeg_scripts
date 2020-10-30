addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% stats for timeslide rsa

%%%%%%%%%%%%%%%%%old, integrated in timeslide for memory issues
% select contrast matrix

% define electrodes



% load rsa mat
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_rsa='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

toi=[2 4];
feature='powlogscale';
norm='z_crosstrials';
nrand=1000;

contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
    'cs_specific','cs_specific_block1','cs_specific_block2',...
    'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};

for c=1%numel(contrasts)
    contrast=contrasts{c};
    folder_out=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
    mkdir(folder_out)
    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        sel_folder=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
        
        
        load(fullfile(path_designmat,strcat(sel_sub,'_contrast_mat')))
        
        contrast_mat=getfield(contrast_def,contrast);
        sortind=contrast_def.sortind_org2usedtrlinfo;
        contrast_vec=reshape(contrast_mat,[],1);
        % select only non nan trials
        sel_ind=~isnan(contrast_vec);
        contrast_vec=contrast_vec(sel_ind);
        
        
        all_chan=dir([sel_folder,'\*rsa_.mat']);
        all_chan={all_chan(:).name}';
        all_chan = erase(all_chan,sel_sub);
        all_chan= erase(all_chan,'_');
        all_chan = erase(all_chan,'rsa.mat');
       for chan=1:numel(all_chan)
            sel_chan=all_chan{chan};
            load(fullfile(sel_folder,strcat(sel_sub,'_',sel_chan,'_rsa_')))
            
            % sort rsa to match contrast_mat
            rsa.trlinfo=rsa.trlinfo(sortind,:);
            rsa.rsa_mat=rsa.rsa_mat(sortind,:,:,:);
            rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:);
            
            
            num_trial=size(rsa.trlinfo,1);
            num_bin=numel(rsa.t1);
            
            
            rsa_tmp=reshape(rsa.rsa_mat,[],num_bin,num_bin);
            rsa_tmp=rsa_tmp(sel_ind,:,:);
            
            % put rsa_mat in ft freq structure
            rsa1_test.powspctrm(:,1,:,:)=rsa_tmp(contrast_vec==1,:,:);
            rsa1_test.label={rsa.label};
            rsa1_test.freq=rsa.time;
            rsa1_test.time=rsa.time;
            rsa1_test.dimord='rpt_chan_freq_time';
            
            rsa2_test=rsa1_test;
            rsa2_test=rmfield(rsa2_test,'powspctrm');
            rsa2_test.powspctrm(:,1,:,:)=rsa_tmp(contrast_vec==0,:,:);
            
            
            cfg=[];
            cfg.latency     = [rsa.time(1),rsa.time(end)];
            cfg.avgoverchan =  'yes';
            cfg.avgovertime =  'no';
            cfg.avgoverfreq =  'no';
            
            % first level
            cfg.method           = 'montecarlo';
            cfg.numrandomization = 1000;
            cfg.correctm         =  'cluster';
            cfg.correcttail      = 'prob';
            cfg.statistic ='indepsamplesT'
            design = zeros(1,size(rsa1_test.powspctrm,1) + size(rsa2_test.powspctrm,1));
            design(1,1:size(rsa1_test.powspctrm,1)) = 1;
            design(1,(size(rsa1_test.powspctrm,1)+1):(size(rsa1_test.powspctrm,1) + size(rsa2_test.powspctrm,1)))= 2;
            cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
            cfg.design = design;
            stat=ft_freqstatistics (cfg,rsa1_test,rsa2_test)
            all_stat{chan}=stat;
            
            clear stat rsa1_test rsa2_test
        end
        save(fullfile(folder_out,strcat(sel_sub,'_rsastat')),'all_stat')
        clear all_stat contrast_mat contrast_vec design nrand rsa rsa1_test sortin rsa_tmp
    end
end