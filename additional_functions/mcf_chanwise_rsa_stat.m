function stat=mcf_chanwise_rsa_stat(sel_chan,data,cfg_freq,cfg_rsa,norm)
   trlinfo=data.trlinfo;


  % define features
            cfg_freq.channel=sel_chan;
            freq=ft_freqanalysis(cfg_freq,data);
            clear data
    
    num_trials=size(trlinfo,1);
    num_freq=numel(freq.freq);
    num_time=numel(freq.time);
    num_chan=numel(freq.label);
    
    switch norm
        case 'z_crosstrials'
            m=nanmean(freq.powspctrm,1);
            s=nanstd(freq.powspctrm,1);
            freq.powspctrm=(freq.powspctrm-repmat(m,num_trials,1,1,1))./repmat(s,num_trials,1,1,1);
            clear m s
    end
    
    %%%%%% here: create time bin vectors
    % reshape in feature vec
    
    % tmp_vec: chan_trialXtimebin
    % average for each bin
    n_bins=cfg_rsa.n_bins;
    num_trl=size(freq.powspctrm,1);
    bin_mat=zeros(num_trl,numel(freq.label),numel(freq.freq),n_bins);
    for t=1:cfg_rsa.n_bins
        bin_mat(:,:,:,t)=nanmean(freq.powspctrm(:,:,:,cfg_rsa.ind_t1(t):cfg_rsa.ind_t2(t)),4);
    end
    tmp_vec=permute(bin_mat,[2,1,4,3]);% chan x trial x bin x freq
    tmp_vec=reshape(tmp_vec,numel(freq.label),[],numel(freq.freq));
    clear bin_mat
    
    rsa.rsa_mat=zeros(num_chan,num_trl,num_trl,n_bins,n_bins,'single');
    % loop over channels to avoid memory overload
    for chan=1:numel(freq.label)
        tmp=corr(squeeze(tmp_vec(chan,:,:))','Type','Spearman');
        % fisher z
        %rsa.rsa_mat=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
        tmp=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
        tmp=reshape(tmp,num_trl,n_bins,num_trl,n_bins);
        rsa.rsa_mat(chan,:,:,:,:)=single(permute(tmp,[1,3,2,4]));
        clear tmp
    end
    
    clear tmp_vec
    %rsa.dim='trial_trial_time_time';
    rsa.dim='chan_trial_trial_time_time';
    
    
    
    rsa.cfg.corr='Spearman';
    rsa.cfg.norm=norm;
    rsa.cfg.pow=cfg_freq;
    rsa.trlinfo=trlinfo;
    rsa.label=freq.label;
    rsa.t1=cfg_rsa.t1;
    rsa.t2=cfg_rsa.t2;
    rsa.time=(cfg_rsa.t1+cfg_rsa.t2)./2;
    clear freq
    % avoid saving because of file size
    %     sel_folder=fullfile(path_out,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
    %     mkdir(sel_folder)
    %     save(fullfile(sel_folder,strcat(sel_sub,'_',sel_chan,'_rsa_')),'rsa','-v7.3')
    %     clear rsa
    
    


 
        
        sortind=cfg_rsa.sortind;
        contrast_vec=reshape(cfg_rsa.contrast_mat,[],1);
        % select only non nan trials
        sel_ind=~isnan(contrast_vec);
        contrast_vec=contrast_vec(sel_ind);
        % sort rsa to match contrast_mat
        rsa.trlinfo=rsa.trlinfo(sortind,:);
        rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:,:);
        rsa.rsa_mat=rsa.rsa_mat(:,:,sortind,:,:);
        num_trial=size(rsa.trlinfo,1);
        num_bin=numel(rsa.t1);
        
        
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


        data_dummy.label={'stat_sel'};
             data_dummy.freq=rsa.time;
             data_dummy.time=rsa.time;
             data_dummy.dimord='rpt_chan_freq_time';    
             
            rsa_tmp=reshape(rsa.rsa_mat,num_chan,[],num_bin,num_bin);
            sel_rsa=single(rsa_tmp(:,sel_ind,:,:));
               label=rsa.label;
               % rotate/flip matrix to get all trial combination
               % (upper/lower diagonal)
            sel_rsa2=permute(rot90(flipud(permute(sel_rsa,[3,4,1,2])),-1),[3,4,1,2]);
            % design needs to account for upper and lower diagonal trials
             design = zeros(1,sum(contrast_vec==1).*2+ sum(contrast_vec==0).*2);
            design(1,1:sum(contrast_vec==1)*2) = 1;
            design(1,(sum(contrast_vec==1)*2+1):(sum(contrast_vec==1)*2 + sum(contrast_vec==0)*2))= 2;
            cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
            cfg.design = design;
            clear rsa_tmp contrast_mat  datainfo design rsa trlinfo
        for chan=1:num_chan
                             
            data1_mat3d=[squeeze(sel_rsa(chan,contrast_vec==1,:,:));squeeze(sel_rsa2(chan,contrast_vec==1,:,:))];
            data2_mat3d=[squeeze(sel_rsa(chan,contrast_vec==0,:,:));squeeze(sel_rsa2(chan,contrast_vec==0,:,:))];
            
            data1=data_dummy;
           data2=data_dummy;
             % put rsa_mat in ft freq structure
             data1.powspctrm=reshape(data1_mat3d,[size(data1_mat3d,1),1,size(data1_mat3d,2),size(data1_mat3d,3)]); 
             data2.powspctrm=reshape(data2_mat3d,[size(data2_mat3d,1),1,size(data2_mat3d,2),size(data2_mat3d,3)]);        
            stat=ft_freqstatistics (cfg,data1,data2)
            stat.label=label(chan);
             % subfunction for parfor loop 
           % all_stat{chan}=mcf_chanwiseclusterstat(cfg,data1_mat3d,data2_mat3d, data_dummy)
        end
            
    end

