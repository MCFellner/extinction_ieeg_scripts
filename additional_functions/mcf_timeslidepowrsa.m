% calculate timeslide rsa using freq and spatial features (concatenates
% across all channels in data fed in, if num chan=1  results are channel
% specific)


   
 %   cfg_rsa.freqdef=pow_feature; def: 'powlogscale' sets definition of ft_freqanalysis
 %   cfg_rsa.norm=norm; def: 'z_crosstrials'normalization across trials?
%cfg_rsa.toi=[2 4]; % timewindow for which rsa is calculated
%cfg_rsa.win_pow=0.05; % in sec, power estimated for every win
%cfg_rsa.win=0.200; % length of sliding window
%cfg_rsa.slide=0.05;  % slide


function rsa=mcf_timeslidepowrsa(cfg_rsa,data)
feature=cfg_rsa.freqdef;   
toi=cfg_rsa.toi;
   win_pow= cfg_rsa.win_pow;
   win= cfg_rsa.win;
slide=cfg_rsa.slide;

    sr_pow=1/(win_pow);
    step=win/(1/sr_pow)+1;
    tois=toi(1):1/sr_pow:toi(2);
    t1=toi(1):slide:(toi(2)-win);
    t2=t1+win;
    ind_t1=1:slide/(1/sr_pow):((numel(tois)-win/(1/sr_pow)));
    ind_t2=ind_t1+win/(1/sr_pow);
    n_bins=numel(t1);

% define features

 switch feature
        case 'powlogscale'
            cfg_freq=[];
            cfg_freq.output='pow';
            %cfg.toi =pre_item:0.1:post_item;
            cfg_freq.keeptrials  = 'yes';
            cfg_freq.pad='nextpow2';
            cfg_freq.foi     = logspace(log10(2),log10(200), 50);
            cfg_freq.toi=toi(1):1/(1/(win_pow)):toi(2);
            cfg_freq.method='wavelet';
            cfg_freq.width = 5;
 end
    freq=ft_freqanalysis(cfg_freq,data);
    clear data
    
    num_trl=size(freq.trialinfo,1);
    num_freq=numel(freq.freq);
    num_chan=numel(freq.label);
    switch cfg_rsa.norm
        case 'z_crosstrials'
            m=nanmean(freq.powspctrm,1);
            s=nanstd(freq.powspctrm,1);
            freq.powspctrm=(freq.powspctrm-repmat(m,num_trl,1,1,1))./repmat(s,num_trl,1,1,1);
            clear m s
    end
    
    %%%%%% here: create time bin vectors
    % reshape in feature vec
        
    
    % tmp_vec: chan_trialXtimebin
    % average for each bin
    bin_mat=zeros(num_trl,num_chan,num_freq,n_bins);
    for t=1:n_bins
        bin_mat(:,:,:,t)=nanmean(freq.powspctrm(:,:,:,ind_t1(t):ind_t2(t)),4);
    end
    tmp_vec=permute(bin_mat,[1,4,2,3]);% trial x bin x chan x freq
    tmp_vec=reshape(tmp_vec,[],num_chan*num_freq);
    clear bin_mat
    
    rsa.rsa_mat=zeros(num_trl,num_trl,n_bins,n_bins,'single');
    tmp=corr(tmp_vec(:,:)','Type','Spearman');
    % fisher z
    %rsa.rsa_mat=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
    tmp=0.5.*log(((ones(size(tmp))+tmp)./(ones(size(tmp))-tmp)));
    tmp=reshape(tmp,num_trl,n_bins,num_trl,n_bins);
    rsa.rsa_mat(:,:,:,:)=single(permute(tmp,[1,3,2,4]));
    clear tmp
    clear tmp_vec
    
    % rsa_mat is symmetric (timextime), set half of the matrix to Nan to
    % avoid spurious large clusters in cluster stat
    nan_mat=repmat(reshape(triu(ones(n_bins),1),1,1,n_bins,n_bins),num_trl,num_trl,1,1);
    rsa.rsa_mat(logical(nan_mat))=NaN;
    %rsa.dim='trial_trial_time_time';
    rsa.dim='trial_time_time';
    rsa.trialinfo=freq.trialinfo;
    %%%%%%%%%% define trial dim properly!!!!!!!!!!!!
    %rsa.trldim=''
    