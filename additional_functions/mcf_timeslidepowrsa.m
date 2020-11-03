  



function rsa=mcf_timeslidepowrsa(cfg_freq,cfg_rsa,data)

% define features
    freq=ft_freqanalysis(cfg_freq,data);
    clear data
    
    num_trials=size(trialinfo,1);
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
    tmp_vec=permute(bin_mat,[1,4,2,3]);% trial x bin x chan x freq
    tmp_vec=reshape(tmp_vec,[],numel(freq.label)*numel(freq.freq));
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
    