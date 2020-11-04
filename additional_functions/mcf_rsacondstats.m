% runs cluster based permutation stats

% cfg_stats.nrand=nrand;
% cfg_stats.permutation='yes';

function [stats]=mcf_rsacondstats(cfg_stats,rsa_ga)
n_bins=numel(rsa_ga.time);
nrand=cfg_stats.nrand;

data_dummy.label={rsa_ga.roi};
data_dummy.freq=rsa_ga.time;
data_dummy.time=rsa_ga.time;
data_dummy.dimord='subj_chan_freq_time';
data1=data_dummy;
data1.powspctrm=rsa_ga.cond_rsa(:,1,:,:);

data2=data_dummy;
data2.powspctrm=rsa_ga.cond_rsa(:,2,:,:);

% define freqstats
cfg=[];
cfg.avgoverchan =  'yes';
cfg.avgovertime =  'no';
cfg.avgoverfreq =  'no';

% first level
cfg.method           = 'montecarlo';
cfg.numrandomization = cfg_stats.nrand;
cfg.correctm         =  'cluster';
cfg.correcttail      = 'prob';
cfg.statistic ='depsamplesT'
% for within-subjects (depsamplesT)
Nsub = size(data1.powspctrm,1);                                       %# of subjects?
design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

cfg.uvar     = 2;
cfg.ivar     = 1;
cfg.design = design;
stats=ft_freqstatistics (cfg,data1,data2);
clear design

switch cfg_stats.permutation
    case 'yes'
        % parfor?
        for i=1:nrand
            % run permutation stats (using_ft_freqstats)
            Nsub = size(rsa_ga.rand_rsa,1);
            data1=data_dummy;
            data1.powspctrm=reshape(squeeze(rsa_ga.rand_rsa(:,1,i,:,:)),Nsub,1,n_bins,n_bins);
            
            data2=data_dummy;
            data2.powspctrm=reshape(squeeze(rsa_ga.rand_rsa(:,2,i,:,:)),Nsub,1,n_bins,n_bins);
            
            cfg=[];
            cfg.avgoverchan =  'yes';
            cfg.avgovertime =  'no';
            cfg.avgoverfreq =  'no';
            
            % first level
            cfg.method           = 'montecarlo';
            cfg.numrandomization = 2;
            cfg.correctm         =  'cluster';
            cfg.correcttail      = 'prob';
            cfg.statistic ='depsamplesT'
            % for within-subjects (depsamplesT)
            Nsub = size(data1.powspctrm,1);
            %# of subjects?
            design=[];
            design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
            design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
            
            cfg.uvar     = 2;
            cfg.ivar     = 1;
            cfg.design = design;
            stat_rand=ft_freqstatistics (cfg,data1,data2);
            % get pos/neg clusters for each random
            rand_pos(i)=0;
            rand_neg(i)=0;
            if isfield(stat_rand,'posclusters')
                if ~isempty(stat_rand.posclusters)
                    rand_pos(i)=stat_rand.posclusters(1).clusterstat;
                end
            end
            if isfield(stat_rand,'negclusters')
                if ~isempty(stat_rand.negclusters)
                    rand_neg(i)=stat_rand.negclusters(1).clusterstat;
                end
            end
        end
        
        
        % check significance
        
        rand_neg=sort(rand_neg,'ascend');
        rand_pos=sort(rand_pos,'descend');
        
        data_pos=0;
        data_neg=0;
        p_pos=1;
        p_neg=1;
        if isfield(stats,'posclusters')
            if ~isempty(stats.posclusters)
                data_pos=[stats.posclusters(:).clusterstat];
                for i=1:numel(data_pos)
                    p_pos(i)=(nearest(rand_pos,data_pos(i))./nrand).*0.5;
                end
            end
        end
        if isfield(stats,'negclusters')
            if ~isempty(stats.negclusters)
                data_neg=[stats.negclusters(:).clusterstat];
                for i=1:numel(data_neg)
                    p_neg(i)=(nearest(rand_neg,data_neg(i))./nrand).*0.5;
                end
            end
        end
        
        
        if isfield(stats,'posclusterslabelmat')
            maskpos=(stats.posclusterslabelmat<=sum(p_pos<0.05)&stats.posclusterslabelmat>0);
        else
            maskpos=zeros(size(stats.stat));
        end
        if isfield(stats,'negclusterslabelmat')
            maskneg=(stats.negclusterslabelmat<=sum(p_neg<0.05)&stats.negclusterslabelmat>0);
        else
            maskneg=zeros(size(stats.stat));
        end
        
        mask=maskpos+maskneg;
        
        
        stats.trial_rand.mask=mask;
        stats.trial_rand.p_pos=p_pos;
        stats.trial_rand.p_neg=p_neg;
        stats.trial_rand.rand_pos=rand_pos;
        stats.trial_rand.rand_neg=rand_neg;
        stats.trial_rand.data_pos=data_pos;
        stats.trial_rand.data_neg=data_neg;
    otherwise
end
