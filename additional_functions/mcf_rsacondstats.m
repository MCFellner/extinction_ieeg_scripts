% runs cluster based permutation stats

% cfg_stats.nrand=nrand;
% cfg_stats.permutation='yes';

function [stats]=mcf_rsacondstats(cfg_stats,rsa_ga)

data_dummy.label={sel_roi};
data_dummy.freq=rsa.time;
data_dummy.time=rsa.time;
data_dummy.dimord='subj_chan_freq_time';
data1=data_dummy;
data1.powspctrm=cond_rsa(:,1,:,:);

data2=data_dummy;
data2.powspctrm=cond_rsa(:,2,:,:);

% define freqstats
cfg=[];
cfg.avgoverchan =  'yes';
cfg.avgovertime =  'no';
cfg.avgoverfreq =  'no';

% first level
cfg.method           = 'montecarlo';
cfg.numrandomization = nrand;
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
stat_data=ft_freqstatistics (cfg,data1,data2);
clear design


% parfor?
parfor i=1:nrand
    % run permutation stats (using_ft_freqstats)
    Nsub = size(rand_rsa,1);
    data1=data_dummy;
    data1.powspctrm=reshape(squeeze(rand_rsa(:,1,i,:,:)),Nsub,1,n_bins,n_bins);
    
    data2=data_dummy;
    data2.powspctrm=reshape(squeeze(rand_rsa(:,2,i,:,:)),Nsub,1,n_bins,n_bins);
    
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
if isfield(stat_data,'posclusters')
    if ~isempty(stat_data.posclusters)
        data_pos=[stat_data.posclusters(:).clusterstat];
        for i=1:numel(data_pos)
            p_pos(i)=(nearest(rand_pos,data_pos(i))./nrand).*0.5;
        end
    end
end
if isfield(stat_data,'negclusters')
    if ~isempty(stat_data.negclusters)
        data_neg=[stat_data.negclusters(:).clusterstat];
        for i=1:numel(data_neg)
            p_neg(i)=(nearest(rand_neg,data_neg(i))./nrand).*0.5;
        end
    end
end


if isfield(stat_data,'posclusterslabelmat')
    maskpos=(stat_data.posclusterslabelmat<=sum(p_pos<0.05)&stat_data.posclusterslabelmat>0);
else
    maskpos=zeros(size(stat_data.stat));
end
if isfield(stat_data,'negclusterslabelmat')
    maskneg=(stat_data.negclusterslabelmat<=sum(p_neg<0.05)&stat_data.negclusterslabelmat>0);
else
    maskneg=zeros(size(stat_data.stat));
end

mask=maskpos+maskneg;

fig=figure
imagesc(stat_data.time,stat_data.time,squeeze(stat_data.stat),[-5 5])
hold on
colormap(jet_grey)
colorbar
title({[contrast,' in ',sel_roi];['pos tsum:',num2str(data_pos(1)),'p=',num2str(p_pos(1))];['neg tsum:',num2str(data_neg(1)),'p=',num2str(p_neg(1))]})
ylabel('t in s')
xlabel('t in s')
contour(stat_data.time,stat_data.time,squeeze(mask),1,'k')
set(gca,'YDir','normal')
% plot results
path_fig=fullfile( folder_out,'fig');
savefig(fig,[path_fig,'\',contrast,'_in_',sel_roi,],'compact')

stat_data.trial_rand.p_pos=p_pos;
stat_data.trial_rand.p_neg=p_neg;

stat_data.trial_rand.rand_pos=rand_pos;
stat_data.trial_rand.rand_neg=rand_neg;
save([path_fig,'\',contrast,'_in_',sel_roi,'.mat'],'stat_data')
close all
clear stat_data stat_rand rand_pos rand_neg p_pos p_neg data_neg data_pos
