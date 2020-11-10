% fix stats in rsa

% p-values were accidently halfed, this error is corrected in the stats
% script, here the already calculated stats are corrected

% new p cluster (twosided)
p_cluster=0.025;
path_stats='D:\Extinction\iEEG\analysis\rsa\powlogscale_timeslide_z_crosstrials_toi2000to4000\stats';

all_files=dir(fullfile(path_stats,'*\fig\*.mat'));


for f=1:numel(all_files)
    sel_path=all_files(f).folder;
    sel_file=all_files(f).name;
  load(fullfile(sel_path,sel_file))  
        
  if exist('stat_data')
      stats=stat_data;
      clear stat_data
  end
         rand_neg=stats.trial_rand.rand_neg;
         rand_pos=stats.trial_rand.rand_pos;
        nrand=numel(rand_neg);
         
        data_pos=0;
        data_neg=0;
        p_pos=1;
        p_neg=1;
        if isfield(stats,'posclusters')
            if ~isempty(stats.posclusters)
                data_pos=[stats.posclusters(:).clusterstat];
                for i=1:numel(data_pos)
                    p_pos(i)=(nearest(rand_pos,data_pos(i))./nrand);
                end
            end
        end
        if isfield(stats,'negclusters')
            if ~isempty(stats.negclusters)
                data_neg=[stats.negclusters(:).clusterstat];
                for i=1:numel(data_neg)
                    p_neg(i)=(nearest(rand_neg,data_neg(i))./nrand);
                end
            end
        end
  
  if isfield(stats,'posclusterslabelmat')
            maskpos=(stats.posclusterslabelmat<=sum(p_pos<p_cluster)&stats.posclusterslabelmat>0);
        else
            maskpos=zeros(size(stats.stat));
        end
        if isfield(stats,'negclusterslabelmat')
            maskneg=(stats.negclusterslabelmat<=sum(p_neg<p_cluster)&stats.negclusterslabelmat>0);
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
        clear maskneg maskpos mask p_pos p_neg rand_neg rand_neg data_pos data_neg
        save(fullfile(sel_path,sel_file),'stats')
        clear stats
end