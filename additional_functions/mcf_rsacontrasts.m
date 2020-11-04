% sub function that uses the contrast_mat as input to calculate rsa
% averages (also randomized)

%     cfg_con.generate_rand='yes';
%     cfg_con.nrand=nrand;
%     cfg_con.contrast_mat=contrast_mat;
%     cfg_con.sortind=contrast_def.sortind_org2usedtrlinfo;


function [rsa_cond]=mcf_rsacontrasts(cfg_con,rsa);
cfg_con.trialinfo=rsa.trialinfo;
sortind=cfg_con.sortind;
contrast_vec=reshape(cfg_con.contrast_mat,[],1);
n_bins=numel(rsa.time);

% sort rsa to match contrast_mat
%     rsa.trialinfo=rsa.trialinfo(sortind,:);
rsa.rsa_mat=rsa.rsa_mat(sortind,:,:,:);
rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:);
rsa.rsa_mat=reshape(rsa.rsa_mat,[],n_bins,n_bins);

% select only non nan trials
sel_ind=~isnan(contrast_vec);
contrast_vec=contrast_vec(sel_ind);
rsa.contrast_vec=contrast_vec;
rsa.rsa_mat=rsa.rsa_mat(sel_ind,:,:);
%%%%%%%%%%%%

% combine data to real avg and random avg

% allocate
cond_rsa=zeros(2,size(rsa.rsa_mat,2),size(rsa.rsa_mat,3));
cond_rsa(1,:,:)=nanmean(rsa.rsa_mat(rsa.contrast_vec==1,:,:));
cond_rsa(2,:,:)=nanmean(rsa.rsa_mat(rsa.contrast_vec==0,:,:));

switch cfg_con.generate_rand
    case 'yes'
        nrand=cfg_con.nrand;
        
        rand_rsa=zeros(2,nrand,size(rsa.rsa_mat,2),size(rsa.rsa_mat,3));
        [~,rand_ind]=sort(rand(numel(rsa.contrast_vec),nrand),1);
        rand_vec=repmat(rsa.contrast_vec,1,nrand);
        rand_vec=rand_vec(rand_ind);
        clear rand_ind
        for i=1:nrand
            rand_rsa(1,i,:,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==1,:,:));
            rand_rsa(2,i,:,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==0,:,:));
        end
        clear rand_vec
        rsa_cond.rand_rsa=rand_rsa;
    otherwise
end
rsa_cond.cond_rsa=cond_rsa;
rsa_cond.time=rsa.time;
rsa_cond.t1=rsa.t1;
rsa_cond.t2=rsa.t2;
rsa_cond.dim_cond='cond_time_time';
rsa_cond.dim_rand='cond_rand_time_time';
rsa_cond.cfg=cfg_con;
rsa_cond.roi=rsa.roi;

