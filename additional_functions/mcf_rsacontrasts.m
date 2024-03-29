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
% select only non nan trials
sel_ind=~isnan(contrast_vec);
contrast_vec=contrast_vec(sel_ind);
rsa.contrast_vec=contrast_vec;

switch rsa.dim
    case 'trial_trial_time_time'
        % sort rsa to match contrast_mat
        %     rsa.trialinfo=rsa.trialinfo(sortind,:);
        n_bins=numel(rsa.time);
        
        rsa.rsa_mat=rsa.rsa_mat(sortind,:,:,:);
        rsa.rsa_mat=rsa.rsa_mat(:,sortind,:,:);
        
        % select only non nan trials
        rsa_mat=double(reshape(rsa.rsa_mat,[],n_bins,n_bins));
        
        rsa_mat=rsa_mat(sel_ind,:,:);
        %%%%%%%%%%%%
        
        % combine data to real avg and random avg
        
        % allocate
        num_cond=max(contrast_vec);
        cond_rsa=zeros(num_cond,size(rsa_mat,2),size(rsa_mat,3));
        for c=1:num_cond
            cond_rsa(c,:,:)=nanmean(rsa_mat(rsa.contrast_vec==c,:,:));
        end
        
        if num_cond==4
            tmp(1,:,:)=cond_rsa(1,:,:)-cond_rsa(2,:,:);
            tmp(2,:,:)=cond_rsa(3,:,:)-cond_rsa(4,:,:);
            cond_rsa=tmp;
            clear tmp
        else
        end
        
        switch cfg_con.generate_rand
            case 'yes'
                nrand=cfg_con.nrand;
                
                
                [~,rand_ind]=sort(rand(numel(rsa.contrast_vec),nrand),1);
                rand_vec=repmat(rsa.contrast_vec,1,nrand);
                rand_vec=rand_vec(rand_ind);
                rand_vec=rand_vec';
                clear rand_ind
                rsa_tmp=reshape(rsa_mat,size(rsa_mat,1),[]);
                rand_rsa=zeros(num_cond,nrand,size(rsa_tmp,2));
                for c=1:num_cond
                    
                    rand_tmp=rand_vec==c;
                    sum_c=sum(rand_tmp(1,:));
                    rand_rsa2(c,:,:)=(rand_tmp*rsa_tmp)./sum_c;
                    
                end
                rand_rsa=reshape(rand_rsa,num_cond,nrand,n_bins,n_bins);
                
                % old slower code
                %         rand_rsa=zeros(num_cond,nrand,size(rsa.rsa_mat,2),size(rsa.rsa_mat,3));
                %         [~,rand_ind]=sort(rand(numel(rsa.contrast_vec),nrand),1);
                %         rand_vec=repmat(rsa.contrast_vec,1,nrand);
                %         rand_vec=rand_vec(rand_ind);
                %         clear rand_ind
                %         for i=1:nrand
                %             for c=1:num_cond
                %                 rand_rsa(c,i,:,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==c,:,:));
                %             end
                %         end
                %         clear rand_vec
                
                
                if num_cond==4
                    tmp(1,:,:,:)=rand_rsa(1,:,:,:)-rand_rsa(2,:,:,:);
                    tmp(2,:,:,:)=rand_rsa(3,:,:,:)-rand_rsa(4,:,:,:);
                    rand_rsa=tmp;
                else
                end
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
        
    case 'freq_trial_trial'
        n_freq=numel(rsa.freq);
        rsa.rsa_mat=rsa.rsa_mat(:,sortind,:);
        rsa.rsa_mat=rsa.rsa_mat(:,:,sortind);
        rsa.rsa_mat=permute(rsa.rsa_mat,[2,3,1]);
        rsa.rsa_mat=reshape(rsa.rsa_mat,[],n_freq);
        
        % select only non nan trials
        rsa.rsa_mat=rsa.rsa_mat(sel_ind,:);
        
        
        % allocate
        num_cond=max(contrast_vec);
        cond_rsa=zeros(num_cond,size(rsa.rsa_mat,2));
        for c=1:num_cond
            cond_rsa(c,:)=nanmean(rsa.rsa_mat(rsa.contrast_vec==c,:,:));
        end
        
        if num_cond==4
            tmp(1,:,:)=cond_rsa(1,:,:)-cond_rsa(2,:,:);
            tmp(2,:,:)=cond_rsa(3,:,:)-cond_rsa(4,:,:);
            cond_rsa=tmp;
            clear tmp
        else
        end
        
        switch cfg_con.generate_rand
            case 'yes'
                nrand=cfg_con.nrand;
                
                rand_rsa=zeros(num_cond,nrand,size(rsa.rsa_mat,2));
                [~,rand_ind]=sort(rand(numel(rsa.contrast_vec),nrand),1);
                rand_vec=repmat(rsa.contrast_vec,1,nrand);
                rand_vec=rand_vec(rand_ind);
                clear rand_ind
                for i=1:nrand
                    for c=1:num_cond
                        rand_rsa(c,i,:)=nanmean(rsa.rsa_mat(rand_vec(:,i)==c,:));
                    end
                end
                clear rand_vec
                if num_cond==4
                    tmp(1,:,:)=rand_rsa(1,:,:)-rand_rsa(2,:,:);
                    tmp(2,:,:)=rand_rsa(3,:,:)-rand_rsa(4,:,:);
                    rand_rsa=tmp;
                else
                end
                rsa_cond.rand_rsa=rand_rsa;
            otherwise
        end
        rsa_cond.cond_rsa=cond_rsa;
        rsa_cond.freq=rsa.freq;
        rsa_cond.dim_cond='cond_freq';
        rsa_cond.dim_rand='cond_rand_freq';
        rsa_cond.cfg=cfg_con;
        rsa_cond.roi=rsa.roi;
        
        
        
        
end
