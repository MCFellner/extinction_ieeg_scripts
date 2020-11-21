 % selects defined clusters form rsa timextimextrialxtrial mat

 % stat cluster definition
% cfg.type='stat';
% cfg.stat_mask=stats.trial_rand.mask; % has dimension tx x ty
% cfg.stat_tx=stats.time; % define time
% cfg.stat_ty=stats.freq; % define time to ensure correct selection
% timewindow
% cfg.type='window';
% cfg.tx=[2 3]; % defined in sec start/end
% cfg.ty=[2 3];
            

function   rsa=    mcf_rsadataselect(cfg,rsa)


switch cfg.type
    case 'stat'
        % check if timing fits
        if ~strcmp(rsa.dim,'trial_trial_time_time')
        error('wrong data dim')
        end
        
       if any(cfg.stat_ty~=rsa.time)
        error(' time in stat mask and rsa do not match')
       elseif any(cfg.stat_tx~=rsa.time)
        error(' time in stat mask and rsa do not match')
       end
        mask(1,1,:,:)=cfg.stat_mask;
        mask=repmat(mask,size(rsa.rsa_mat,1),size(rsa.rsa_mat,1),1,1);
        mask(mask==0)=NaN;
        % select data
        rsa_tmp=squeeze(nanmean(nanmean(rsa.rsa_mat.*mask,3),4));
        rsa.dim='trial_trial';
        rsa.rsa_mat=rsa_tmp;
        cfg.previous=rsa.cfg;
        cfg.previous.time=rsa.time;
        cfg.previous.t1=rsa.t1;
        cfg.previous.t2=rsa.t2;

        rsa.cfg=cfg;
        rsa=rmfield(rsa,'t1');
        rsa=rmfield(rsa,'t2');
        rsa=rmfield(rsa,'time');

    case 'window'
        error('not implemented yet')
end

            