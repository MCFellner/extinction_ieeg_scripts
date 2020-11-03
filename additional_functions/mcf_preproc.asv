% standard preproc

% cfg=[]
%             cfg.montage=datainfo.elec_info.bipolar.montage_withoutartichan;
%             cfg.resamplefs      = sr;
%             cfg.trl_def=trl;
%             cfg.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
%             cfg.channel     = all_roi_labels{sub_ind,r};

function data=mcf_preproc(cfg_all, data)
            
            
            % apply bipolar montage
            montage=cfg_all.montage;
            data = ft_apply_montage(data,montage);
            clear montage
            
            % cut trials

            cfg=[];
            cfg.trl=cfg_all.trl_def;
            data=ft_redefinetrial(cfg,data);
            clear trl
            data.trialinfo=cfg_all.trialinfo;
            % downsample to common sampling rate
            cfg=[];
            cfg.resamplefs      = cfg_all.resamplefs;
            cfg.detrend='yes';
            data=ft_resampledata(cfg,data);
            
            % only select artifree trials in trialinfo

            cfg=[]
            cfg.trials=cfg_all.trials;
            cfg.channel     = cfg_all.channel;
            data=ft_selectdata(cfg,data);
            
