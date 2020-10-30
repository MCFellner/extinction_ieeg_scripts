% automatic artifact marking in continuous data

% do automatic checks
% chop data in arbitrary trials (1 sec, fake trigger every 0.8 s, pre 0.1, post 0.9)

% metrics: same as rejectvisual

% cfg.segmentlength=1; in sec
% cfg.segmentoverlap=0.2, in sec
% cfg.zcutoff=3;
% cfg.iteration=number or until 'no_rejected' values% how many iteration of rejecting trial with z>cutoff
% cfg.metric={ 'var','maxabs','range','kurtosis','zvalue','var_over_one'};
% cfg.hpfilter='yes';
% cfg.hpfreq=0.1;

% select epochs in channel
% or select epochs for all channels


% output: struct with sample points for segment exceeding z cutoff for
% every metric

function [auto_artifact]=mcf_continuous_autoarti(cfg, data)
slide=round(cfg.segmentoverlap.*data.fsample);
trig_int=round(cfg.segmentlength.*data.fsample-slide);
seg=round(cfg.segmentlength.*data.fsample);
metric=cfg.metric;
zcutoff=cfg.zcutoff;
iteration=cfg.iteration;

auto_artifact.cfg=cfg;

sp_1=data.sampleinfo(1);
sp_end=data.sampleinfo(2);
all_endsp=sp_1+seg:trig_int:sp_end;
trl(:,1)=all_endsp-seg;
trl(:,2)=all_endsp;
trl(:,3)=ones(size(trl,1),1).*round(double(slide).*-0.5);
clear sp_1 sp_end all_endsp

% run the hp filter
cfg.hpinstabilityfix ='reduce';
data=ft_preprocessing(cfg,data);

% chop data in fake trials
cfg=[];
cfg.trl=trl;
datain=ft_redefinetrial(cfg,data)
dat=vertcat(datain.trial{:}); % (chan*trial) x time

if strcmp(iteration, 'no_rejected')
    threshold=-1;
else
    threshold=iteration;
end


for m=1:numel(metric)
    sel_metric=metric{m};
    counter=0;
    nan_trial=nan(size(dat,2),1);
    dat=vertcat(datain.trial{:}); % (chan*trial) x time
    z_threshold=zeros(1,numel(datain.label),numel(datain.trial));
    while counter~=threshold
        counter=counter+1;
        % calculate metrics (value for each channel x trial)
        switch sel_metric
            case 'var'
                metric_result = nanstd(dat, [], 2).^2;
            case 'maxabs'
                metric_result = nanmax(abs(dat), [], 2);
            case 'range'
                metric_result = nanmax(dat, [], 2) - nanmin(dat, [], 2);
            case 'kurtosis'
                metric_result = kurtosis(dat, [], 2);
            case 'zvalue'
                % for z value
                tmp_dat=reshape(dat, numel(datain.label),[]);
                mval=nanmean(tmp_dat,2);
                sd=nanstd(tmp_dat,0,2);
                metric_result = nanmean((dat-repmat(mval, numel(datain.trial),size(dat,2) ))./repmat(sd, numel(datain.trial),size(dat,2)),2);
            case 'var_over_one'
                metric_result = 1./(nanstd(dat, [], 2).^2);

            otherwise
                sel_metric
        end
        metric_result = reshape(metric_result, numel(datain.label),[]); % format chan*trial
        % nan function replace data with 0, set all arti trial nan
        mask=find(sum(z_threshold));
        metric_result(mask)=nan;
        % remove all trial with z value of metric higher than zcutoff
         zmap_tmp=nanzscore(metric_result,0,2);
        z_threshold_tmp=(zmap_tmp)>zcutoff;
        z_threshold(counter,:,:)=z_threshold_tmp;
        mask=find(sum(z_threshold));
        
        % replace above threshold trials with nan trial
        nan_ind=reshape(z_threshold_tmp,[],1);
        sum(nan_ind)
        dat(nan_ind,:)=repmat(nan_trial',sum(nan_ind),1);
        tmp_artis.iteration=counter;
        if strcmp(iteration, 'no_rejected')& sum(sum(z_threshold_tmp))==0
            counter=-1;
        end
    end
    zmap_tmp=nanzscore(metric_result);
    z_threshold_tmp=(zmap_tmp)>zcutoff;
    tmp_artis.arti_z_across_chan=z_threshold_tmp;
    tmp_artis.arti_pertrial_chan=squeeze(sum(z_threshold));
    auto_artifact = setfield(auto_artifact,sel_metric,tmp_artis)
    clear z_threshold z_threshold_tmp zmap_tmp  metric_result
end
auto_artifact.labels=datain.label;
auto_artifact.samples=trl;