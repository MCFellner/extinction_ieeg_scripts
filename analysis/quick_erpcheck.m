

% quick and dirty amydala erp for one patient
addpath('D:\matlab_tools\fieldtrip-20200130')
addpath('D:\matlab_tools\)
ft_defaults
'
% load data
load('D:\Extinction\iEEG\data\preproc\ieeg\readin\p_sub03_data.mat')
load('D:\Extinction\iEEG\data\preproc\ieeg\datainfo\p_sub03_datainfo.mat')


% check atlas label for localization
amy_channel=37;



cfg=[];
cfg.channel=datainfo.elec_info.ana_labels.labels{amy_channel};
cfg.lpfilter='yes';
cfg.lpfreq=15;
sel_data=ft_preprocessing(cfg,data)

% The trial definition "trl" is an Nx3 matrix, N is the number of trials.
% The first column contains the sample-indices of the begin of each trial
% relative to the begin of the raw data, the second column contains the
% sample-indices of the end of each trial, and the third column contains
% the offset of the trigger with respect to the trial. An offset of 0
% means that the first sample of the trial corresponds to the trigger. A
% positive offset indicates that the first sample is later than the trigger,
% a negative offset indicates that the trial begins before the trigger.

trl(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*0.5);
trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*5);
trl(:,3)=-(data.fsample.*0.5);
cfg=[];
cfg.trl=trl;
sel_data=ft_redefinetrial(cfg,sel_data);


sel_trials1=datainfo.trialinfo(:,9)==1;
sel_trials2=datainfo.trialinfo(:,9)==0&datainfo.trialinfo(:,2)<=2;

cfg=[];
cfg.trials=find(sel_trials1);
cfg.keeptrials='yes';
erp1=ft_timelockanalysis(cfg,sel_data);
cfg=[];
cfg.keeptrials='yes';
cfg.trials=find(sel_trials2);
erp2=ft_timelockanalysis(cfg,sel_data);
cfg=[];
cfg.baseline=[-0.2 0];
erp1=ft_timelockbaseline(cfg,erp1)
erp2=ft_timelockbaseline(cfg,erp2)



[h,p,~,stat]=ttest2(squeeze(erp1.trial),squeeze(erp2.trial))
[pthr,pcor,padj]=fdr(p);
sig_sample=find(pcor<0.05);
figure
plot(padj)

figure
hold on
for i=1:numel(sig_sample)
plot([erp1.time(sig_sample(i)),erp1.time(sig_sample(i))],[-60 60],'Color',[0.9 0.9 0.9],'LineWidth',2);
end

plot(erp1.time, squeeze(nanmean(erp1.trial)),'r')
plot(erp2.time, squeeze(nanmean(erp2.trial)),'g')
legend('us trial', 'no us trial')
xlabel('time relative to video onset')
ylabel('amplitude')

