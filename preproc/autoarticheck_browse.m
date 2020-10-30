%% autoartifact check and browse

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_artifactinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\auto_arti\';
mkdir(path_artifactinfo)
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

 sel_sub='p_sub08';
 
 sel_ref='bip'; % or 'org'
 all_metric={ 'range','kurtosis'};
 
 
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    auto_arti_file=strcat(path_artifactinfo,sel_sub,'_autoarti_segment100ms')
    load(auto_arti_file)   
    load(strcat(path_preproc,sel_sub,'_data.mat'))

      % select ieeg channels & artifactfree channels
    firstchannelcheck=datainfo.artifact_info.firstchannelcheck;
    all_elec=[];
    all_out=[];
    for i=1:numel(firstchannelcheck.elec)
        sel_elec=firstchannelcheck.labels{i};
        out_elec_ind=[firstchannelcheck.out_electrode{i},firstchannelcheck.artifact_electrode{i}]
        out_elec=sel_elec(out_elec_ind);
        sel_elec(out_elec_ind)=[];
        all_elec=[all_elec;sel_elec];
        all_out=[all_out;out_elec];
    end
    
    cfg=[];
    cfg.channel= all_elec;
    data=ft_preprocessing(cfg,data);
    
    % bipolar referencing
    montage=datainfo.elec_info.bipolar.montage;
    for e=1:numel(all_out)
      sel_elec=all_out{e};
        ind_old=strcmp(montage.labelold,sel_elec);
        montage.tra(:,ind_old)=[];
        montage.labelold(ind_old)=[];
    end
    labelnew_ind=sum(montage.tra~=0,2)==2;
    montage.labelnew=montage.labelnew(labelnew_ind);
    montage.tra=montage.tra(labelnew_ind,:);
    data_bip = ft_apply_montage(data,montage);
    
    switch sel_ref
        case 'org'
            sel_auto_artifact=(auto_artifact.dataref);
            ref='org';
        case 'bip'
            sel_auto_artifact=(auto_artifact.databip);
            ref='bip';
            data=data_bip;
        end
        % evaluate detected artifacts
        for m=1:numel(all_metric)
            sel_metric=all_metric{m};
            detected_artis=getfield( sel_auto_artifact,sel_metric);
            artis=detected_artis.arti_pertrial_chan;
            % relative number of artis
            rel_chan= sum(artis,2)./size(artis,2);
            rel_trial=sum(artis,1)./size(artis,1);
            figure
            subplot(3,2,3)
            imagesc(artis)
            xlabel('segments')
            ylabel('channel')
            title(strcat(ref,': ',sel_metric))
            subplot(3,2,1)
            plot(rel_trial)
            subplot(3,2,4)
            plot(rel_chan,1:size(artis,1))
            subplot(3,2,5)
            hist(rel_trial)
            title('distribution trial')
             subplot(3,2,6)
            hist(rel_chan)
            title('distribution chan')           
            % relativ for electrodes
            % get measures for each electrode
            artis_all(m,:,:)=artis;
            arti_chan_all(m,:,:)=detected_artis.arti_z_across_chan;

        end
         artis=squeeze(sum(artis_all))>0;
         artis_chan=squeeze(sum(arti_chan_all))>0;
            % relative number of artis
            rel_chan= sum(artis,2)./size(artis,2);
            rel_trial=sum(artis,1)./size(artis,1);
            figure
            subplot(3,2,3)
            imagesc(artis)
            xlabel('segments')
            ylabel('channel')
            title(strcat(ref,': combined'))
            subplot(3,2,1)
            plot(rel_trial)
            subplot(3,2,4)
            plot(rel_chan,1:size(artis,1))
            subplot(3,2,5)
            hist(rel_trial)
            title('distribution trial')
             subplot(3,2,6)
            hist(rel_chan)
            title('distribution chan')           
    
%sugest channels to remove
z_chan=zscore(sum(artis_chan,2));
rel_arti_chan=sum(artis_chan,2)./size(artis_chan,2);
reject_channel=[sel_auto_artifact.labels(z_chan>3)',num2cell(rel_arti_chan(z_chan>3))]

% check for dead channels/low variance
z_chan_novar=nanzscore(sum(sel_auto_artifact.var_over_one.arti_z_across_chan,2));
rel_novar_chan=sum(sel_auto_artifact.var_over_one.arti_z_across_chan,2)./size(sel_auto_artifact.var_over_one.arti_z_across_chan,2);
dead_channel=[sel_auto_artifact.labels(z_chan_novar>3)',num2cell(rel_novar_chan(z_chan_novar>3))]

% select threshold to mark segment as artifact
relchan_thres=0.1; 
rel_trial= sum(artis)./size(artis,1);
 sel_artis=rel_trial>relchan_thres;
arti_sp=sel_auto_artifact.samples(sel_artis,1:2);

% show info regarding channels in a figure
figure
subplot(1,2,1)
for e=1:numel(datainfo.elec_info.ana_labels.labels)
text(0,numel(datainfo.elec_info.ana_labels.labels)-e,strcat(datainfo.elec_info.ana_labels.labels{e},': ',datainfo.elec_info.ana_labels.freesurferDK{e,1}))
ylim([-2 numel(datainfo.elec_info.ana_labels.labels)])
end
axis off
subplot(1,2,2)
axis off
text(0,10,strcat('artifact elec:',[reject_channel{:}]))
hold on
text(0,5,strcat('no var elec:',[dead_channel{:}]))
ylim([0 12])

% check for channels to exclude & mark undetected artifact & remove wrong
% artifacts
cfg=[];
cfg.viewmode = 'vertical';
cfg.ylim =  'maxmin';
cfg.blocksize =30;
cfg.continuous  = 'yes';
cfg.artfctdef.xxx.artifact  =arti_sp;
artif_check=ft_databrowser(cfg,data_bip)

% add breakpoint for easier channel selection
sel_spikechan=input('spike channels? write as [{''}]')
sel_artichan=input('arti channels? write as [{''}]')
sel_arti=artif_check.artfctdef.xxx;

browsercheck.spikechan=sel_spikechan;
browsercheck.artichan=sel_artichan;
browsercheck.arti_sp=sel_arti;

datainfo.artifact_info=setfield(datainfo.artifact_info,strcat('browsercheck_',sel_ref),browsercheck);
datainfo.artifact_info.autoartifact=sel_auto_artifact;

save(info_file,'datainfo')
clear
close all
delete(gcf)