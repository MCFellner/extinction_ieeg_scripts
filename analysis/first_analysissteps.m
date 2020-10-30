addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')
% first analysis steps

%% trialinfo coding as reminder
% 1. what trial number (position in presentation)?
% 2. which Phase?
% 3. which context was used?
% 4. what was the role of the video (A,B,C1,C2)
% 5. which item was shown?
% 6. which type of item was shown? % cs+/cs+=1;cs+/cs-=2;cs-/cs-=3;
% 7. what response was given?
% 8. cs (0/1) current cs+/cs-
% 9. us 0/1 (y/n)
%%%% SR logfile 10000
% 10. sample point trialonset
% 11. sample point videoonset
% 12. sample point cueonset
% 13. sample point us onset
% 14. sample point response 

%later added
%15. number of item rep in each block
%16 number of us in total

%% erps: check for significant erp electrodes: difference during item/
% difference anticipating us
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\erp\';
mkdir(path_out)

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
 
% segment data in different trial parts
% item window: -1 to 4 (imp
pre_item=1;
post_item=5.5;
stat_windows=[2 2.5;2.5 3;3 3.5;3.5 4;4 4.5]; % add rows for more windows

% downsample to smallest sr
sr=1000;

conditions={'B_switch','B_plus'};
cond1_def=[{'2'},{'==1'};{'6'},{'==1'}];% column, value (through eval also <= or ~=), definition across columns combined with &  
cond2_def=[{'2'},{'==1'};{'6'},{'==2'}];

for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)  
    load(strcat(path_preproc,sel_sub,'_data.mat'))

    % apply bipolar montage
    montage=datainfo.elec_info.bipolar.montage_withoutartichan;
    data = ft_apply_montage(data,montage);

    
    % cut trials
    trl(:,1)=datainfo.trigger.trigger_sp-(data.fsample.*pre_item);
    trl(:,2)=datainfo.trigger.trigger_sp+(data.fsample.*post_item);
    trl(:,3)=ones(numel(datainfo.trigger.trigger_sp),1).*(-1.*pre_item.*data.fsample);    
    cfg=[];
    cfg.trl=trl;
    data=ft_redefinetrial(cfg,data);
    
   % downsample to common sampling rate
    cfg=[];
    cfg.resamplefs      = sr;
    cfg.detrend='yes';
    data=ft_resampledata(cfg,data);
    % lowpass filter 
    cfg=[];
    cfg.lpfilter='yes';
    cfg.lpfreq=15;
    data=ft_preprocessing(cfg,data);
    
    % only select artifree trials in trialinfo
    trlinfo=datainfo.trialinfo;
    trlinfo=trlinfo(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree,:);
    
    cfg=[]
    cfg.trials=find(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree);
    data=ft_selectdata(cfg,data);
    
    
    % define two conditions
    for d=1:size(cond1_def,1)
    eval(strcat('tmp(:,d)=trlinfo(:,',cond1_def{d,1},')',cond1_def{d,2}));
    end
    sel_trials1=find(sum(tmp,2)==size(cond1_def,1));
    clear tmp
    
    for d=1:size(cond2_def,1)
    eval(strcat('tmp(:,d)=trlinfo(:,',cond2_def{d,1},')',cond2_def{d,2}));
    end
    sel_trials2=find(sum(tmp,2)==size(cond2_def,1));
    clear tmp
    
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=sel_trials1;
    erp1=ft_timelockanalysis(cfg,data);
    cfg.trials=sel_trials2;
    erp2=ft_timelockanalysis(cfg,data);
    
    cfg=[];
    cfg.baseline=[-0.2 0];
    erp1=ft_timelockbaseline(cfg,erp1);
    erp2=ft_timelockbaseline(cfg,erp2);

    % loop over channels
  for t=1:size(stat_windows,1) 
   sel_window=   stat_windows(t,:);
  parfor e=1:numel(erp1.label)
    sel_elec=erp1.label{e};
    % contrast
   cfg=[];
    cfg.latency     = sel_window;
   cfg.channel   = erp1.label{e};
   cfg.avgoverchan =  'yes';                   
   cfg.avgovertime =  'no';  
   
   % first level 
   cfg.method           = 'montecarlo';
   cfg.numrandomization = 10000;
   cfg.correctm         =  'cluster';
   cfg.correcttail      = 'prob';
   cfg.statistic ='indepsamplesT'
    design = zeros(1,size(erp1.trial,1) + size(erp2.trial,1));
    design(1,1:size(erp1.trial,1)) = 1;
    design(1,(size(erp1.trial,1)+1):(size(erp1.trial,1) + size(erp2.trial,1)))= 2;
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
    cfg.design = design;    
    stat=ft_timelockstatistics (cfg,erp1,erp2)   
 
    if isfield(stat,'posclusters')
        if ~isempty(stat.posclusters)
        pos_prob(e)=stat.posclusters(1).prob;
        else
        pos_prob(e)=1;
        end
    else
    pos_prob(e)=1;   
    end
    
   if isfield(stat, 'negclusters')
        if ~isempty(stat.negclusters)
        neg_prob(e)=stat.negclusters(1).prob;
        else
        neg_prob(e)=1;
        end
    else
    neg_prob(e)=1;   
   end   
  end
   erpstat.label=erp1.label;
   erpstat.cfg=cfg;
   erpstat.prob_pos=pos_prob;
   erpstat.prob_neg=neg_prob;
   erpstat.time_window= sel_window;
   erpstat.conditions=conditions;
   erpstat.cond1_def=cond1_def;
   erpstat.cond2_def=cond2_def;
   erpstat.sel_trials1=sel_trials1;
   erpstat.sel_trials2=sel_trials2;


  % save in specific folder
  folder_out=fullfile(path_out,strcat(conditions{1},'_vs_',conditions{2},'_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));
  mkdir(folder_out)
  save(fullfile(folder_out,strcat(sel_sub,'erpstat')),'erpstat')
  clear erpstat pos_prob neg_prob  trl
  end  
end

%% plot sig electrodes on surface brain

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_stat='D:\Extinction\iEEG\analysis\erp\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};


%conditions={'A_cs_2half','A_nocs_2half'};
%conditions={'A_cs_2half','A_nocs_2half'};
conditions={'B_switch','B_plus'};
sel_window=[3.5 4]; % add rows for more windows
alpha_def=0.05;

roi.vmpfc={'ctx-lh-lateralorbitofrontal','ctx-lh-medialorbitofrontal','ctx-rh-lateralorbitofrontal','ctx-rh-medialorbitofrontal'};
roi.acc= {'ctx-lh-caudalanteriorcingulate','ctx-rh-caudalanteriorcingulate', 'ctx-lh-rostralanteriorcingulate','ctx-rh-rostralanteriorcingulate'};
roi.ifg_r={'ctx-rh-parstriangularis','ctx-rh-parsopercularis','ctx-rh-parsorbitalis'};
roi.ifg_l={'ctx-lh-parstriangularis','ctx-lh-parsopercularis','ctx-lh-parsorbitalis'};
roi.dm_pfc_r ={'ctx-rh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal'};
roi.dm_pfc_l={'ctx-lh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal'};

roi.amy_r={'Right-Amygdala'};
roi.amy_l={'Left-Amygdala'};
roi.hip_l={'Left-Hippocampus'};
roi.hip_r={'Right-Hippocampus'};

roi.ventraltempocci_l={'ctx-lh-fusiform','ctx-lh-inferiortemporal','ctx-lh-lateraloccipital','ctx-lh-lingual','ctx-lh-middletemporal','ctx-lh-parahippocampal','ctx-lh-temporalpole'};
roi.ventraltempocci_r={'ctx-rh-fusiform','ctx-rh-inferiortemporal','ctx-rh-lateraloccipital','ctx-rh-lingual','ctx-rh-middletemporal','ctx-rh-parahippocampal','ctx-rh-temporalpole'};

for sub=1:numel(allsubs)
    sel_sub=allsubs{sub};
 info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)  
folder_in=fullfile(path_stat,strcat(conditions{1},'_vs_',conditions{2},'_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));
load(fullfile(folder_in,strcat(sel_sub,'erpstat')))
stat=erpstat;
sig_def{sub}=stat.prob_pos<alpha_def|stat.prob_neg<alpha_def;
all_elec{sub}=stat.label;
        
% get positions        
 [~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');     
all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK(ind,1);
end

% get labels and relative count in each area
[region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def, roi)

save(fullfile(folder_in,'results_table'),'region_count','subject_count','roi_count')
%
num_all=sum([region_count.absnumelecinregion{:}])
num_sig=sum([region_count.absnumsigelec_pos{:}])
rel_sig=sum([region_count.absnumsigelec_pos{:}])/sum([region_count.absnumelecinregion{:}])


% plot electrodes and count electrodes per region/pat

% elec to plot
elec_toplot.unit ='mm';
elec_toplot.coordsys ='mni';
load('D:\Extinction\iEEG\scripts\additional_functions\sel_colorseries.mat')

sig_ind=[];
all_label=[];
all_pos=[];
for n=1:numel(all_elec)
   sig_ind=[sig_ind;sig_def{n}'];
   all_label=[all_label;all_elec{n}];
   all_pos=[all_pos;all_elec_pos{n}];
end


views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
mesh.coordsys = 'mni';
hemispheres={'left','right'};
elec_def=[-1,1];
type_elec=[1,0];
sel_col=[1,0,0;0.5 0.5 0.5];
trans_sphere=[1 0.2];
for h=1:numel(hemispheres)
    sel_hemi=hemispheres{h};
    sel_elec_def=elec_def(h);
    load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_',sel_hemi,'.mat')));
    figure
    ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
    hold on
    % elec to plot
    elec_toplot.unit ='mm';
    elec_toplot.coordsys ='mni';
    % all elecs in selected hemisphere
   sel_elec_hemi=sign(all_pos(:,1))==sel_elec_def;

    for t=1:numel(type_elec)
    sel_elec=sel_elec_hemi&(sig_ind==type_elec(t));
    sel_trans=trans_sphere(t);
    elec_toplot.label=all_label(sel_elec,:);
    elec_toplot.elecpos=all_pos(sel_elec,:);
    elec_toplot.chanpos=all_pos(sel_elec,:);
    sel_color=sel_col(t,:);
    
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color,'facealpha',sel_trans);
    end
    view([-90 20]);
    material dull;
    view(squeeze(views(h,1,:))');
    c1=camlight(0,0);
    set(c1, 'style', 'infinite');

    view(squeeze(views(h,2,:))');
    c2=camlight(0, 0);
    set(c2, 'style', 'infinite');
    view(squeeze(views(h,3,:))');
    print('-f1','-r600','-dtiff',fullfile(folder_in,strcat(sel_hemi,'_lat.tiff'))) 
    view(squeeze(views(h,4,:))');
    print('-f1','-r600','-dtiff',fullfile(folder_in,strcat(sel_hemi,'_med.tiff'))) 
    clear c1 c2 
    close all
end



% %%
% 
% % time frequency for each electrode: logspaced
% 
% 
% % rsa: powerspctzrm
% %phase A:cs plus more similar, phase B: cs minus more similar?
% % learning curve rsa?
% %