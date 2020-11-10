% stats plotter

%plot results of rsa stats in organized plots

% file definition
path_stats='D:\Extinction\iEEG\analysis\rsa\powlogscale_timeslide_z_crosstrials_toi2000to4000\stats';
path_out='D:\Extinction\iEEG\analysis\rsa\powlogscale_timeslide_z_crosstrials_toi2000to4000\';

% adapt file definition to load the wanted files

% define whether to plot multiplots or single plots
cfg_plot.multiplot='yes';
load('D:\matlab_tools\jet_grey.mat')
cfg_plot.def_colormap=jet_grey;

all_rois={'hip_l','hip_r','vmpfc','ifg','dm_pfc','amy_r','amy_l','ventraltempocci'};
all_contrasts=dir(path_stats);
all_contrasts={all_contrasts(:).name}';
all_contrasts=all_contrasts(cellfun(@numel,all_contrasts)>2);
for r=1:numel(all_rois)
cfg_plot.contrast='*'; %* if irrelevant
cfg_plot.roi=all_rois{r}; %* if irrelevant
cfg_plot.path_stats=path_stats;
fig=mcf_rsastatsmultiplot(cfg_plot)

    path_fig=fullfile(path_out,'summary_fig',all_rois{r})
    mkdir(path_fig)

savefig(fig,fullfile(path_fig,['summary.fig']),'compact')
close all
end


for c=1:numel(all_contrasts)
cfg_plot.contrast=all_contrasts{c}; %* if irrelevant
cfg_plot.roi='*'; %* if irrelevant
cfg_plot.path_stats=path_stats;

fig=mcf_rsastatsmultiplot(cfg_plot)

    path_fig=fullfile(path_out,'summary_fig',all_contrasts{c})
    mkdir(path_fig)
savefig(fig,fullfile(path_fig,['summary.fig']),'compact')
close all
end