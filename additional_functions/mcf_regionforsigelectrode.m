% subfunction, returns table with electrode counts for each region

% input (all organized as cell (e.g. each cell a subject)
% all_elec_pos: all electrode postions
% all_elec_label: anatomical labels for each electrode
% sig_def: vec with 0/1 definition sig elec
% roi: struct with definition of rois e.g.  roi.dm_pfc_l={'ctx_lh_rostralmiddlefrontal','ctx_lh_caudalmiddlefrontal'};

%%% fix: add roi option
%%% also: count for roi per sub

function [region_count,subject_count,roi_count]=mcf_regionforsigelectrode(all_elec_pos,all_elec_label,sig_def,roi_def)

% get labels and relative count in each area
all_label=[];
all_sig=[];
for n=1:numel(all_elec_label)
all_label=[all_label;all_elec_label{n}];
all_sig=[all_sig;sig_def{n}'];

% sig elec in each sub
pos_sig(n)=sum(sig_def{n}==1);
neg_sig(n)=sum(sig_def{n}==-1);
total_elec(n)=numel(all_elec_label{n});
end
subject_count=table(total_elec',pos_sig',neg_sig','VariableNames',{'numelec','possig','negsig'});

all_label=[all_label{:}];
all_regions=unique(all_label)';

for r=1:numel(all_regions)
sel_region=all_regions{r};
ind= strcmp(all_label,sel_region);
rel_sig_pos{r,1}=sum(all_sig(ind)==1)./sum(ind);  
rel_sig_neg{r,1}=sum(all_sig(ind)==-1)./sum(ind);  
absolute_num{r,1}= sum(ind);
sig_num_pos{r,1}= sum(all_sig(ind)==1);
sig_num_neg{r,1}= sum(all_sig(ind)==-1);

end

region_count=table(all_regions,rel_sig_pos,rel_sig_neg,absolute_num, sig_num_pos,sig_num_neg,'VariableNames',{'analabel','relsigelec_pos','relsigelec_neg','absnumelecinregion','absnumsigelec_pos','absnumsigelec_neg'})


% get count for each roi
all_rois=fieldnames(roi_def);

for r=1:numel(all_rois)
    sel_def=all_rois{r};
    sel_label=getfield(roi_def,sel_def)
    ind=[];
    for l=1:numel(sel_label)
    [ind_tmp]=find(strcmp(sel_label{l},all_label ));
    ind=[ind,ind_tmp];
    end
roi_rel_sig_pos{r,1}=sum(all_sig(ind)==1)./numel(ind);  
roi_rel_sig_neg{r,1}=sum(all_sig(ind)==-1)./numel(ind);  
roi_absolute_num{r,1}= numel(ind);
roi_sig_num_pos{r,1}= sum(all_sig(ind)==1);
roi_sig_num_neg{r,1}= sum(all_sig(ind)==-1);

end
roi_count=table(all_rois,roi_rel_sig_pos,roi_rel_sig_neg,roi_absolute_num, roi_sig_num_pos,roi_sig_num_neg,'VariableNames',{'roilabel','relsigelec_pos','relsigelec_neg','absnumelecinregion','absnumsigelec_pos','absnumsigelec_neg'})

