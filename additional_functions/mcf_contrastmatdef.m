% function to define contrastmatrizes to feed into mcf_rsacontrasts
% input contrastdef as defined by rsa_contrastmatrix
% contrast: name of selected contrast, either field in contrastdef or
% multiplication contrastx_mask_maskx
% also for setting up interaction: contrast_vs_contrast

function    contrast_mat=mcf_contrastmatdef(contrast_def,contrast)
vs_con=strfind(contrast,'_interaction_');

% seperate contrasts incase of vs
if isempty(vs_con)
    all_contrasts{:}=contrast;
else
    % seperate contrasts
    all_contrasts{1}=contrast(1:vs_con-1);
    all_contrasts{2}=contrast(vs_con+13:end);
end


for c=1:numel(all_contrasts)
    sel_def=all_contrasts{c};
    mask_def=strfind(sel_def,'_mask_');
    
    % multiply contrast with mask if defined
    if isempty(mask_def)
        all_contrast_mat{c}=getfield(contrast_def,sel_def);
    else
        % seperate contrasts
        num_mask=numel(mask_def);
        mask_def=sort(mask_def,'descend')
        if contains(sel_def,'trial_slidingavg_def')
            tmp_mask=ones(size(contrast_def.trial_slidingavg_def));
        else
            tmp_mask=ones(size(contrast_def.item_specific));
        end
        tmp_def=sel_def;
        for m=1:num_mask
            sel_mask=tmp_def(mask_def(m)+6:end);
            tmp_def=tmp_def(1:mask_def(m)-1);
            tmp_mask=tmp_mask.*getfield(contrast_def,sel_mask);
        end
        sel_contrast=getfield(contrast_def,tmp_def);
        
        if contains(sel_def,'trial_slidingavg_def')
            sel_contrast=repmat(reshape(sel_contrast,[1,size(sel_contrast)]),size(tmp_mask,1),1,1);
        else
        end
        
        all_contrast_mat{c}=sel_contrast.*tmp_mask;
    end
    
    
    % add contrasts
    
    if c==2
        all_contrast_mat{c}=4-all_contrast_mat{c};
    elseif c>2
        error('contrasts between more than two contrasts not defined')
    else
        all_contrast_mat{c}=2-all_contrast_mat{c};
        
    end
    if contains(sel_def,'trial_slidingavg_def')
        contrast_mat(c,:,:,:)=all_contrast_mat{c};
        nan_masks(c,:,:,:)=isnan(all_contrast_mat{c});
    else
        contrast_mat(c,:,:)=all_contrast_mat{c};
        nan_masks(c,:,:)=isnan(all_contrast_mat{c});
    end
    
end

if c>1
    nan_masks=squeeze(sum(nan_masks))==c;
    contrast_mat=squeeze(nansum(contrast_mat));
else
    nan_masks=squeeze(nan_masks);
    contrast_mat=squeeze(contrast_mat);
end
contrast_mat(nan_masks)=nan;
