%% design matrizes for rsa contrasts


% general structure:
% 1 cond1, 0 cond2, Nan unused
% or for learningcurves rep numbers
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_out='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';
mkdir(path_out)
path_figs=fullfile(path_out,'figs')
mkdir(path_figs)
contrasts={'item_specific',...
    'cs_specific',...
    'type1to2_vs_type2to3',...
    'learn_cond1',...
    'learn_cond2',...
    'learn_cond3',...
    'in_block','all',...
    'block1','block2','block3',...
    'first_half_eachblock','second_half_eachblock',...
        'first_half_block1','second_half_block1',...
        'first_half_block2','second_half_block2',...
        'first_half_block3','second_half_block3',...
    'first2second_half_inblock',...
    'no_ustrials',...
    'trial_slidingavg_def'};

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
sym='sym'; % or no_sym defines symmetric or non symmetric contrasts
for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    contrast_def=struct;
    
    % only select artifree trials in trialinfo
    trlinfo=datainfo.trialinfo;
    trlinfo=trlinfo(datainfo.artifact_info.clean_trials.item.artifactfree&datainfo.artifact_info.clean_trials.us.artifactfree,:);
    
    % get sort ind for later viz of the matrizes
    trlinfo_org=trlinfo;
    [trlinfo,sort_ind]=sortrows(trlinfo,[2,6,15]);
    fig= figure
    
    % axis labels
    for b=1:3
        blockaxis(b)=round(mean(find((trlinfo(:,2)==b))));
        for cond=1:3
            typeaxis(cond,b)=round(mean(find((trlinfo(:,2)==b)&(trlinfo(:,6)==cond))));
        end
    end
    
    for c=1:numel(contrasts)
        sel_contrast=contrasts{c};
        %setup empty tmp matrix
        num_trial=size(trlinfo,1);
        contrast_mat=zeros(num_trial);
        switch sym
            case 'no_sym'
                % remove upper triangle (no same trial, no doubeling)
                upper=triu(ones(size(contrast_mat)));
                contrast_mat(logical(upper))=NaN;
            case 'sym'
                diag_def=diag(ones(size(contrast_mat,1),1));
                contrast_mat(logical(diag_def))=NaN;
                
        end
        switch sel_contrast
            case 'item_specific'
                
                sel_col=trlinfo(:,6);
                all_cat=unique(sel_col);
                for i=1:numel(all_cat)
                    tmp=sel_col==all_cat(i);
                    tmp_mat=tmp*tmp';
                    contrast_mat=contrast_mat+tmp_mat;
                end
                clear tmp tmp_mat
                
         
                
            case  'cs_specific'
                % exclude other block 3
                nan_col=trlinfo(:,2)==3;
                nan_mat=contrast_mat;
                nan_mat(nan_col,:)=nan;
                nan_mat(:,nan_col)=nan;
                
                sel_col=trlinfo(:,8)+1;
                all_cat=unique(sel_col);
                for i=1:numel(all_cat)
                    tmp=sel_col==all_cat(i);
                    tmp_mat=tmp*tmp';
                    contrast_mat=contrast_mat+tmp_mat;
                end
                contrast_mat=contrast_mat+nan_mat;
                clear tmp tmp_mat nan_mat nan_col
                
           
                
            case  {'type1to2_vs_type2to3'}
%                 nan_col=trlinfo(:,2)~=str2double(sel_contrast(end));
%                 nan_mat=contrast_mat;
%                 nan_mat(nan_col,:)=nan;
%                 nan_mat(:,nan_col)=nan;
                
                col_1=trlinfo(:,6)==1;
                col_2=trlinfo(:,6)==2;
                tmp1=col_2*col_1';
                col_3=trlinfo(:,6)==3;
                tmp2=(col_3*col_2').*2;
                
                tmp=tmp1+tmp2;
                tmp(tmp==0)=NaN;
                tmp(tmp==2)=0;
                
                contrast_mat=tmp+contrast_mat;%+nan_mat;
                clear tmp1 tmp2 tmp3  col_1 col_2 col_3  nan_col
                
            case  {'learn_cond1',...
                    'learn_cond2',...
                    'learn_cond3',}
                contrast_mat=nan(num_trial);
                % matrix for learning curves: correlate in each cond consecutive
                % trials (i.e. one of the diagonal
                  trlinfo(:,17)=NaN;
                trlinfo(trlinfo(:,2)==1,17)=trlinfo(trlinfo(:,2)==1,15);
                trlinfo(trlinfo(:,2)==2,17)=trlinfo(trlinfo(:,2)==2,15)+24;
                trlinfo(trlinfo(:,2)==3,17)=trlinfo(trlinfo(:,2)==3,15)+48;
                sel_trial=trlinfo(:,6)==str2double(sel_contrast(11));
                num_sel_trials=trlinfo(sel_trial,17);
                
                
                tmp=diag(num_sel_trials(2:end),-1);
                contrast_mat(sel_trial,sel_trial)=tmp;
                contrast_mat(contrast_mat==0)=NaN;

                
                clear tmp sel_trial num_sel_trials
             case 'in_block'
                def_vec=trlinfo(:,2)==1;
                contrast_mat1=def_vec*def_vec';
                def_vec=trlinfo(:,2)==2;
                contrast_mat2=def_vec*def_vec';
                   def_vec=trlinfo(:,2)==3;
                contrast_mat3=def_vec*def_vec';
                contrast_mat=contrast_mat1+contrast_mat2+contrast_mat3;
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec contrast_mat1 contrast_mat2 contrast_mat3
            case 'all'
                contrast_mat=contrast_mat+(ones(num_trial,num_trial));
                
            case {'block1',  'block2', 'block3'}
                def_vec=trlinfo(:,2)==str2double(sel_contrast(end));
                contrast_mat=def_vec*def_vec';
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec
                
            case {'first_half_eachblock'}
                def_vec=(trlinfo(:,2)<=2 & trlinfo(:,15)<=12)|(trlinfo(:,2)==3 & trlinfo(:,15)<=8);
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec
            case {'first_half_block1','first_half_block2'}
                 def_vec=(trlinfo(:,2)==str2double(sel_contrast(end)) & trlinfo(:,15)<=12);
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec 
 
            case {'second_half_eachblock'}
                def_vec=(trlinfo(:,2)<=2 & trlinfo(:,15)>12)|(trlinfo(:,2)==3 & trlinfo(:,15)>8);
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec
                
                
            case {'second_half_block1','second_half_block2'}
                 def_vec=(trlinfo(:,2)==str2double(sel_contrast(end)) & trlinfo(:,15)>12);
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec  
            case {'first_half_block3'}
                 def_vec=(trlinfo(:,2)==3) & trlinfo(:,15)<=8;
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec     
            case {'second_half_block3'}
                 def_vec=(trlinfo(:,2)==3) & trlinfo(:,15)>8;
                contrast_mat=contrast_mat+(def_vec*def_vec');
                contrast_mat(contrast_mat==0)=NaN;
                clear def vec      
                
            case {'first2second_half_inblock'}
                def_vec1=(trlinfo(:,2)<=2 & trlinfo(:,15)<=12)|(trlinfo(:,2)==3 & trlinfo(:,15)<=8);
                def_vec2=(trlinfo(:,2)<=2 & trlinfo(:,15)>12)|(trlinfo(:,2)==3 & trlinfo(:,15)>8);
                
                block_mat=((trlinfo(:,2)==1)*(trlinfo(:,2)==1)')+...
                    ((trlinfo(:,2)==2)*(trlinfo(:,2)==2)')+...
                    ((trlinfo(:,2)==3)*(trlinfo(:,2)==3)');
                contrast_mat=contrast_mat+((def_vec1*def_vec2').*block_mat);
                
                contrast_mat(contrast_mat==0)=NaN;
                clear def_vec1 def_vec2 block_mat
            case 'no_ustrials'    
               def_vec=trlinfo(:,9)==0;
                contrast_mat=contrast_mat+(def_vec*def_vec');
                
                contrast_mat(contrast_mat==0)=NaN;
                clear def_vec
            case 'trial_slidingavg_def'
                % here a contrast_mat for every rep is created
                trlinfo(:,17)=NaN;
                trlinfo(trlinfo(:,2)==1,17)=trlinfo(trlinfo(:,2)==1,15);
                trlinfo(trlinfo(:,2)==2,17)=trlinfo(trlinfo(:,2)==2,15)+24;
                trlinfo(trlinfo(:,2)==3,17)=trlinfo(trlinfo(:,2)==3,15)+48;
                
                contrast_mat3d=nan(max(trlinfo(:,17))-4,num_trial,num_trial);
                for i=1:(max(trlinfo(:,17))-4)
                    min_rep=i;
                    max_rep=i+4;
                    def_vec=trlinfo(:,17)>=min_rep & trlinfo(:,17)<=max_rep;
                    contrast_mat3d(i,:,:)=contrast_mat+(def_vec*def_vec');
                end
                contrast_mat3d(contrast_mat3d==0)=NaN;

                contrast_mat=squeeze(nansum(contrast_mat3d,1));
                
        end
        
        % add field to contrast_def
        switch sel_contrast
            case 'trial_slidingavg_def'
                contrast_def=setfield(contrast_def,sel_contrast,contrast_mat3d);
                
            otherwise
                contrast_def=setfield(contrast_def,sel_contrast,contrast_mat);
                
        end
        % plot design mats
        subplot(4,9,c)
        imagesc(contrast_mat,[-1,max(max(contrast_mat))])
        title(sel_contrast)
        xticks(blockaxis)
        xticklabels({'acq','ext','test'})
        yticks(reshape(typeaxis,1,[]))
        yticklabels(repmat({'+','+/-','-'},1,3))
        
        clear contrast_mat
    end
    
    subplot(4,9,c+1)
    imagesc(((trlinfo(:,2).*10)+trlinfo(:,6))*((trlinfo(:,2).*10)+trlinfo(:,6))');
    title('trial combis')
    colormap('jet')
    xticks(blockaxis)
    xticklabels({'acq','ext','test'})
    yticks(reshape(typeaxis,1,[]))
    yticklabels(repmat({'+','+/-','-'},1,3))
    
    
    savefig(fig,fullfile(path_figs,strcat(sel_sub,'_mat')),'compact')
    
    contrast_def.used_trlinfo=trlinfo;
    contrast_def.org_trlinfo=trlinfo_org;
    contrast_def.sortind_org2usedtrlinfo=sort_ind;
    
    % save
    save(fullfile(path_out,strcat(sel_sub,'_contrast_mat_',sym)),'contrast_def')
    clear  contrast_def trlinfo trlinfo_org sort_ind blockaxis typeaxis num_trial sel_col upper
    close all
end