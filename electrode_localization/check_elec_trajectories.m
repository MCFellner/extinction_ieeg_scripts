% check electrode position 
%do tracetories and distances make sense?

path_data='D:\Extinction\iEEG\';
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
        'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};
path_fig='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\check_electrode_pos\';
      mkdir(path_fig)
      
 for i=1:numel(allsubs)
% save labels in file & as excel sheet
sel_sub=allsubs{i};
info_file=strcat(path_info,sel_sub,'_datainfo');
load(info_file)

 elecs=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));
  fig=  figure
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        sel_channel=datainfo.elec_info.reordered.sorted_ieeg(strcmp(datainfo.elec_info.reordered.cell_ieeg(:,1),sel_elec));
  
        % check for each electrode the distance between elec matrix
   [~,~,ind]=intersect(sel_channel,datainfo.elec_info.elec_ct_mr.label,'stable');
   sel_pos=datainfo.elec_info.elec_ct_mr.elecpos(ind,:);
   tmp = pdist(sel_pos);
   dist_mat = squareform(tmp);
   ind_e2e=logical(diag(ones(1,numel(sel_channel)-1),-1));
   dist_e2e=dist_mat(ind_e2e);
   mean_dist=mean(dist_e2e);
   std_dist=std(dist_e2e);
   dist_diffmean_mat=abs(diff(dist_mat))-mean_dist;
   
   subplot(ceil(numel(elecs)./3),6,(e*2)-1)
   imagesc(dist_mat);
   %xticklabels(sel_channel);
   yticks(1:numel(sel_channel))
   yticklabels(sel_channel);
   title(strcat(sel_elec,': dist all elecs'))
  colorbar
   subplot(ceil(numel(elecs)./3),6,(e*2))
   imagesc(dist_diffmean_mat,[0 3]);
  % xticklabels(sel_channel);
  yticks(1:numel(sel_channel))
   yticklabels(sel_channel);
   title(strcat(sel_elec,'mean=',num2str(mean_dist),'std=',num2str(std_dist)))
   colorbar
    end
    savefig(fig,strcat(path_fig,sel_sub,'electrode_distance_check'))
    close all
keep path_fig path_data path_info allsubs

 end