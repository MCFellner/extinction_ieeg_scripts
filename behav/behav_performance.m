addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')
addpath('D:\matlab_tools\LCM-master')
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


%% analyze behav data


path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';
path_out='D:\Extinction\iEEG\analysis\erp\';
mkdir(path_out)

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
% 
% missing responses in p_sub07
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub08'};
 for sub=1:length(allsubs)
    sel_sub=allsubs{sub};
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)  
 
    trlinfo=datainfo.trialinfo;
    %   get response curves for each condition
%     cond1=trlinfo(trlinfo(:,6)==1,:);
%     cond2=trlinfo(trlinfo(:,6)==2,:);
%     cond3=trlinfo(trlinfo(:,6)==3,:);

% % replace NaNs with prior response
% for c=1:3
%     nan_ind=0;
%     count=0;
%    while ~isempty(nan_ind)
% sel_trials=find(trlinfo(:,6)==c);
% sel_response=trlinfo(sel_trials,7);
% nan_ind=find(isnan(sel_response));
% if any(nan_ind==1)
% sel_response(nan_ind(1))=sel_response(nan_ind(1)+1);
% sel_response(nan_ind(2:end))=sel_response(nan_ind(2:end)-1);
% else
% sel_response(nan_ind)=sel_response(nan_ind-1);
% end
% 
% if count==100
%     count=0;
% if any(nan_ind==1)
% sel_response(nan_ind(1))=sel_response(nan_ind(1)+1);
% sel_response(nan_ind(2:end))=sel_response(nan_ind(2:end)+1);
% else
% sel_response(nan_ind)=sel_response(nan_ind+1);
% end
% end
%     
% trlinfo(sel_trials,7)=sel_response;
% count=count+1;
%    end
% end

   %   get response curves for each condition
  conditions={'cs+cs+','cs+cs-','cs-cs-'};
   figure
   for c=1:3 
   cond{c}=trlinfo(trlinfo(:,6)==c,:);
   subplot(2,3,c)
      hold on 
      y=cond{c}(:,7);
      y=smoothdata(y,'movmean',3);
     plot(y,'--x')
     y_all{sub,c}=cond{c}(:,7);
     x_all{sub,c}=1:numel(y_all{sub,c});
    %plot dots for every us
    if c<=2
    ind_us=find(cond{c}(:,9)==1);
    scatter(ind_us,ones(size(ind_us)),'r','o')
    end 
    title([sel_sub, conditions{c}])
   end
     
   for c=1:3 
   cond{c}=trlinfo(trlinfo(:,6)==c,:);
   subplot(2,3,4)
      hold on 
      x=find(cond{c}(:,9)==0);
      x_nous{sub,c}=x;
      y=cond{c}(x,7);
      y_nous{sub,c}=y;
      y=smoothdata(y,'movmean',3);
     plot(x,y,'--x')
    legend(conditions)
    title([sel_sub, conditions{c}])
   end  
  
 end

 % add nans for missing trials in last sub
 x_all{end,1}=[x_all{end,1},NaN];
  y_all{end,1}=[y_all{end,1};NaN];
 x_all{end,2}=[x_all{end,2},NaN];
  y_all{end,2}=[y_all{end,2};NaN];
  

 x_nous{end,2}=[x_nous{end,2};NaN];
  y_nous{end,2}=[y_nous{end,2};NaN]; 
  
 % plot average response with and without
 
 for c=1:size(y_all,2)
 for sub=1:size(y_all,1)
 % for all use linspaced x
 y_all_mat(sub,:)=y_all{sub,c};
 sm_y_all_mat(sub,:) = smoothdata(y_all{sub,c},'movmean',3);

 % for no use use average x
 y_nous_mat(sub,:)=y_nous{sub,c};
  sm_y_nous_mat(sub,:)=smoothdata(y_nous{sub,c},'movmean',3);

  x_nous_mat(sub,:)=x_nous{sub,c};

 end    
 y_all_mean(c,:)=nanmean(y_all_mat);
  y_all_std(c,:)=nanstd(y_all_mat);
x_nous_mean{c}=nanmean(x_nous_mat);
 y_nous_mean{c}=nanmean(y_nous_mat);
  y_nous_std{c}=nanstd(y_nous_mat);
  
   sm_y_all_mean(c,:)=nanmean(sm_y_all_mat);
  sm_y_all_std(c,:)=nanstd(sm_y_all_mat);
x_nous_mean{c}=nanmean(x_nous_mat);
 sm_y_nous_mean{c}=nanmean(sm_y_nous_mat);
  sm_y_nous_std{c}=nanstd(sm_y_nous_mat);
 clear y_all_mat y_nous_mat x_nous_mat  sm_y_all_mat sm_y_nous_mat
 end
 
 
 % plot boundedline plots
 figure
hold on
fig_stuff=subplot(2,2,1)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=1:size(y_all_mean,2);
    y1=y_all_mean(c,:);
    b1=y_all_std(c,:)./sqrt(size(y_all,2));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
 title('all responses')
 end
 
 fig_stuff=subplot(2,2,2)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=x_nous_mean{c};
    y1=y_nous_mean{c};
    b1=y_nous_std{c}./sqrt(size(y_all,2));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
  title('no post us responses')

 end
 
 fig_stuff=subplot(2,2,3)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=1:size(y_all_mean,2);
    y1=sm_y_all_mean(c,:);
    b1=sm_y_all_std(c,:)./sqrt(size(y_all,2));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
 title('all responses smoothed')
 end
 
 fig_stuff=subplot(2,2,4)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=x_nous_mean{c};
    y1=sm_y_nous_mean{c};
    b1=sm_y_nous_std{c}./sqrt(size(y_all,2));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
  title('no post us responses smoothed')

 end
 
 % 