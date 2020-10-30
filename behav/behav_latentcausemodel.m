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

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
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

% replace NaNs with prior response
for c=1:3
    nan_ind=0;
    count=0;
   while ~isempty(nan_ind)
sel_trials=find(trlinfo(:,6)==c);
sel_response=trlinfo(sel_trials,7);
nan_ind=find(isnan(sel_response));
if any(nan_ind==1)
sel_response(nan_ind(1))=sel_response(nan_ind(1)+1);
sel_response(nan_ind(2:end))=sel_response(nan_ind(2:end)-1);
else
sel_response(nan_ind)=sel_response(nan_ind-1);
end

if count==100
    count=0;
if any(nan_ind==1)
sel_response(nan_ind(1))=sel_response(nan_ind(1)+1);
sel_response(nan_ind(2:end))=sel_response(nan_ind(2:end)+1);
else
sel_response(nan_ind)=sel_response(nan_ind+1);
end
end
    
trlinfo(sel_trials,7)=sel_response;
count=count+1;
   end
end


CR=(trlinfo(:,7)-1)./4;
all_data(sub).CR=CR;
all_data(sub).US=trlinfo(:,9);
CS(:,1)= trlinfo(:,6)==1;
CS(:,2)= trlinfo(:,6)==2;
CS(:,3)= trlinfo(:,6)==3;
CS(:,4)= trlinfo(:,2)==1;
CS(:,5)= trlinfo(:,2)==2;
CS(:,6)= trlinfo(:,2)==3;
all_data(sub).CS=CS;
 clear CS CR
 
 end
 opts = LCM_opts([])
 opts.M=1;
  opts.K=10;
  opts.stickiness = 20
 results = LCM_fit(all_data,opts)
  % check out latent cause model

 %data - [nSubjects x 1] structure containing the following fields:
    %           .CR - [nTrials x 1] conditioned response
    %           .CS - [nTrials x nCues] conditioned stimului
    %           .US - [nTrials x 1] unconditioned response