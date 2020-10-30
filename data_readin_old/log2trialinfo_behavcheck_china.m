%% behavdata readin

%% preprocessing script EXTINCTION money loss data

% first: add toolboxes to your path

% add path with your fieldtrip toolbox
addpath ('D:\matlab_tools\fieldtrip-20160122'); 
% add path with additional functions
addpath ('D:\Extinction\iEEG\scripts\additional_functions');
ft_defaults % adds the right folders of fieldtrip

            
%% read in log file information (for trialinfo & for behavioral data analysis)


%  we have read in the eeg files, however we need more information: 
% 12. what trial number (position in presentation)?
% 13. which Phase?
% 14. which context was used?
% 15. what was the role of the video (A,B,C1,C2)
% 16. which item was shown?
% 17. which type of item was shown?
% 18. what response was given?
% 19. RT?
% 20. cs
% 21. us (y/n)

% this info we can get from the presentation logfiles

path_data='D:\Extinction\iEEG\rawdata\China\QichenGroup\';
path_out='D:\Extinction\iEEG\data\eeg\China\QichenGroup\preproc\trialinfo\';

mkdir(path_out);
all_dir=dir(path_data);
all_datasets=[{all_dir.name}];
tmp=cellfun(@numel,all_datasets)>2;
all_datasets=all_datasets(tmp);

for f=1:numel(all_datasets)
    sel_sub=all_datasets{f};
    path_sub=strcat(path_data,sel_sub,'\',sel_sub,'\');
    subdir=dir(strcat(path_sub,'*behavioral'));  
   
    path_log=strcat(path_sub,subdir.name,'\');
    logdir=dir(strcat(path_log,'logfile*'));
    
path_in=strcat(path_log,logdir.name,'\');

allfiles={'A&B','C'};
for f=1:numel(allfiles) 
sel_file=allfiles{f}
switch sel_file
    case 'A&B'
        allruns={'A','B'};
        allblocks=[1,2];
    case 'C'
         allruns={'C'};
         allblocks=3;
end

for r=1:numel(allruns)
sel_phase=allruns{r};
sel_block=allblocks(r);

    sel_log= dir(strcat(path_in,'*',sel_file,'_chinese.log'));

    file_name=strcat(path_in,sel_log.name);
    
    log = log2matSHORT(file_name); % info from logfile_x is now in this variable;  
    clear file_ind 

      

    
            start_trial = strncmp(strcat(sel_phase,'trial'),log{1,3},6);
            id_trials = find(start_trial);
            log_time_trials=log{1,4}(start_trial);
         phase_ind = ((r-1)*numel(id_trials))+1:r*numel(id_trials); % 1-72 or 73-144

            item_trials = find(strncmp('Item',log{1,3},4));
            item_trials = item_trials(phase_ind);            
            trlinfo(:,1) = log{1,1}(id_trials); % logfile trl number in first colum of trialinfo
            trlinfo(:,2)=ones(size(trlinfo)).*sel_block; % % number of task block
             
            % indices contexts, converts context into double, writes in trlinfo
            z = find(strncmp('context',log{1,3}, 7));
            z = z(phase_ind);

            vid1 = cell(1, length(z));
            vid2 = cell(1, length(z));
            
                for i=1:numel(z)
                    % context first index (context1_4)
                    vid1{i}=  log{1,3}{z(i)}(8);
                    % context second index
                    vid2{i}=  log{1,3}{z(i)}(10);
            
                end
                
            vid=cellfun(@str2num,strcat(vid1,vid2),'UniformOutput',false);
            trlinfo(:,3) = [vid{:}];
            clear z vid vid1 vid2 
            
            
            % indices items, converts items into double, writes in trlinfo
            % indices types, converts types into double, writes in trlinfo
            % z = find(strncmp('Item',log{1,3}, 4));
            % z = z(phase_ind);
            % z not needed, we've got item_trials
            item = zeros(length(item_trials), 1);
            type = zeros(length(item_trials), 1);
            
            for i=1:numel(item_trials)
                item(i)=str2double(log{1,3}{item_trials(i)}(5));
            %    type(i)=str2double(log{1,3}{item_trials(i)}(end));
             type(i)=0; % type coding is off, recode at the end
            end
            trlinfo(:,5) = item;
            trlinfo(:,6) = type;
                clear item type i 

        % --------------------
        % until here, same A-D
        % --------------------
        
        if f == 1 % only in Phase A and B and all response trials now, therefore...
            
            % responses (1,2,3 or 4 ); participant probably didn't respond! (or responded with SPACE or responded too early (before response screen) or responded no matter what....
            % if no answer to response trial: NaN
            % ... probably switched responses...
            % ... responses after us_event/no_us not included
            % always take first response
            y = zeros(length(trlinfo),1);
            RT = zeros(length(id_trials),1);
            for i = 1:length(id_trials)
                if (id_trials(i)+7)<numel(log{1,3})
                
               tmp_tr= log{1,3}(id_trials(i)+3:id_trials(i)+7);
               tmp_rt= log{1,4}(id_trials(i)+3:id_trials(i)+7);
                else (id_trials(i)+7)>numel(log{1,3})
               tmp_tr= log{1,3}(id_trials(i)+3:end);
               tmp_rt= log{1,4}(id_trials(i)+3:end);
                end
               
               tmp_allres=str2double(tmp_tr);
               res_ind= find(~isnan(tmp_allres));
               if isempty(res_ind)
                    y(i)= NaN; % didn't respond when response_trial -> missing
                    RT(i)=NaN;
               else
                    y(i)=tmp_allres(res_ind(1)); 
                    RT(i)=(tmp_rt(res_ind(1))-log{1,4}(item_trials(i)))./10;
               end
               
            end
            trlinfo(:,7) = y;
            
            % RT (response - item trial)
            % NaN: missing
            trlinfo(:,8) = RT;
            
            % cs1 vs. cs0
            cs = strncmp('cs',log{1,3}, 2);
            cs = find(cs);
             cs = cs(phase_ind);
            u = zeros(length(id_trials),1);
            for i=1:numel(cs)
                u(i)=str2double(log{1,3}{cs(i)}(3));
            end
            trlinfo(:,9) = u;         
            % keep y!
            % us_event (1) or no_us(0)
            w = zeros(length(trlinfo),1);
            a1 = strncmp(log{1,3},'us_event',7);
            a2 = strncmp(log{1,3},'no_us',7);
            b = find(a1+a2);
            b = b(phase_ind);

            for i = 1:length(b)
                if strcmp('no_us',log{1,3}(b(i)))
                    w(i) = 0;
                elseif strcmp('us_event',log{1,3}(b(i)))
                    w(i) = 1;
                end
            end
            trlinfo(:,10) = w;
           clear w b a1 a2 u cs
            
        elseif f == 2
            % responses and RTs
            response = strncmp(strcat(sel_phase,'response'),log{1,3}, 5); 
            id_response = find(response);
            resp = zeros(length(id_response), 1);
            %first_button_press = zeros(length(log{1,3}));
            RT = zeros(length(id_response), 1);
            y = zeros(length(trlinfo),1);
            for i = 1:length(id_response)
               tmp_tr= log{1,3}(id_response(i)-1:id_response(i)+1);
               tmp_rt= log{1,4}(id_response(i)-1:id_response(i)+1);
               tmp_allres=str2double(tmp_tr);
               res_ind= find(~isnan(tmp_allres));
               if isempty(res_ind)
                    y(i)= NaN; % didn't respond when response_trial -> missing
                    RT(i)=NaN;
               else
                    y(i)=tmp_allres(res_ind(1)); 
                    RT(i)=(tmp_rt(res_ind(1))-log{1,4}(item_trials(i)))./10;
               end
            end
            trlinfo(:,7) = y;
            
            % RT in trlinfo
            trlinfo(:,8) = RT;
            clear RT y 
            % all no_us -> not in trialinfo
            % pattern ABC -> A=1; B =2; C=3;
            h = find(strncmp('A',log{1,3}, 1));
            test = cell(length(id_trials),1);
            type = cell(length(id_trials),1);
            for i=1:numel(h)
              test{i}=log{1,3}{h(i)}(1:3);
              type{i}=log{1,3}{h(i)}(end);
            end          
            test = strrep(test, 'A', '1');
            test = strrep(test, 'B', '2');
            test = strrep(test, 'C', '3');
            test = (cellfun(@str2num,test))-120;
            type = cellfun(@str2double, type);
            trlinfo(:,4) = test;
            %trlinfo(:,8) = type;
           
            clear h test %type
        % (for A and B: 6 = response trial y/n deleted) 
        % [1:trialnumber, 2:phase index, 3:contexts, 4: role video (extinction pattern), 5:items, 6:types, 7:response, 8: RT, 9:cs, 10:us (y/n) 11: time vid start in logfile]
        % for C and D:
        end
    trlinfo(:,11)=log_time_trials;
    clear log_time_trials;
        
        
   log_allblock{sel_block}=trlinfo;
   clear trlinfo  
        
end
end


% concatenate all logs to define different types:
% cs+/cs+=1;cs+/cs-=2;cs-/cs-=3;

trlinfo=vertcat(log_allblock{:});

type_sum=[24,12,0];

for i=1:3
sum_us(i)=sum(trlinfo(trlinfo(:,5)==i,10));
sel_type=find(sum_us(i)==type_sum);
trlinfo(trlinfo(:,5)==i,6)=sel_type;
end


if strncmp(sel_sub,'1-DBX',5)
   % patient switched button codes in phase A & B
   trlinfo(trlinfo(:,2)<3,7)=5-trlinfo(trlinfo(:,2)<3,7);
   
end

save(strcat(path_out,sel_sub,'_trlinfo'),'trlinfo')
clear trlinfo

end




%% get some FIGURES for behavioral data per participant
path_info='D:\Extinction\iEEG\data\eeg\China\QichenGroup\preproc\data_check\behav\'; % define where you want to save your file:
path_in='D:\Extinction\iEEG\data\eeg\China\QichenGroup\preproc\trialinfo\';

mkdir(path_info)

allvps = {'1-DBX','2-LSY','3-JJL'};%, '02'}; 

for vp = 1%:numel(allvps)
    sel_vp=allvps{vp};
    load(strcat(path_in,sel_vp,'_trlinfo'))
    
    % organize trialinfo in matrix
    trialinfo=trlinfo;
    
    res = 4; % rating 1-4
    label = {'Dangerous','Safe'}; % 
    types = {'cs+/cs+','cs+/cs-','cs-/cs-'}; % only 3 types in EXT loss (cs+/cs+ = type1, cs+/cs- = type 2,cs-/cs- = type3) 
    items = 1:3; % only 3 items in EXT loss (2, 4, 6)
    block_def=[1 24;25 48;49 64]; % #of trials in blocks
    
    % which item is which type?
    type_item = zeros(1,length(items));
    responses_item = zeros(3, block_def(end,end)); 
    RT_item_alltrials = zeros(3, block_def(end,end));
    for i = 1:length(items)
        type_item(i) = mean(trialinfo(trialinfo(:,5) == items(i),6)); % mean gives exact value of itemtype
        responses_item(i,:) = trialinfo(trialinfo(:,5) == i,7); % give response from column 18 at (column 17 equals 1 or 2 or 3 (itemtype)
        tmp = trialinfo(trialinfo(:, 5) == i, 8); % gives RTs for every type
        RT_item_alltrials(i, :) = tmp';
    end
    
    % compute average over types
    for ty=1:numel(types)
    RT_type_alltrials(ty,:)= RT_item_alltrials(type_item==ty,:);
    responses(ty,:)=responses_item(type_item==ty,:);
    end
    
    
    % gives for every itemtype the responses throughout the experiment 
    
    % 1st line: type 1, 2nd line: type 2, 3rd line: type 3
    % gives reaction time for every itemtype (1st line: type 1, 2nd line: type 2, 3rd line: type 3)
    
    avg_response_type = nanmean(responses,2); % mean resonses per item type
    % 1st line: type 1, 2nd line: type 2, 3rd line: type 3
    % replace following line above if interested in median responses per item_type
    % avg_response_type = median(responses,2, 'omitnan');
    avg_RT_type_alltrials = nanmean(RT_type_alltrials,2); % mean RT per item type
    std_RT_type_alltrials = nanstd(RT_type_alltrials,0,2); % std RT per item type
    
    % median and mean responses and RT per BLOCK for every itemtype
    response_medianblock = zeros(length(items),size(block_def, 1));
    response_avgblock = zeros(length(items),size(block_def, 1));
    response_stdblock = zeros(length(items),size(block_def, 1));
    RT_avgblock = zeros(length(items),size(block_def, 1));
    RT_stdblock = zeros(length(items),size(block_def, 1));
    
    for b = 1:length(types)
        for i = 1:size(block_def, 1) % 4columns
            response_medianblock(b,i) = median(responses(b,(block_def(i,1):block_def(i,2))), 2, 'omitnan');
            response_avgblock(b,i) = nanmean(responses(b,(block_def(i,1):block_def(i,2))), 2);
            response_stdblock(b,i)= nanstd(responses(b,(block_def(i,1):block_def(i,2))));
            RT_avgblock(b,i) = nanmean(RT_type_alltrials(b,(block_def(i,1):block_def(i,2))), 2);
            RT_stdblock(b,i)= nanstd(RT_type_alltrials(b,(block_def(i,1):block_def(i,2))));
        end
    end
    
    clear tmp
    
   h= figure % initializes figure
    
        for i = 1:length(types)
            subplot(3,4,i)
            hold on 
            rectangle('Position',[0 0 16 res],'FaceColor',[1 .9 .9])
            rectangle('Position',[16 0 16 res],'FaceColor',[0.9 .9 .9])
            rectangle('Position',[32 0 16 res],'FaceColor',[1 .8 1])
            rectangle('Position',[48 0 16 res],'FaceColor',[1 1 .8])
        
        
            stem(responses(i,:))
            title(types{i})
            h = gca;
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
        end % creates first 3 figures with colored stem plot
    
        % mean responses per itemtype
        subplot(3,3,4)
        plot(avg_response_type)
        title('mean responses per type')
        hold on
        h = gca;
        h.XTick = 1:3;
        h.XTickLabel = types;
        h.YLim = [0 4];
        h.YTick = 0:res;
        h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
        h.YTickLabelRotation = 90;
        
        % median responses per itemtype per block
        for i = 1:length(types)
            subplot(3,3,5)
            hold on 
            plot(response_medianblock(i,:))
            title('median rating/block')
            h = gca;
            h.XTick = 1:4;
            h.YLim = [0 4];
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
            %errorbar(response_avgblock(i,:),response_stdblock(i,:))
        end
        legend(types)
    
        % mean rating for each block
        for ty = 1:length(types)
            subplot(3,3,6)
            hold on
            plot(response_avgblock(ty,:));
            %errorbar(1:4,response_avgblock(ty,:),response_stdblock(ty,:))
            title('mean rating/block')
            h = gca;
            h.XTick = 1:4;
            h.YLim = [0 4];
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
        end
        legend(types)
        
        % mean RT per itemtype
        subplot(3,3,7)
        hold on
        bar(1:length(types),avg_RT_type_alltrials)
        errorbar(1:length(types),avg_RT_type_alltrials,std_RT_type_alltrials,'.')
        title('mean RT per type')
        h = gca;
        h.XTick = 1:3;
        h.XTickLabel = types;
        
        % mean RT per itemtype
        for ty = 1:length(types)
            subplot(3,3,8)
            hold on
            plot(RT_avgblock(ty,:));
            %errorbar(1:4,RT_avgblock(ty,:),RT_stdblock(ty,:))
            title('mean RT/block')
            h = gca;
            h.XTick = 1:4;
        end
        legend(types)
      
        
        
        
    savefig (strcat(path_info,  sel_vp, '_rating_summary.fig'));   
    
    figure
    % plot rating for each item
    for i=1:numel(items)
    subplot(4,2,i)
     hold on 
            rectangle('Position',[0 0 16 res],'FaceColor',[1 .9 .9])
            rectangle('Position',[16 0 16 res],'FaceColor',[0.9 .9 .9])
            rectangle('Position',[32 0 16 res],'FaceColor',[1 .8 1])
            rectangle('Position',[48 0 16 res],'FaceColor',[1 1 .8])
        
        
            stem(responses_item(i,:))
            title(types{type_item(i)})
            h = gca;
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
    end
    
end







