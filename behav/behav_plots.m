%% analyse behav: individual responses

 path_trlinfo='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\'

selvps = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08',...
           'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
          'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20'};


for vp=1:numel(selvps)
% plot responses for each type across experiment


load(strcat(path_trlinfo,selvps{vp},'_datainfo.mat'))
trialinfo=datainfo.trialinfo;
label={'Dangerous','Safe'};


types={'cs+/cs+','cs+/cs-','cs-/cs-'};
items=1:3;


for i=1:3
type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
tmp=trialinfo(trialinfo(:,5)==i,7);
responses(i,:)=[tmp',nan(1,length(responses)-numel(tmp))];
tmp_us=trialinfo(trialinfo(:,5)==i,9);
us_vec(i,:)=[tmp_us',nan(1,length(tmp_us)-numel(tmp_us))]
end

figure
for i=1:3
subplot(2,3,i)
stem(responses(i,:))
title(types{type_item(i)})
end
figure
for ty=1:3
avg_response_type(ty,:)=nanmean(responses(type_item==ty,:),1) 
avg_us(ty,:)=us_vec(type_item==ty,:);
subplot(2,2,ty)
hold on
rectangle('Position',[0 0 24 4],'FaceColor',[1 .9 .9])
rectangle('Position',[24 0 24 4],'FaceColor',[0.9 .9 .9])
rectangle('Position',[48 0 16 4],'FaceColor',[1 .8 1])
rectangle('Position',[64 0 16 4],'FaceColor',[1 1 .8])


stem(avg_response_type(ty,:))
title(types{ty})
h=gca;
h.YTick=[0.1;4-(0.2*4)];
h.YTickLabel=label;
h.YTickLabelRotation=90;

end
end
clear avg_response_type tmp tmp_us

% close all

%% analyse behavior: average responses

 path_trlinfo='D:\Extinction\iEEG\data\preproc\trialinfo\'
% path_trlinfo='/Volumes/Untitled/Extinction/EEG_noES/data_mishal/trialinfo/';


for vp=1:numel(selvps)
% plot responses for each type across experiment
load(strcat(path_trlinfo,selvps{vp},'_trlinfo.mat'))
trialinfo=trlinfo;
types={'cs+/cs+','cs+/cs-','cs-/cs-'};
items=1:3;




for i=1:3
type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
tmp=trialinfo(trialinfo(:,5)==i,7);
responses(i,:)=[tmp',nan(1,64-numel(tmp))];
end

% average responses for each  n x type x response

for ty=1:3
  avg_response_type(vp,ty,:)=nanmean(responses(type_item==ty,:),1) 
    
end
end
%%
figure
hold on
fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;


for ty=1:3
    x1=1:size(avg_response_type,3)
    y1=squeeze(nanmean(avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(avg_response_type(:,ty,:),1))./sqrt(numel(selvps));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');

%plot(squeeze(nanmean(avg_response_type(:,ty,:))))
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
h.YTickLabelRotation=90;
end
legend(reshape(repmat(types,2,1),[],1))
plot([24 24],[1 4],':k')
plot([48 48],[1 4],':k')

%%
window=2;
% smooth trajectories
move_avg=ones(1,window)/window;
for vp=1:numel(selvps)
    for ty=1:3
     tmp_avg_response_type(vp,ty,:)=avg_response_type(vp,ty,:);
           
     % interpolate NaNs
     nan_ind=  find(isnan(tmp_avg_response_type(vp,ty,:)))
     
     while ~isempty(nan_ind)  
         for ind=1:numel(nan_ind)
             j=nan_ind(ind);
             if j>2 & j<(size(tmp_avg_response_type,3)-1)
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j-1:j+1)));
             elseif j<=2
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j+1)));
             elseif j>=(size(avg_response_type,3)-1)
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j-1)));   
             end
         end
        nan_ind=  find(isnan(tmp_avg_response_type(vp,ty,:)))
     end
     
filt_avg_response_type(vp,ty,:) = filtfilt(move_avg,1,squeeze(tmp_avg_response_type(vp,ty,:)))
% plot average trajectories
end
end


figure
hold on
fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;
cmap_default=cmap_default([1 2 4],:)


for ty=1:3
    x1=1:size(avg_response_type,3)
    y1=squeeze(nanmean(filt_avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(filt_avg_response_type(:,ty,:),1))./sqrt(numel(selvps));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');

%plot(squeeze(nanmean(avg_response_type(:,ty,:))))
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
h.YTickLabelRotation=90;
end
legend(reshape(repmat(types,2,1),[],1))
plot([24 24],[1 4],':k')
plot([48 48],[1 4],':k')



for t=1:size(filt_avg_response_type,3)
[htest(t),p(t)]=ttest(squeeze(filt_avg_response_type(:,1,t)),squeeze(filt_avg_response_type(:,2,t)))
end

figure
for ty=1:3
subplot(2,2,ty)
hold on
rectangle('Position',[0 0 24 4],'FaceColor',[1 .9 .9])
rectangle('Position',[24 0 24 4],'FaceColor',[0.9 .9 .9])
rectangle('Position',[48 0 16 4],'FaceColor',[1 .8 1])


stem(squeeze(nanmean(filt_avg_response_type(:,ty,:))))
title(types{ty})
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
h.YTickLabelRotation=90;


end

figure

contrasts=[1,2;1,3;2,3];
contrast_label={'cs+/cs+ vs cs+/cs-','cs+/cs+ vs cs-/cs-','cs+/cs- vs cs-/cs-'};
for con=1:3
    con1=contrasts(con,1);
    con2=contrasts(con,2);

subplot(2,3,con)
hold on
rectangle('Position',[0 -2 24 4],'FaceColor',[1 .9 .9])
rectangle('Position',[24 -2 24 4],'FaceColor',[0.9 .9 .9])
rectangle('Position',[48 -2 16 4],'FaceColor',[1 .8 1])
[htest,p,~,tstat]=ttest(squeeze(filt_avg_response_type(:,con1,:)),squeeze(filt_avg_response_type(:,con2,:)))


h=gca;
stem(squeeze(nanmean(filt_avg_response_type(:,con1,:)-filt_avg_response_type(:,con2,:))))
hold on
plot(htest)
title(types{ty})
h=gca;
h.YLim=[-2,2];
h.YTick=[-2;2];
h.YTickLabel={'less safe','safer'};
title(contrast_label{con})


end
%% average rating each block

block_def=[1 24;25 48;49 64];
figure
hold  on
for b=1:3
    response_avgblock(:,b)=nanmean(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),1),3);
   response_stdblock(:,b)=(nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1))./sqrt(numel(selvps));
    response_avgblocksub(:,:,b)=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);

end
for ty=1:3
  errorbar(1:4,response_avgblock(ty,:),response_stdblock(ty,:))

end
legend('cs+/cs+','cs+/cs-','cs-/cs-')
h=gca;
h.XLim=[0.5,4.5];
h.XTick=[1:4];
h.XTickLabel={'learn1','learn2','test1','test2',}
%% split performance block 3

block_def=[1 24;25 48;49 64];
figure
hold  on
 b=3
 response_avgsumblock3=nanmean(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1);
 response_stdsumblock3=nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1)./sqrt(numel(selvps));

all_response_avgblock3=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);

figure
subplot(1,3,1)

bar(response_avgsumblock3)
hold on
errorbar(1:3,response_avgsumblock3,response_stdsumblock3)
scatter(reshape(repmat((1:3),numel(selvps),1),numel(selvps)*3,1),reshape(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),numel(selvps)*3,1))
%plot(nansum(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3)')
h=gca;

h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
title('ABC')

avg_response_typeC=(avg_response_type(:,:,block_def(b,1):block_def(b,2)));
%% split performance block 4 in A & B



 path_trlinfo='D:\Extinction\iEEG\data\preproc\trialinfo\'
 %path_trlinfo='/Volumes/Untitled/Extinction/EEG_noES/data_mishal/trialinfo/';

%selvps={'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'};


for vp=1:numel(selvps)
% plot responses for each type across experiment

load(strcat(path_trlinfo,selvps{vp},'_trlinfo.mat'))
trialinfo=trlinfo;




trialinfoA=trialinfo(trialinfo(:,2)==4&trialinfo(:,4)==1,:);
trialinfoB=trialinfo(trialinfo(:,2)==4&trialinfo(:,4)==2,:);
 

for i=1:3
type_itemA(i)=nanmean(trialinfoA(trialinfoA(:,5)==i,6));
responsesA(i,:)=trialinfoA(trialinfoA(:,5)==i,7);

type_itemB(i)=nanmean(trialinfoB(trialinfoB(:,5)==i,6));
responsesB(i,:)=trialinfoB(trialinfoB(:,5)==i,7);
end

% average responses for each  n x type x response

for ty=1:3
  avg_response_typeB(vp,ty,:)=nanmean(responsesB(type_itemB==ty,:),1) 
  avg_response_typeA(vp,ty,:)=nanmean(responsesA(type_itemA==ty,:),1) 

end

end

all_response_avgblock4A=nanmean(avg_response_typeA,3);
all_response_avgblock4B=nanmean(avg_response_typeB,3);

 response_avgsumblock4A=nanmean(nanmean(avg_response_typeA(:,:,:),3),1);
 response_stdsumblock4A=nanstd(nanmean(avg_response_typeA(:,:,:),3),1)./sqrt(numel(selvps));

 response_avgsumblock4B=nanmean(nanmean(avg_response_typeB(:,:,:),3),1);
 response_stdsumblock4B=nanstd(nanmean(avg_response_typeB(:,:,:),3),1)./sqrt(numel(selvps));


subplot(1,3,2)
bar(response_avgsumblock4A)
hold on
errorbar(1:3,response_avgsumblock4A,response_stdsumblock4A)
% individual color and size of scatter dots
%sz = repmat(1:numel(selvps),1,4).*10;
%c = repmat(linspace(1,10,numel(selvps)),1,4);

scatter(reshape(repmat((1:3),numel(selvps),1),numel(selvps)*3,1),reshape(nanmean(avg_response_typeA(:,:,:),3),numel(selvps)*3,1))
%plot(nansum(avg_response_typeA(:,:,:),3)')

h=gca;
h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
title('ABA')

subplot(1,3,3)
bar(response_avgsumblock4B)
hold on
scatter(reshape(repmat((1:3),numel(selvps),1),numel(selvps)*3,1),reshape(nanmean(avg_response_typeB(:,:,:),3),numel(selvps)*3,1))
errorbar(1:3,response_avgsumblock4B,response_stdsumblock4B)
%plot(nansum(avg_response_typeB(:,:,:),3)')

h=gca;
h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
title('ABB')


figure
hold on
errorbar(1:3,response_avgsumblock4A,response_stdsumblock4B)
errorbar(1:3,response_avgsumblock4B,response_stdsumblock4B)
errorbar(1:3,response_avgsumblock3,response_stdsumblock4B)
h=gca;
h.XTick=[1:3];
h.XLim=[0.5,3.5];

h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
legend('ABA','ABB','ABC')
ylabel({'average response';'safe                                                  dangerous'})

%% plot time course test
%addpath(genpath('E:\matlab_tools\boundedline\kakearney-boundedline-pkg-8179f9a\boundedline'))


figure
fig_stuff=subplot(1,3,1)
cmap_default=fig_stuff.ColorOrder;
for ty=1:3
    x1=1:size(avg_response_typeA,3)
    y1=squeeze(nanmean(avg_response_typeA(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeA(:,ty,:),1))./sqrt(numel(selvps));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end
title('Test in A')

h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};

subplot(1,3,2)
for ty=1:3
    x1=1:size(avg_response_typeB,3)
    y1=squeeze(nanmean(avg_response_typeB(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeB(:,ty,:),1))./sqrt(numel(selvps));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end


h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
title('Test in B')

h.YTickLabel={'Dangerous','Safe'};
subplot(1,3,3)

for ty=1:3
    x1=1:size(avg_response_typeC,3)
    y1=squeeze(nanmean(avg_response_typeC(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeC(:,ty,:),1))./sqrt(numel(selvps));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
title('Test in C')
h.YTickLabel={'Dangerous','Safe'};
legend('cs+/cs+','','cs+/cs-','','cs-/cs-','');

%%

figure
fig_stuff=subplot(1,4,1)
cmap_default=fig_stuff.ColorOrder;

contexts={'A','B','C'};
conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
    figure
for ty=1:3
fig_stuff=subplot(1,3,ty)

cmap_default=fig_stuff.ColorOrder;
    
for con=1:3
    
    eval(strcat('tmp=avg_response_type',contexts{con}),';')
    x1=1:size(tmp,3)
    y1=squeeze(nanmean(tmp(:,ty,:),1));
    b1=squeeze(nanstd(tmp(:,ty,:),1))./sqrt(numel(selvps)); 
boundedline(x1, y1, b1, 'cmap',cmap_default(con,:),'transparency',0.1,'alpha');

end
title(conds{ty})

h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
end

legend('Test in A','','Test in B','','Test in C','');

%% plot distributions


% difference A vs b in each subject

    figure

diffab=all_response_avgblock4A-all_response_avgblock4B;

conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
    subplot(1,3,1)
    hist(diffab)
    legend (conds)
    title('difference A vs B')
    
    
diffac=all_response_avgblock4A-all_response_avgblock3;
    subplot(1,3,2)
    
hist(diffac)
    legend (conds)
    title('difference A vs C')
    
    diffbc=all_response_avgblock4B-all_response_avgblock3;
    subplot(1,3,3)
    
hist(diffbc)
    legend (conds)
    title('difference B vs C')

%%

    figure

diffab=all_response_avgblock4A-all_response_avgblock4B;
diffac=all_response_avgblock4A-all_response_avgblock3;
diffbc=all_response_avgblock4B-all_response_avgblock3;




conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
diff_cond={'difference A vs B','difference A vs C','difference B vs C'};
bin_vec=-3:0.5:3;
    subplot(1,3,1)
    hist([diffab(:,1),diffac(:,1),diffbc(:,1)],bin_vec)


    legend (diff_cond)
    title(conds{1})
    h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    
    subplot(1,3,2)
    hist([diffab(:,2),diffac(:,2),diffbc(:,2)],bin_vec)

    legend (diff_cond)
    title(conds{2})
        h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    subplot(1,3,3)
    
    hist([diffab(:,3),diffac(:,3),diffbc(:,3)],bin_vec)

    legend (diff_cond)
    title(conds{3})
    h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    
    
    %% bounded lines for learning -extinction - test
    mean_response_type=squeeze(nanmean(filt_avg_response_type(:,:,1:block_def(3,2))));
    std_response_type=squeeze(nanstd(filt_avg_response_type(:,:,1:block_def(3,2))))./sqrt(numel(selvps));
    
    figure
   hold on
   
rectangle('Position',[0 0 24 4])
rectangle('Position',[24 0 24 4])
rectangle('Position',[48 0 16 4])

for con=1:3
    
    x1=1:size(mean_response_type,2);
    y1=mean_response_type(con,:);
    b1=std_response_type(con,:);
boundedline(x1, y1, b1, 'cmap',cmap_default(con,:),'transparency',0.1,'alpha');

end
    h=gca;
    h.XLim=[1 block_def(3,2)];
    h.YLim=[1 4];
    h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
xlabel('Trial No')
legend('cs+/cs+','','cs+/cs-','','cs-/cs-','');    

%%
% learning
conds={'plusplus','plusminus','minusminus'};
export_spss=[squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,2,1)),squeeze(response_avgblocksub(:,3,1)),...
    squeeze(response_avgblocksub(:,1,2)),squeeze(response_avgblocksub(:,2,2)),squeeze(response_avgblocksub(:,3,2)),...
     squeeze(response_avgblocksub(:,1,3)),squeeze(response_avgblocksub(:,2,3)),squeeze(response_avgblocksub(:,3,3))];
names={strcat('learn1_',conds{1}),strcat('learn1_',conds{2}),strcat('learn1_',conds{3}),...
    strcat('learn2_',conds{1}),strcat('learn2_',conds{2}),strcat('learn2_',conds{3}),...
    strcat('test_',conds{1}),strcat('test_',conds{2}),strcat('test_',conds{3})}';



% test

    