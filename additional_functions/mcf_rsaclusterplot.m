
%rsa_mat % matrix sub x values in each sub cond
% factors: cell array with factor names
%cat_mat % cell, each column labels of a factor
% colorscheme=repmat([1 0 0; 0 0 1]) % vector

function [rmtable]=mcf_rsaclusterplot(rsa_mat,factors,cat_mat,colorscheme)


co = colorscheme;
set(groot,'defaultAxesColorOrder',co)
cmap_default=co;

num_cond=size(rsa_mat,2)
for i=1:num_cond
cat{i}=[cat_mat{i,:}];
end

subplot(1,2,1)
violinplot( rsa_mat,cat,'ViolinAlpha',0.5)
% bar(reshape(squeeze(nanmean(all_contrast_rsa.item_specific_x_block))',1,[]))
% hold on
% scatter(reshape(repmat([1 2 3 4 5 6],numel(sel_subs),1),1,[]),...
% reshape(permute(all_contrast_rsa.item_specific_x_block,[1,3,2]),1,[]))
% xticks(1:6)
% xticklabels({'wi block1','bi block1','wi block2','bi block2','wi block3','bi block3'})
xtickangle(45)
set(groot,'defaultAxesColorOrder','remove')

for i = 1 :   num_cond
    v = strcat('V',num2str(i));    
    varNames{i,1} = v;
end
% Create a table storing the respones
rsa_val = array2table(rsa_mat, 'VariableNames',varNames);
within = array2table(cat_mat, 'VariableNames', factors);
% fit the repeated measures model
rm = fitrm(rsa_val,[varNames{1},'-',varNames{end},'~1'],'WithinDesign',within);
model=[];
for m=1:(numel(factors)-1)
model=[model,factors{m},'*'];
end
model=[model,factors{end}];

[rmtable] = ranova(rm, 'WithinModel',model);

subplot(1,2,2)
ylim([0 1])
axis off
hold on
rmtbl=table2cell(rmtable);

% delete unneeded info from table
rownames=rmtable.Properties.RowNames;
varnames=rmtable.Properties.VariableNames;
%first two rows
rmtbl(1:2,:)=[];
rownames(1:2)=[];
% delete sums, alternative p values etc
varnames=varnames([2,4,5]);
rmtbl=rmtbl(:,[2,4,5]);
% add error df 
error_ind=strncmp(rownames,'Error',5);
error_df=rmtbl(error_ind,1);
rmtbl=rmtbl(error_ind==0,:);
rownames=rownames(error_ind==0);
rownames=strrep(rownames,'(Intercept):','')
rmtbl=[rmtbl(:,1),error_df,rmtbl(:,2:3)];
% convert table to str
rmtbl=cellfun(@num2str,rmtbl,'UniformOutput',0);

% loop through rows to write results as text
num_con=size(rmtbl,1);
spc=1/num_con;
for i=1:num_con
   text(0,1-(spc*i),[rownames{i},': F(',rmtbl{i,1},',',rmtbl{i,2},')=',rmtbl{i,3},' p=',rmtbl{i,4}],'FontSize',8) 
end

