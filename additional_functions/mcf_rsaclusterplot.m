
%rsa_mat % matrix sub x values in each sub cond

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
violinplot( rsa_mat,cat)
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


within = table(cat_mat, 'VariableNames', factors);
% fit the repeated measures model
rm = fitrm(rsa_val,'V1-V6~1','WithinDesign',within);
model=[];
for m=1:(numel(factors)-1)
model=[model,factors{m},'*'];
end
model=[model,factors{end}];

[rmtable] = ranova(rm, 'WithinModel','type*block');


subplot(1,2,2)
ylim([0 1])
axis off
hold on
rmtbl=rmtable;

% delete unneeded info from table

% intercept=strcmp('');
% error_rows=strcmp('Error');

% delete headings
% delete alternative p values etc

% add error df 


% loop through rows to write results as text

text(0,1,['anova: type x block'])
text(0,0.8,['type:','F(',num2str(ranovatblb.DF('(Intercept):type')),',',...
    num2str(ranovatblb.DF('Error(type)')),')=',num2str(ranovatblb.F('(Intercept):type')),...
    'p=',num2str(ranovatblb.pValue('(Intercept):type'))],'FontSize',8)
text(0,0.6,['block:','F(',num2str(ranovatblb.DF('(Intercept):block')),',',...
    num2str(ranovatblb.DF('Error(block)')),')=',num2str(ranovatblb.F('(Intercept):block')),...
    'p=',num2str(ranovatblb.pValue('(Intercept):block'))],'FontSize',8)
text(0,0.4,['typexblock:','F(',num2str(ranovatblb.DF('(Intercept):type:block')),',',...
    num2str(ranovatblb.DF('Error(type:block)')),')=',num2str(ranovatblb.F('(Intercept):type:block')),...
    'p=',num2str(ranovatblb.pValue('(Intercept):type:block'))],'FontSize',8)

clear varNames factorNames
