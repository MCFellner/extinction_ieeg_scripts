%% plot sig electrodes on surface brain

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_stat='D:\Extinction\iEEG\analysis\erp\';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};
 
conditions={'A_us','A_nous'};
sel_window=[4 4.5]; % add rows for more windows
alpha_def=0.05;


for sub=1:numel(allsubs)
    sel_sub=allsubs{sub};
 info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)  
folder_in=fullfile(path_stat,strcat(conditions{1},'_vs_',conditions{1},'_toi',num2str(sel_window(1)*1000),'to',num2str(sel_window(2)*1000),'sec'));
load(fullfile(folder_in,strcat(sel_sub,'erpstat')))
erpstat=stat;
sig_def{sub}=stat.prob_pos<alpha_def|stat.prob_neg<alpha_def;
all_elec{sub}=stat.label;
        
% get positions        
 [~,~,ind]=intersect(all_elec{sub},datainfo.elec_info.bipolar.elec_mni.label,'stable');     
all_elec_pos{sub}=datainfo.elec_info.bipolar.elec_mni.elecpos(ind,:);
all_elec_label{sub}=datainfo.elec_info.bipolar.ana_labels.freesurferDK(ind,1);
end

% get labels and relative count in each area
all_label=[];
all_sig=[];
for n=1:numel(all_elec_label)
all_label=[all_label;all_elec_label{n}];
all_sig=[all_sig;sig_def{n}'];
end
all_label=[all_label{:}];
all_regions=unique(all_label)';

for r=1:numel(all_regions)
sel_region=all_regions{r};
ind= strcmp(all_label,sel_region);
rel_sig{r,1}=sum(all_sig(ind))./sum(ind);  
absolute_num{r,1}= sum(ind);
sig_num{r,1}= sum(all_sig(ind));
end

region_count=table(all_regions,rel_sig,absolute_num, sig_num,'VariableNames',{'analabel','relsigelec','absnumelecinregion','absnumsigelec',})


save(fullfile(folder_in,'results_table'),'region_count')
%%

% plot electrodes and count electrodes per region/pat

% elec to plot
elec_toplot.unit ='mm';
elec_toplot.coordsys ='mni';
load('D:\Extinction\iEEG\scripts\additional_functions\sel_colorseries.mat')

sig_ind=[];
all_label=[];
all_pos=[];
for n=1:numel(all_elec)
   sig_ind=[sig_ind;sig_def{n}'];
   all_label=[all_label;all_elec{n}];
   all_pos=[all_pos;all_elec_pos{n}];
end


views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90;90 -40;];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90;-90 -40];
mesh.coordsys = 'mni';
hemispheres={'left','right'};
elec_def=[-1,1];
type_elec=[1,0];
sel_col=[1,0,0;0.5 0.5 0.5];
trans_sphere=[1 0.2];
for h=1:numel(hemispheres)
    sel_hemi=hemispheres{h};
    sel_elec_def=elec_def(h);
    load(fullfile('D:\matlab_tools\fieldtrip-20200130\template\anatomy',strcat('surface_pial_',sel_hemi,'.mat')));
    figure
    ft_plot_mesh(mesh,'facealpha',0.2,  'edgealpha',0.2);
    hold on
    % elec to plot
    elec_toplot.unit ='mm';
    elec_toplot.coordsys ='mni';
    % all elecs in selected hemisphere
   sel_elec_hemi=sign(all_pos(:,1))==sel_elec_def;

    for t=1:numel(type_elec)
    sel_elec=sel_elec_hemi&(sig_ind==type_elec(t));
    sel_trans=trans_sphere(t);
    elec_toplot.label=all_label(sel_elec,:);
    elec_toplot.elecpos=all_pos(sel_elec,:);
    elec_toplot.chanpos=all_pos(sel_elec,:);
    sel_color=sel_col(t,:);
    
    
    ft_plot_sens(elec_toplot,'elec','true','elecshape','sphere','facecolor',sel_color,'facealpha',sel_trans);
    end
    view([-90 20]);
    material dull;
    view(squeeze(views(h,1,:))');
    c1=camlight(0,0);
    set(c1, 'style', 'infinite');

    view(squeeze(views(h,2,:))');
    c2=camlight(0, 0);
    set(c2, 'style', 'infinite');

    view(squeeze(views(h,3,:))');
    print('-f1','-r600','-dtiff',fullfile(folder_in,strcat(sel_hemi,'_lat.tiff'))) 

    view(squeeze(views(h,4,:))');
    print('-f1','-r600','-dtiff',fullfile(folder_in,strcat(sel_hemi,'_med.tiff'))) 
    clear c1 c2 
    close all
end
sum(sig_ind)./numel(sig_ind)

