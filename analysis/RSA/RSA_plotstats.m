% stats plotter

%plot results of rsa stats in organized plots

% file definition
path_stats='D:\Extinction\iEEG\analysis\rsa\powlogscale_timeslide_z_crosstrials_toi2000to4000\stats';

% adapt file definition to load the wanted files
contrast='*'; %* if irrelevant
roi='amy_l'; %* if irrelevant

% define whether to plot multiplots or single plots
multiplot='yes'

load('D:\matlab_tools\jet_grey.mat')

if strcmp(contrast,'*')&strcmp(roi,'*')
    all_files=dir(fullfile(path_stats,[contrast,'\fig\',roi,'.mat']));
    
else
    all_files=dir(fullfile(path_stats,[contrast,'\fig\*',roi,'.mat']));
end

switch multiplot
    case 'no'
        for f=1:numel(all_files)
            sel_path=all_files(f).folder;
            sel_file=all_files(f).name;
            load(fullfile(sel_path,sel_file))
            % find seperator in str
            sep_ind=strfind(sel_file,'_in_');
            sel_roi=sel_file(sep_ind+4:end-4);
            sel_contrast=sel_file(1:sep_ind-1);
            fig=figure
            imagesc(stats.time,stats.time,squeeze(stats.stat),[-5 5])
            hold on
            colormap(jet_grey)
            colorbar
            title({[sel_contrast];[sel_roi];['pos tsum:',num2str(stats.trial_rand.data_pos(1)),'p=',num2str(stats.trial_rand.p_pos(1))];['neg tsum:',num2str(stats.trial_rand.data_neg(1)),'p=',num2str(stats.trial_rand.p_neg(1))]})
            
            ylabel('t in s')
            xlabel('t in s')
            contour(stats.time,stats.time,squeeze(stats.trial_rand.mask),1,'k')
            set(gca,'YDir','normal')
            savefig(fig,fullfile(sel_path,[sel_file,'.fig']),'compact')
            close all
        end
        
        
    case 'yes'
        fig=figure
        num_col=ceil(numel(all_files)./3);
        if strcmp(contrast,'*')
            sgtitle(roi)
        elseif strcmp(roi,'*')
            sgtitle(contrast)
        end
        
        for f=1:numel(all_files)
            sel_path=all_files(f).folder;
            sel_file=all_files(f).name;
            load(fullfile(sel_path,sel_file))
            % find seperator in str
            sep_ind=strfind(sel_file,'_in_');
            sel_roi=strrep(sel_file(sep_ind+4:end-4),'_',' ');
            sel_contrast=strrep(sel_file(1:sep_ind-1),'_',' ');
            sel_contrast=strrep(sel_contrast,'mask','x');
             sel_contrast=strrep(sel_contrast,'type','');

            % breakup the long string for interactions
            int_ind=strfind(sel_contrast,'interaction');
            if ~isempty(int_ind)
                tmp_contrast={[sel_contrast(1:int_ind-1),' x '];[sel_contrast(int_ind+11:end)]};
            else
                tmp_contrast{1}=sel_contrast;
                tmp_contrast{2}=' ';
            end
            sel_contrast=tmp_contrast;
            
            subplot(3,num_col,f)
            imagesc(stats.time,stats.time,squeeze(stats.stat),[-5 5])
            hold on
            colormap(jet_grey)
            colorbar
            if strcmp(contrast,'*')
                
                title({sel_contrast{1};sel_contrast{2};...
                    ['pos tsum:',num2str(stats.trial_rand.data_pos(1)),'p=',num2str(stats.trial_rand.p_pos(1))];...
                    ['neg tsum:',num2str(stats.trial_rand.data_neg(1)),'p=',num2str(stats.trial_rand.p_neg(1))]},...
                    'FontSize',8,'FontWeight','normal')
                
                path_fig=fullfile(path_stats,'summary_fig',roi)
                mkdir(path_fig)
            elseif strcmp(roi,'*')
                title({[sel_roi];['pos tsum:',num2str(stats.trial_rand.data_pos(1)),'p=',num2str(stats.trial_rand.p_pos(1))];['neg tsum:',num2str(stats.trial_rand.data_neg(1)),'p=',num2str(stats.trial_rand.p_neg(1))]})
                path_fig=fullfile(path_stats,'summary_fig',contrast)
                mkdir(path_fig)
            end
            ylabel('t in s')
            xlabel('t in s')
            contour(stats.time,stats.time,squeeze(stats.trial_rand.mask),1,'k')
            set(gca,'YDir','normal')
            
        end
        savefig(fig,fullfile(path_fig,['summary.fig']),'compact')
        
        close all
        
end
