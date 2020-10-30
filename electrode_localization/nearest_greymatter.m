addpath('D:\matlab_tools\fieldtrip-20200130')
ft_defaults
addpath('D:\Extinction\iEEG\scripts\additional_functions')

%% check for electrodes in wm, no label or ventricle for nearest grey matter

path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};

to_check={'no_label_found','ctx-rh-unknown',...
    'WM-hypointensities','Right-Cerebral-White-Matter','Left-Cerebral-White-Matter',...
    'Right-choroid-plexus','Left-choroid-plexus','Right-Inf-Lat-Vent','Right-Lateral-Ventricle','Left-Lateral-Ventricle','Left-Inf-Lat-Vent'}

rad=2;% maximalvoxel rad to check for grey matter
all_ref={'bipolar','no_reref'};
for r=1:numel(all_ref)
    ref=all_ref{r};
    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
        
        switch ref
            case 'bipolar'
                elecinfo=datainfo.elec_info.bipolar;
            case 'no_reref'
                elecinfo=datainfo.elec_info;
        end
        max_col=nearest(elecinfo.ana_labels.freesurferDK_def.sphereradius_mm,rad);
        
        % check ana_labels for to_check labels, for these labels check up to
        % max col for a no to_check label        
        sel_labels=[elecinfo.ana_labels.freesurferDK{:,1}]';
        sel_next_labels=elecinfo.ana_labels.freesurferDK(:,1:max_col);
        for lab=1:numel(sel_labels)
            sel_label=sel_labels{lab};
            check=any(strcmp(sel_label,to_check));
            if check
                % select first no_check label (labels are ordered for distance)
                tmp_labels= [sel_next_labels{lab,2:max_col}];
                [~,ind]=setdiff(tmp_labels,to_check,'stable');
                if isempty(ind)
                    new_labels{lab,1}=sel_label;
                else
                    new_labels{lab,1}=tmp_labels(ind(1));
                    
                end
            else
                new_labels{lab,1}=sel_label;
                
            end
        end
               
        switch ref
            case 'bipolar'
                datainfo.elec_info.bipolar.ana_labels.nearestGMlabelfreesurferDK=new_labels;
            case 'no_reref'
                datainfo.elec_info.ana_labels.nearestGMlabelfreesurferDK=new_labels;
        end
        save(info_file,'datainfo')
        clear new_labels sel_labels sel_next_labels tmp_labels ind sel_label
    end
end
