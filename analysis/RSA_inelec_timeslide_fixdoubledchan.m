%% fix falsely saved channels



path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_rsa='D:\Extinction\iEEG\analysis\rsa\';
path_designmat='D:\Extinction\iEEG\analysis\rsa\contrast_mat\';


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
    'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
    'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub08'};




contrasts={'item_specific','item_specific_block1','item_specific_block2','item_specific_block3',...
    'cs_specific','cs_specific_block1','cs_specific_block2',...
    'type1to2_vs_type2to3_block1','type1to2_vs_type2to3_block2'};
feature='powlogscale';
norm='z_crosstrials';
toi=[2 4];
win_pow=0.05;
for c=1:numel(contrasts)
    contrast=contrasts{c};
    folder_in=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),'stats',contrast);
    path_out=fullfile(folder_in,'fig');
    mkdir(path_out)
    
    for sub=1:length(allsubs)
        sel_sub=allsubs{sub};
        sel_folder=fullfile(path_rsa,strcat(feature,'_timeslide_',norm,'_toi',num2str(toi(1)*1000),'to',num2str(toi(2)*1000)),sel_sub);
        % electrodeinfo
        info_file=strcat(path_info,sel_sub,'_datainfo');
        load(info_file)
        load(fullfile(folder_in,strcat(sel_sub,'_rsastat')))
        correct_chan=datainfo.elec_info.bipolar.montage_withoutartichan.labelnew;
        all_chan=1:numel(correct_chan);
        all_stat=all_stat(all_chan);
        save(fullfile(folder_in,strcat(sel_sub,'_rsastat')),'all_stat')
    end
end