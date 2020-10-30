%% build montages for rereferencing
%check data avg reference, bipolar, or nearest white matter
% average reference: use ft_preprocessing
path_info='D:\Extinction\iEEG\data\preproc\ieeg\datainfo\';
path_preproc='D:\Extinction\iEEG\data\preproc\ieeg\readin\';

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%     'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18','c_sub20',...
%     'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07'};
allsubs = {'c_sub20'};
% load datainfo for each subject
for n=1:numel(allsubs)
    sel_sub=allsubs{n};
    
    % electrodeinfo
    info_file=strcat(path_info,sel_sub,'_datainfo');
    load(info_file)
    
    % datainfo.elec_info.reordered.cell_ieeg: defines every electrode
    % get for every electrode a rereference scheme
    elecs=unique(datainfo.elec_info.reordered.cell_ieeg(:,1));
    
    %%%%%bipolar
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        sel_chan=datainfo.elec_info.reordered.sorted_ieeg(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1)));
        sel_num=[datainfo.elec_info.reordered.cell_ieeg{(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1))),2}];
        %check for missing number
        if any(diff(sel_num)~=1)
            bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
            bip_elec{e}(diff(sel_num)~=1,:)=[];
        else
            bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
        end
        bip_elec{e}=[sel_chan(1:end-1),sel_chan(2:end)];
    end
    
    bip_reref=vertcat(bip_elec{:});
    bip_labelold=datainfo.elec_info.reordered.sorted_ieeg;
    
    check=setdiff(bip_labelold,datainfo.elec_info.orderindata.label);
    if ~isempty(check)
        error('labels in elec_info.orderindata and reordered do not match')
    end
    
    for chan=1:size(bip_reref,1)
        % get new name
        bip_labelnew{chan}=strcat(bip_reref{chan,1},'-',bip_reref{chan,2});
        % define tra mat (1,-1)
        tra(chan,:)=strcmp(bip_reref{chan,1},bip_labelold)-strcmp(bip_reref{chan,2},bip_labelold);
        % recalculate elec pos
        elec1_ind=strcmp(bip_reref{chan,1},datainfo.elec_info.elec_ct_mr.label);
        elec2_ind=strcmp(bip_reref{chan,2},datainfo.elec_info.elec_ct_mr.label);
        bip_pos_ct_mr(chan,:)=(datainfo.elec_info.elec_ct_mr.elecpos(elec1_ind,:)+datainfo.elec_info.elec_ct_mr.elecpos(elec2_ind,:)).*0.5;
        
        elec1_ind=strcmp(bip_reref{chan,1},datainfo.elec_info.elec_mni.label);
        elec2_ind=strcmp(bip_reref{chan,2},datainfo.elec_info.elec_mni.label);
        bip_pos_mni(chan,:)=(datainfo.elec_info.elec_mni.elecpos(elec1_ind,:)+datainfo.elec_info.elec_mni.elecpos(elec2_ind,:)).*0.5;
    end
    
    % add new info to datainfo
    datainfo.elec_info.bipolar.elec_ct_mr.label=bip_labelnew;
    datainfo.elec_info.bipolar.elec_ct_mr.chanpos=bip_pos_ct_mr;
    datainfo.elec_info.bipolar.elec_ct_mr.elecpos=bip_pos_ct_mr;
    datainfo.elec_info.bipolar.elec_ct_mr.unit='mm';
    datainfo.elec_info.bipolar.elec_ct_mr.coordsys='acpc';
    datainfo.elec_info.bipolar.elec_ct_mr.tra=diag(ones(1,numel(bip_labelnew)));
    
    datainfo.elec_info.bipolar.elec_mni.label=bip_labelnew;
    datainfo.elec_info.bipolar.elec_mni.chanpos=bip_pos_mni;
    datainfo.elec_info.bipolar.elec_mni.elecpos=bip_pos_mni;
    datainfo.elec_info.bipolar.elec_mni.unit='mm';
    datainfo.elec_info.bipolar.elec_mni.coordsys='mni';
    datainfo.elec_info.bipolar.elec_mni.tra=diag(ones(1,numel(bip_labelnew)));
    
    % build montage for each way of rereferncing
    %    montage.tra      = MxN matrix
    %    montage.labelold = Nx1 cell-array
    %    montage.labelnew = Mx1 cell-array
    
    datainfo.elec_info.bipolar.montage.tra=tra;
    datainfo.elec_info.bipolar.montage.labelold=bip_labelold;
    datainfo.elec_info.bipolar.montage.labelnew=bip_labelnew;
    
    % check bipolar montage
    %   [data_bip]    = ft_apply_montage(data,datainfo.elec_info.bipolar.montage)
     
    %%%%%% white matter contact
    % load data
    load(strcat(path_preproc,sel_sub,'_data.mat'))
    
    % check for each electrode for white matter contacts (in ct_mr)
    for e=1:numel(elecs)
        sel_elec=elecs{e};
        datainfo.elec_info.whitematter.electrodes=elecs;
        sel_chan=datainfo.elec_info.reordered.sorted_ieeg(strcmp(sel_elec,datainfo.elec_info.reordered.cell_ieeg(:,1)));
        datainfo.elec_info.whitematter.sel_chan_per_electrodes=sel_chan;

        for chan=1:numel(sel_chan)
            % check for white matter label 'White-Matter'
            all_labels=datainfo.elec_info.ana_labels.freesurferDK(strcmp(sel_chan{chan},datainfo.elec_info.ana_labels.labels),:);
            for l=1:numel(all_labels)
                sel_labels=all_labels{l};
                % get numbers of found regions for each
                num_labels(chan,l)=numel(sel_labels);
                check=(strfind(sel_labels, 'White-Matter'));
                wm_labels(chan,l)=~isempty([check{:}]);
            end
        end
       relation_WM2other_freesurferDK=wm_labels./num_labels;
        datainfo.elec_info.whitematter.relation_WM2other_freesurferDK{e}=relation_WM2other_freesurferDK; 
        WM_score_freesurferDK=cumsum(relation_WM2other_freesurferDK==1,2,'reverse');
        datainfo.elec_info.whitematter.WM_score_freesurferDK{e}=WM_score_freesurferDK;
        datainfo.elec_info.whitematter.possible_WM_reference{e}=sel_chan(max(WM_score_freesurferDK(:,1))==WM_score_freesurferDK(:,1));        
        clear num_labels wm_labels
        [~,ind]=intersect(sel_chan,data.label,'stable')
        %check for matching data
        if numel(ind)~=numel(sel_chan)
            error('datainfo channels and data channels do not match')
        end        
        sel_data=data.trial{1,1}(ind,:);
        datainfo.elec_info.whitematter.corr_mat{e}=corr(sel_data');
        datainfo.elec_info.whitematter.cov_mat{e}=cov(sel_data');
        datainfo.elec_info.whitematter.var_mat{e}=var(sel_data');
        datainfo.elec_info.whitematter.mean_mat{e}=mean(sel_data');
    end
    save(info_file,'datainfo')
  keep n allsubs path_info path_preproc
end
