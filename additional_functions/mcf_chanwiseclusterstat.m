% parfor cluster stats for RSA_inelec_timeslide

%cfg= cfg config for ft_clusterstats
%data1_mat3d = 3d matrix (trial*dim1*dim2)

function stat=mcf_chanwiseclusterstat(cfg,data1_mat3d,data2_mat3d, data_dummy)


            data1=data_dummy;
           data2=data_dummy;
             % put rsa_mat in ft freq structure
             data1.powspctrm=reshape(data1_mat3d,[size(data1_mat3d,1),1,size(data1_mat3d,2),size(data1_mat3d,3)]); 
             data2.powspctrm=reshape(data2_mat3d,[size(data2_mat3d,1),1,size(data2_mat3d,2),size(data2_mat3d,3)]);
%         
            stat=ft_freqstatistics (cfg,data1,data2)