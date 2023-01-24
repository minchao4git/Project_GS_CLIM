function [] = preproc_GIMMS_NDVI3g_Global(pp_type)
    %% NDVI3g
    global data_src;
    global domain_def;
    global chosen_prd;
    global DATA_VI_05rs_out DATA_VI_05rs_sm_out;
    global DATA_VI_05rs_bw_out DATA_VI_05rs_bw_sm_out;
    global nds nmd;
    
    % Datasource ID
    ds_id=1;

    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'NDVI3g', pp_type);
    vi_scalar=10000;

    % period for complete annual data
    % yr_e=1981;
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    
    if strcmp(pp_type,'RS_ORG')
        
        dir_in=sprintf('%s/GIMMS/NDVI3g/0.083deg',data_src);

        % Bi-weekly data
        DATA_NDVI3g_bw_out = nan(24, dmn_lon_n_g, dmn_lat_n_g, (yr_e-yr_s)+1);
        % will change back to format(dmn_lon_n_g, dmn_lat_n_g, nmon, nyr) in the
        % end
        % Monthly data
        DATA_NDVI3g_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1);

        clearvars DATA_NDVI3g_out DATA_NDVI3g_bw_out;
        for y=yr_s:yr_e
            % ---- Read data ---
            % mmdd either 0106: the first half year 
            %          or 0712: the second half year

            % --> I: read in first half year
            mmdd='0106';
            file_name=sprintf('ndvi3g_geo_v1_%4d_%s.nc4',y,mmdd);
            nc_var = ncgeodataset(sprintf('%s/%s', dir_in, file_name));
            fprintf(sprintf('--> Processing file : %s\n',file_name));

            % ---- Vegetation indices  ---
            % ndvi (time, lat, lon)
            data_tmp=squeeze(double(nc_var.data('ndvi')));
            data_tmp(data_tmp<double(0.0))=nan;
            data_tmp=squeeze(data_tmp(:, y_s_e_g(1):y_s_e_g(2), x_s_e_g(1):x_s_e_g(2)));
            for bw=1:12
                DATA_NDVI3g_bw_out(bw,:,:,(y-yr_s)+1)=squeeze(data_tmp(bw,:,:))';
            end
            clearvars nc_var data_tmp;

            % --> II: read in first half year
            mmdd='0712';
            file_name=sprintf('ndvi3g_geo_v1_%4d_%s.nc4',y,mmdd);
            nc_var = ncgeodataset(sprintf('%s/%s', dir_in, file_name));
            fprintf(sprintf('--> Processing file : %s\n',file_name));

            % ---- Vegetation indices  ---
            % ndvi (time, lat, lon)
            data_tmp=squeeze(double(nc_var.data('ndvi')));
            data_tmp(data_tmp<double(0.0))=nan;
            data_tmp=squeeze(data_tmp(:, y_s_e_g(1):y_s_e_g(2), x_s_e_g(1):x_s_e_g(2)));
            for bw=13:24
                DATA_NDVI3g_bw_out(bw,:,:,(y-yr_s)+1)=squeeze(data_tmp(bw-12,:,:))';
            end
            clearvars nc_var data_tmp;
        end

        DATA_NDVI3g_bw_out=permute(DATA_NDVI3g_bw_out,[2 3 1 4]);
        
        for m=1:12
            DATA_NDVI3g_out(:,:,m,:)=squeeze(nanmean(DATA_NDVI3g_bw_out(:,:,(m*2-1):m*2,:),3));
        end
        DATA_NDVI3g_bw_out=DATA_NDVI3g_bw_out/double(10000); % 10000 is a scalar for original data;
        DATA_NDVI3g_out=DATA_NDVI3g_out/double(10000);
    end % RS_ORG
    
    if strcmp(pp_type,'RS_05')
        
        dir_in=sprintf('%s/GIMMS/NDVI3g/0.5deg',data_src);
        
        % check if the array DATA_VI_05rs has created
        if size(DATA_VI_05rs_out,2) <= 1
            % Create 
            % Monthly data          nlon         nlat   nweek/nmonth  nyear   ndatasource
            DATA_VI_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1,  nds);
            % Monthly data smoothed
            DATA_VI_05rs_sm_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1,  nds, nmd);
            
            % Bi-weekly data
            DATA_VI_05rs_bw_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1, nds);
            % Bi-weekly data smoothed
            DATA_VI_05rs_bw_sm_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1, nds, nmd);
        end
        
        tmpdata_bw = nan(24, dmn_lon_n_g, dmn_lat_n_g, (yr_e-yr_s)+1);
        
        for y=yr_s:yr_e
            % ---- Read data ---
            % mmdd either 0106: the first half year 
            %          or 0712: the second half year

            % --> I: read in first half year
            mmdd='0106';
            file_name=sprintf('ndvi3g_geo_v1_%4d_%s_0.5deg_remapbil.nc',y,mmdd);
            nc_var = ncgeodataset(sprintf('%s/%s', dir_in, file_name));
            fprintf(sprintf('--> Processing file : %s\n',file_name));

            % ---- Vegetation indices  ---
            % ndvi (time, lat, lon)
            data_tmp=squeeze(double(nc_var.data('ndvi')));
            data_tmp(data_tmp<double(0.0))=nan;
            data_tmp=squeeze(data_tmp(:, y_s_e_g(1):y_s_e_g(2), x_s_e_g(1):x_s_e_g(2)));
            for bw=1:12
                tmpdata_bw(bw,:,:,(y-yr_s)+1)=squeeze(data_tmp(bw,:,:))';
            end
            clearvars nc_var data_tmp;

            % --> II: read in first half year
            mmdd='0712';
            file_name=sprintf('ndvi3g_geo_v1_%4d_%s_0.5deg_remapbil.nc',y,mmdd);
            nc_var = ncgeodataset(sprintf('%s/%s', dir_in, file_name));
            fprintf(sprintf('--> Processing file : %s\n',file_name));

            % ---- Vegetation indices  ---
            % ndvi (time, lat, lon)
            data_tmp=squeeze(double(nc_var.data('ndvi')));
            data_tmp(data_tmp<double(0.0))=nan;
            data_tmp=squeeze(data_tmp(:, y_s_e_g(1):y_s_e_g(2), x_s_e_g(1):x_s_e_g(2)));
            for bw=13:24
                tmpdata_bw(bw,:,:,(y-yr_s)+1)=squeeze(data_tmp(bw-12,:,:))';
            end
            clearvars nc_var data_tmp;
        end
        
        ts_bw_org=permute(tmpdata_bw,[2 3 1 4])/double(10000); % 10000 is a scalar for original data;
        [s1 s2 s3 s4]=size(ts_bw_org);
        
        ts_bw_org_rs=reshape(ts_bw_org,[s1 s2 s3*s4]);
        
        % Smooth the bi-weekly data based on different methods
        w=7; % windows size
        ts_sm1=smoothdata(ts_bw_org_rs,3,'sgolay',w, 'Degree',2);
        ts_sm2=smoothdata(ts_bw_org_rs,3,'movmean',w);
        ts_sm3=smoothdata(ts_bw_org_rs,3,'gaussian',w);
        
        DATA_VI_05rs_bw_out(:,:,:,:,ds_id)=ts_bw_org;
        
        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,1)=reshape(ts_sm1, [s1 s2 s3 s4]);
        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,2)=reshape(ts_sm2, [s1 s2 s3 s4]);
        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,3)=reshape(ts_sm3, [s1 s2 s3 s4]);
        % mean of the smoothed values from different methods
        
        % Aggregate to monthly
        for m=1:12
            DATA_VI_05rs_out(:,:,m,:,ds_id)=squeeze(nanmean(DATA_VI_05rs_bw_out(:,:,(m*2-1):m*2,:,ds_id),3));
            DATA_VI_05rs_sm_out(:,:,m,:,ds_id,1:nmd)=squeeze(nanmean(DATA_VI_05rs_bw_sm_out(:,:,(m*2-1):m*2,:,ds_id,1:nmd),3));
        end
        
    end % RS_05
    
end


