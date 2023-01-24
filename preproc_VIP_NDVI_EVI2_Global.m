function [] = preproc_VIP_NDVI_EVI2_Global(pp_type)
    % pp_type: RS_ORG, RS_05
    %% VIP15_v004, EVI2 and NDVI data
    global data_src;
    global domain_def;
    global chosen_prd;
    global DATA_VI_05rs_out;
    global DATA_VI_05rs_sm_out;
    global DATA_VI_05rs_bw_out;
    global DATA_VI_05rs_bw_sm_out;
    global nds nmd;
    
    % Datasource ID
    ds_id=2;

    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'VIP', pp_type);
    
    dtype={'EVI2'}; % EVI2, NDVI
    vi_scalar=10000;

    % period for complete annual data
    % yr_e=1981;
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);

    if strcmp(pp_type,'RS_ORG')
        clearvars DATA_VIP_NDVI_bw_out DATA_VIP_EVI2_bw_out DATA_VIP_NDVI_out DATA_VIP_EVI2_out;
        if strcmp(dtype,'EVI2')
            % Bi-weekly data
            DATA_VIP_EVI2_bw_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1);
            % Monthly data
            DATA_VIP_EVI2_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1);
        end
        if strcmp(dtype,'NDVI')
            % Bi-weekly data
            DATA_VIP_NDVI_bw_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1);
            % Monthly data
            DATA_VIP_NDVI_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1);
        end

        import matlab.io.hdfeos.*

        for y=yr_s:yr_e
            % ---- Read data ---
            % mmdd: dd changes by month, therefore not fixed date format can be
            % used

            yi=y-yr_s+1;
            wi=0;
            for mm=1:12
                for dd=1:31

                    date_fmt=sprintf('%d.%s.%s',y, num2str(mm,'%02d'), num2str(dd,'%02d'));
                    dir_in=sprintf('%s/VIP/VIP15.004/%s/',data_src, date_fmt);
                    fname=dir(sprintf('%s/*.hdf',dir_in));

                    if size(fname,1)==1

                        fprintf(sprintf('--> Processing file : %s\n',fname.name));

                        wi=wi+1;
                        fullname=sprintf('%s/%s',dir_in, fname.name);
                        gfid = gd.open(fullname,'read');

                        gridID = gd.attach(gfid,'VIP_CMG_GRID');
                        cornerlat = [domain_def(3) domain_def(4)];
                        cornerlon = [domain_def(1) domain_def(2)];
                        regionID = gd.defBoxRegion(gridID,cornerlat,cornerlon);

                        % Field name can be found from hdfinfo (e.g. finfo.Vgroup.Vgroup(1).SDS(2).Name)
                        % data valid range: [-1*vi_scalar 1*vi_scalar], for both EVI2 and NDVI
                        if strcmp(dtype,'EVI2')
                            dtmp = double(gd.extractRegion(gridID,regionID, 'CMG 0.05 Deg 15DAYS EVI2'));
                            dtmp(dtmp<-vi_scalar)=nan; dtmp(dtmp>vi_scalar)=nan;
                            DATA_VIP_EVI2_bw_out(:,:,wi,yi)=dtmp;
                        end

                        if strcmp(dtype,'NDVI')
                            dtmp = double(gd.extractRegion(gridID,regionID, 'CMG 0.05 Deg 15DAYS NDVI'));
                            dtmp(dtmp<-vi_scalar)=nan; dtmp(dtmp>vi_scalar)=nan;
                            DATA_VIP_NDVI_bw_out(:,:,wi,yi)=dtmp;
                        end

                        gd.detach(gridID);
                        gd.close(gfid);
                    end
                end %% dd
                % Aggregrated to monthly values
                if strcmp(dtype,'EVI2')
                    DATA_VIP_EVI2_out(:,:,mm,yi)=squeeze(nanmean(DATA_VIP_EVI2_bw_out(:,:,(mm-1)*2+1:mm*2, yi),3));
                end
                if strcmp(dtype,'NDVI')
                    DATA_VIP_NDVI_out(:,:,mm,yi)=squeeze(nanmean(DATA_VIP_NDVI_bw_out(:,:,(mm-1)*2+1:mm*2, yi),3));
                end

            end %% mm
            clearvars dtmp date_fmt dir_in fname;
        end

        % 10000 is a scalar for original data;
        if strcmp(dtype,'EVI2')
            DATA_VIP_EVI2_bw_out=DATA_VIP_EVI2_bw_out/double(vi_scalar);
            DATA_VIP_EVI2_out=DATA_VIP_EVI2_out/double(vi_scalar);
        end

        if strcmp(dtype,'NDVI')
            DATA_VIP_NDVI_bw_out=DATA_VIP_NDVI_bw_out/double(vi_scalar);
            DATA_VIP_NDVI_out=DATA_VIP_NDVI_out/double(vi_scalar);
        end
    end % RS_ORG
    
    if strcmp(pp_type,'RS_05')
        
        % check if the array DATA_VI_05rs has created
        if sum(size(DATA_VI_05rs_out)) == 0 || sum(size(DATA_VI_05rs_sm_out)) == 0 ...
            || sum(size(DATA_VI_05rs_bw_out)) == 0 || sum(size(DATA_VI_05rs_bw_sm_out)) == 0
            
            % Create 
            % Monthly data          nlon         nlat   nweek/nmonth  nyear   ndatasource
            DATA_VI_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1,  nds);
            % Monthly data smoothed
            DATA_VI_05rs_sm_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1,  nds, nmd);
            
            % Bi-weekly data
            DATA_VI_05rs_bw_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1, nds);
            % Bi-weekly data smoothed
            DATA_VI_05rs_bw_sm_out = nan(dmn_lon_n_g, dmn_lat_n_g, 24, (yr_e-yr_s)+1, nds, nmd);
            
            fprintf(('Arrays created! \n'));
        end
        
        for y=yr_s:yr_e
            % ---- Read data ---
            % mmdd: dd changes by month, therefore not fixed date format can be
            % used

            yi=y-yr_s+1;
            wi=0;
            
            for mm=1:12
                for dd=1:31

                    date_fmt=sprintf('%d.%s.%s',y, num2str(mm,'%02d'), num2str(dd,'%02d'));
                    dir_in=sprintf('%s/VIP/VIP15.004/netcdf/0.5deg/',data_src);
                    fname=dir(sprintf('%s/VIP15.EVI2_NDVI_%s_0.5deg_remapbil.nc',dir_in, date_fmt));

                    if size(fname,1)==1

                        fprintf(sprintf('--> Processing file : %s\n',fname.name));

                        wi=wi+1;
                        fullname=sprintf('%s/%s',dir_in, fname.name);
                        nc_var = ncgeodataset(fullname);
                        
                        % ---- Vegetation indices  ---
                        % EVI2 (lat, lon)
                        dtmp=squeeze(double(nc_var.data('EVI2')));
                        dtmp(dtmp<-vi_scalar)=nan; dtmp(dtmp>vi_scalar)=nan;
                        
                        % Converstion with scalar has been applied when
                        DATA_VI_05rs_bw_out(:,:,wi, yi, ds_id)=squeeze(dtmp(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
                    end
                        
                end %% dd
                
                % Aggregrated to monthly values
                DATA_VI_05rs_out(:,:,mm,yi, ds_id)=squeeze(nanmean(DATA_VI_05rs_bw_out(:,:,(mm-1)*2+1:mm*2, yi, ds_id),3));
            end %% mm
            
            clearvars dtmp date_fmt dir_in fname;
        end
        
        
        [s1 s2 s3 s4 s5]=size(DATA_VI_05rs_bw_out);
        data_bw_rs=reshape(DATA_VI_05rs_bw_out, [s1 s2 s3*s4 s5]);

        % Smooth the bi-weekly data based on different methods
        w=7; % windows size
        ts_sm1=smoothdata(data_bw_rs(:,:,:,ds_id),3,'sgolay',w, 'Degree',2);
        ts_sm2=smoothdata(data_bw_rs(:,:,:,ds_id),3,'movmean',w);
        ts_sm3=smoothdata(data_bw_rs(:,:,:,ds_id),3,'gaussian',w);

        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,1)=reshape(ts_sm1, [s1 s2 s3 s4]);
        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,2)=reshape(ts_sm2, [s1 s2 s3 s4]);
        DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,3)=reshape(ts_sm3, [s1 s2 s3 s4]);
        % mean of the smoothed values from different methods
%         DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,4)=nanmean(DATA_VI_05rs_bw_sm_out(:,:,:,:,ds_id,1:3),6);            

        % Aggregate to monthly
        for mm=1:12
            DATA_VI_05rs_sm_out(:,:,mm,:,ds_id,1:nmd)=squeeze(nanmean(DATA_VI_05rs_bw_sm_out(:,:,(mm*2-1):mm*2,:,ds_id,1:nmd),3));
        end
    end
end







