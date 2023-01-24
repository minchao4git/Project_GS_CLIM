function [] = preproc_ERA5()
    
    vars_in={'t2m',      'ssrd',   'pr'};
    vars_conv={'- 273.15',   '',   '*days_of_mon(m)*1000'};
    %             1  2  3  4  5  6  7  8  9 10 11,12
    days_of_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
    
    % Local climate data from ERA5
    global data_src;
    global domain_def;
    global chosen_prd;
    global DATA_ERA5_05rs_out;
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'ERA5', 'RS_05');
    
    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
        
    % check if the array dataset has created
    if size(DATA_ERA5_05rs_out,2) <= 1
        % Create 
        % Monthly data:              nlon         nlat   nweek/nmonth  nyear   nvariable
        DATA_ERA5_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1, size(vars_in,2));
    end

    for v=1:size(vars_in,2)
        for y=yr_s:yr_e
            % ---- Read data ---
            yi=y-yr_s+1;

            dir_in=sprintf('%s/ECMWF/ERA5-Land/remap/0.5deg/%s',data_src,vars_in{v});
            fname=sprintf('%s_ECMWF-ERA5-Land_rean_mon_%d-%d_0.5deg_remapycon.nc',vars_in{v},y,y);

            fprintf(sprintf('--> Processing file : %s\n',fname));

            fullname=sprintf('%s/%s',dir_in, fname);
            nc_var = ncgeodataset(fullname);

            % ---- Climate variables  ---
            % vars (month, lat, lon)
            dtmp=squeeze(double(nc_var.data(vars_in{v})));
            [s1 s2 s3]=size(dtmp);
            dtmp_uconv=nan(s2,s3);

            for m=1:12
                % Unit conversion
                dtmp_uconv=squeeze(eval(['dtmp(m,:,:)' vars_conv{v}]));
                
                DATA_ERA5_05rs_out(:,:, m, yi, v)=squeeze(dtmp_uconv(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
            end

            clearvars dtmp dir_in fname;
        end
    end
end