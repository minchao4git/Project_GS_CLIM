function [] = preproc_GLEAM()
    % ---- Landcover ----
    global data_src;
    global domain_def;
    global x_s_e_g;
    global y_s_e_g;
    global dmn_lon_n_g;
    global dmn_lat_n_g;
    global chosen_prd;
    global DATA_SMRoot_05rs_out;

    % Data source I: Europe domain (get: x_s_e, y_s_e, ind_dmn_mask)
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'GLEAM_SM', 'RS_05');

    % Available data period
    data_sy=1980; % data start year (from Jan.)
    
    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    nyr=yr_e-yr_s+1;

    % check if the array dataset has created
    if size(DATA_SMRoot_05rs_out,2) <= 1
        % Create 
        %                        nlon         nlat          nyear*nmonth     # of soil moisture var.
        DATA_SMRoot_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, ((yr_e-yr_s)+1)*12, 1);
    end
    

    % ---- Read data ---
    dir_in=sprintf('%s/GLEAM/3.3a/monthly/0.5deg',data_src);
    fname='SMroot_1980_2018_GLEAM_v3.3a_MO_0.5deg_remapcon.nc';

    fprintf(sprintf('--> Processing file : %s\n',fname));

    fullname=sprintf('%s/%s',dir_in, fname);
    nc_var = ncgeodataset(fullname);

    % Var. dimension in netcdf file
    % (time (month x year), lon, lat)
    dtmp=double(squeeze(nc_var.data('SMroot')));
    DATA_SMRoot_05rs_out=permute(dtmp(((yr_s-data_sy)*12+1):((yr_e-data_sy+1)*12), x_s_e_g(1):x_s_e_g(2), y_s_e_g(1):y_s_e_g(2)),[2 3 1]);
    
%     clearvars dtmp;
end
