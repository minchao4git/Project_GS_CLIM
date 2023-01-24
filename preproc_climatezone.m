function [] = preproc_climatezone()
    
    % Local climate data from ERA5
    global data_src;
    global domain_def;
    global DATA_KG_CLMZ_05rs_out;
    
    % NOTE: the lon and lat order are different from other nc file
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'KG_CLMZ', 'RS_05');
    
    % check if the array dataset has created
    if size(DATA_KG_CLMZ_05rs_out,2) <= 1
        % Create 
        %                               nlon         nlat
        DATA_KG_CLMZ_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g);
    end

    % ---- Read data ---
    dir_in=sprintf('%s/Koppen-Geiger_Climate_Classification/obs/netcdf/',data_src);
    fname='1976-2000.nc';

    fprintf(sprintf('--> Processing file : %s\n',fname));

    fullname=sprintf('%s/%s',dir_in, fname);
    nc_var = ncgeodataset(fullname);

    % ---- Vegetation indices  ---
    % vars (lat, lon)
    dtmp=squeeze(double(nc_var.data('clsid')))';
    dtmp(dtmp<0)=nan;
    DATA_KG_CLMZ_05rs_out=dtmp(x_s_e_g(1):x_s_e_g(2),y_s_e_g(1):y_s_e_g(2));
end
