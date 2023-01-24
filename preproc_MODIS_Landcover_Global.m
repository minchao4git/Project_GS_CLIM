function [] = preproc_MODIS_Landcover_Global()
    % ---- Landcover ----
    global data_src;
    global domain_def;
    global x_s_e_g;
    global y_s_e_g;
    global dmn_lon_n_g;
    global dmn_lat_n_g;
    global chosen_prd;
    global DATA_LC_out;

    % Data source I: Europe domain (get: x_s_e, y_s_e, ind_dmn_mask)
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'MODIS_LC', 'RS_05');

    % period for complete annual data
    yr_s=2001;
    yr_e=chosen_prd(end);

    % check if the array dataset has created
    if size(DATA_LC_out,2) <= 1
        % Create 
        %                   nlon         nlat          nyear     nigbp land class 
        DATA_LC_out = nan(dmn_lon_n_g, dmn_lat_n_g, (yr_e-yr_s)+1, 17);
    end
    
    var_name={
        'water_igbp',                                % 1
        'evergreen_needleleaf_forest_igbp',          % 2
        'evergreen_broadleaf_forest_igbp',           % 3  
        'deciduous_needleleaf_forest_igbp',          % 4
        'deciduous_broadleaf_forest_igbp',           % 5
        'mixed_forest_igbp',                         % 6
        'closed_shrublands_igbp',                    % 7
        'open_shrublands_igbp',                      % 8
        'woody_savannas_igbp',                       % 9
        'savannas_igbp',                             % 10
        'grasslands_igbp',                           % 11
        'permanent_wetlands_igbp',                   % 12
        'croplands_igbp',                            % 13
        'urban_and_builtup_igbp',                    % 14
        'cropland_natural_vegetation_mosaic_igbp',   % 15
        'snowandice_igbp',                           % 16
        'barren_sparsely_vegetated_igbp'             % 17
    };

    for y=yr_s:yr_e

        % ---- Read data ---
        dir_in=sprintf('%s/MODIS/UHAM-ICDC/LandCoverType/0.5deg',data_src);
        fname=sprintf('MODIS-C006_MCD12C1_landcover__LPDAAC__0.5deg__%d_fv0.02_remapcon.nc',y);
        
        fprintf(sprintf('--> Processing file : %s\n',fname));

        fullname=sprintf('%s/%s',dir_in, fname);
        nc_var = ncgeodataset(fullname);

        % ---- IGBP land type  ---
        for v=1:size(var_name,1)
            dtmp=double(squeeze(nc_var.data(var_name{v}))');
            DATA_LC_out(:,:,(y-yr_s)+1,v)=dtmp(x_s_e_g(1):x_s_e_g(2),y_s_e_g(1):y_s_e_g(2));
        end
    end

end
