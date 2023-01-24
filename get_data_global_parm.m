function [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, dtype_str, pp_type)

global data_src;

    if any(strcmp(dtype_str,'SPOT_LAI'))
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/SPOT-PROBA/UHAM-ICDC/LAI/2000',data_src);
            file_name='c_gls_LAI__RT6_global_V2.0.1__0.5deg__20000110__UHAM-ICDC__v01.0.nc';

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end

    if any(strcmp(dtype_str,'GLEAM_SM'))
        if strcmp(pp_type,'RS_ORG')
            dir_in=sprintf('%s/GLEAM/3.3a/monthly/',data_src);
            file_name='SMsurf_1980_2018_GLEAM_v3.3a_MO.nc';
        end
        
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/GLEAM/3.3a/monthly/0.5deg',data_src);
            file_name='SMroot_1980_2018_GLEAM_v3.3a_MO_0.5deg_remapcon.nc';
        end

        % NetCDF format
        nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

        fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

        % ---- LON and LAT ---
        lons_tmp=nc_var.data('lon')';
        lats_tmp=nc_var.data('lat')';

        lon_n = length(lons_tmp);
        lat_n = length(lats_tmp);
        lons=repmat(lons_tmp', [1, lat_n]);
        lats=repmat(lats_tmp,  [lon_n, 1]);

        % Mask data for data this data format
        fla_mask=0;
        [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
        dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
        dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
    end

    if any(strcmp(dtype_str,'MODIS_LC'))
        if strcmp(pp_type,'RS_ORG')
            dir_in=sprintf('%s/MODIS/UHAM-ICDC/LandCoverType/',data_src);
            file_name='MODIS-C006_MCD12C1_landcover__LPDAAC__0.05deg__2002_fv0.02.nc';
        end
        
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/MODIS/UHAM-ICDC/LandCoverType/0.5deg',data_src);
            file_name='MODIS-C006_MCD12C1_landcover__LPDAAC__0.5deg__2001_fv0.02_remapcon.nc';
        end
        
        % NetCDF format (Resolution?)
        nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

        fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

        % ---- LON and LAT ---
        lons_tmp=nc_var.data('lon')';
        lats_tmp=nc_var.data('lat')';

        lon_n = length(lons_tmp);
        lat_n = length(lats_tmp);
        lons=repmat(lons_tmp', [1, lat_n]);
        lats=repmat(lats_tmp,  [lon_n, 1]);

        % Mask data for data this data format
        fla_mask=0;
        [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
        dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
        dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
    end

    if any(strcmp(dtype_str,'NDVI3g'))
        if strcmp(pp_type,'RS_ORG')
            file_name='ndvi3g_geo_v1_1982_0106.nc4';
            dir_in=sprintf('%s/GIMMS/NDVI3g/',data_src);

            % NetCDF format ("1/12 x 1/12 degrees", lonxlat: 4320x2160)
            nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
        
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/GIMMS/NDVI3g/0.5deg/',data_src);
            file_name='ndvi3g_geo_v1_1981_0712_0.5deg_remapbil.nc';

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end

    if any(strcmp(dtype_str,'VIP'))
        if strcmp(pp_type,'RS_ORG')
            dir_in=sprintf('%s/VIP/VIP15.004/1981.01.01',data_src);
            file_name='VIP15.A1981001.004.2016170095151.hdf';

            % NetCDF format ("1/12 x 1/12 degrees", lonxlat: 4320x2160)
            % HDF file does not need to get the sub-region manually with dm_mask.
            % the package matlab.io.hdfeos.* can automatically fulfill this.

            import matlab.io.hdfeos.*
            gfid = gd.open(sprintf('%s/%s',dir_in, file_name),'read');
            gridID = gd.attach(gfid,'VIP_CMG_GRID');

            cornerlat = [domain_def(3) domain_def(4)];
            cornerlon = [domain_def(1) domain_def(2)];
            regionID = gd.defBoxRegion(gridID,cornerlat,cornerlon);

            % Field name can be found from hdfinfo (e.g. finfo.Vgroup.Vgroup(1).SDS(2).Name)
            dtmp = double(gd.extractRegion(gridID,regionID, 'CMG 0.05 Deg 15DAYS EVI2'));
            dmn_lon_n_g=size(dtmp,1);
            dmn_lat_n_g=size(dtmp,2);

            clearvars dtmp;
        end
        
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/VIP/VIP15.004/netcdf/0.5deg',data_src);
            file_name='VIP15.EVI2_NDVI_1981.01.01_0.5deg_remapbil.nc';

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end
    
    if any(strcmp(dtype_str,'CRU'))
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/CRU/CRU-TS404/mon/orig/frs',data_src);
            file_name=sprintf('%s/cru_ts4.04.1901.1910.frs.dat.nc',dir_in);

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(file_name);

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end
    
    if any(strcmp(dtype_str,'ERA5'))
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/ECMWF/ERA5-Land/remap/0.5deg/t2m',data_src);
            file_name=sprintf('%s/t2m_ECMWF-ERA5-Land_rean_mon_1987-1987_0.5deg_remapycon.nc',dir_in);

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(file_name);

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end
    
    if any(strcmp(dtype_str,'KG_CLMZ'))
        if strcmp(pp_type,'RS_05')
            dir_in=sprintf('%s/Koppen-Geiger_Climate_Classification/obs/netcdf/',data_src);
            file_name=sprintf('%s/1976-2000.nc',dir_in);

            % NetCDF format ("0.5 degrees", lonxlat: 720x360)
            nc_var = ncgeodataset(file_name);

            fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

            % ---- LON and LAT ---
            lons=nc_var.data('lon')';
            lats=nc_var.data('lat')';

            % Mask data for data this data format
            fla_mask=0;
            [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
            dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
            dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
        end
    end
    
    if any(strcmp(dtype_str,'FLUXCOM'))

        dir_in=sprintf('%s/FluxCom/CarbonFluxes/RS_METEO/member/CRUJRA_v1/',data_src);
        file_name='GPP.RS_METEO.FP-DT.MLM-RF.METEO-CRUJRA_v1.720_360.monthly.1982.nc';

        % NetCDF format ("0.5 degrees", lonxlat: 720x360)
        nc_var = ncgeodataset(sprintf('%s/%s',dir_in,file_name));

        fprintf(sprintf('Data mask --> Processing file : %s\n',file_name));

        % ---- LON and LAT ---
        lons_tmp=nc_var.data('lon')';
        lats_tmp=nc_var.data('lat')';

        lon_n = length(lons_tmp);
        lat_n = length(lats_tmp);
        lons=repmat(lons_tmp', [1, lat_n]);
        lats=repmat(lats_tmp,  [lon_n, 1]);

        % Mask data for data this data format
        fla_mask=0;
        [x_s_e_g, y_s_e_g, ind_dmn_mask_g] = dm_mask(domain_def, lons, lats, fla_mask);
        dmn_lon_n_g=x_s_e_g(2)-x_s_e_g(1)+1;
        dmn_lat_n_g=y_s_e_g(2)-y_s_e_g(1)+1;
    end
end



