function [] = preproc_calcLocalClim(temp_int,pwi,use_gs)
    
    global DATA_ERA5_05rs_out DATA_SMRoot_05rs_out;
    global DATA_CLM_05rs_out_gs DATA_CLM_05rs_out_gs1;
    global DATA_CLM_05rs_dtds_out DATA_CLM_05rs_dtds_out_gs1;
    global DATA_CLM_05rs_dtdspw_out;
    global nyr_gs;
    global sos_map lgs_map;
    
	[s1 s2 s3 s4 s5]=size(DATA_ERA5_05rs_out);
    nv=3; % number of variables
    
    % ---- Get growing season data based on the defined growing season
	%                   lon/lat,lon/lat,month*year, nvar
    data_ts=nan(s1,s2,s3*s4,nv);
    DATA_CLM_05rs_out_gs=nan(s1,s2,s3,nyr_gs,nv);
    
    % t2m, ssrd: from ERA5
    data_ts(:,:,1:s3*s4,1:2)=reshape(DATA_ERA5_05rs_out(:,:,:,:,1:2),[s1 s2 s3*s4 2]);

    % SMRoot: from GLEAM
    data_ts(:,:,1:s3*s4,3)=DATA_SMRoot_05rs_out(:,:,:,1);
    
    if use_gs

        dsid=2; % growing season definition using EVI2
        for v=1:nv
            for i=1:s1
                for j=1:s2
                    % Growing season definition based on astronomic month starting from December
                    % gssn_s=12;

                    gssn_s=sos_map(i,j,dsid);

                    sosmi=2; % sos month index in array
                    gssn_s1=gssn_s-(sosmi-1);
                    if gssn_s1<1
                        gssn_s1=gssn_s1+12;
                    end

                    if ~isnan(gssn_s)
%                       m_s=(y-1)*12+gssn_s;
                        for y=1:nyr_gs
                            % start from 1 month before the sos
                            m_s=(y-1)*12+gssn_s1;
                            m_e=m_s+11;

                            % adjusted 12-month growing season datasets on the shift of the SOS
                            DATA_CLM_05rs_out_gs(i,j,1:12,y,v)=squeeze(data_ts(i,j,m_s:m_e,v));
                        end % year
                    end 
                end % lon/lat
            end % lon/lat
        end % variables
        
        % season 1 growing season dataset
        DATA_CLM_05rs_out_gs1=nan(s1,s2,12,nyr_gs,3);
        for v=1:nv
            for i=1:s1
                for j=1:s2

                    m_lgs=lgs_map(i,j,dsid);
                    if ~isnan(m_lgs)

                        if m_lgs < 12
                            gs_prd=(1:m_lgs)+(sosmi-1);
                        else
                            gs_prd=1:12;
                        end

                        % mean of actual growing season (season 1)
                        DATA_CLM_05rs_out_gs1(i,j,gs_prd,:,v)=DATA_CLM_05rs_out_gs(i,j,gs_prd,:,v); 
                    end
                end
            end
        end
        
        ny=nyr_gs;
        dtmp=DATA_CLM_05rs_out_gs1;
    end
   
    % --- Detrend and deseasonalize growing season dataset ---
    if strcmp(temp_int,'monthly_dsclim')

        % lon/lat lon/lat yr ds
        dtmp_movm=movmean(squeeze(nanmean(dtmp,3)),5, 3); % 5-year running mean of annual mean

        % Detrend the monthly data with the detrended annual mean
        atmp=repmat(dtmp_movm,[1 1 1 1 12]);
        dtmp=dtmp-permute(atmp,[1 2 5 3 4]);
        
        % Deseasonalize with the climatology seasonal cycle
        for v = 1:nv

            climssn=squeeze(nanmean(dtmp(:,:,:,:,v),4));
            dtmp(:,:,:,:,v)=dtmp(:,:,:,:,v)-repmat(climssn,[1 1 1 ny]);
            
            % Prewithening along the interannual direction
            nts=0;
            if pwi==1
                nhi=2;    % Highest AR order to consider
                k2=[1 2]; % set k2(1):1 to fit model from 1 to the defined nhi, 
                          % k2(2): 2 accept the lowest AIC
                          
                for i=1:s1
                    for j=1:s2
                        for m=1:12
                            ts_tmp=squeeze(dtmp(i,j,m,:,v));
                            
                            if sum(isnan(ts_tmp))==0
                                nts=nts+1;
                                fprintf(sprintf('Prewhitening DS: %d, %d time series data ... \n',v, nts));
                                [DATA_CLM_05rs_dtdspw_out(i,j,m,:,v),dump1,dump2,dump3] = whit1(ts_tmp,nhi,k2);
                            end
                        end
                    end
                end
            end    

        end
        
        if use_gs
            DATA_CLM_05rs_dtds_out_gs1=dtmp;
        else
            DATA_CLM_05rs_dtds_out=dtmp;
        end
        clearvars dtmp climssn dtmp_movm atmp;
    end

end

