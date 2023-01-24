function [] = preproc_calcVI_I(temp_int,pwi,use_gs)
    
    % Input
    global DATA_VI_05rs_sm_out 
    % Output
    global DATA_VI_05rs_out_gs DATA_VI_05rs_out_gs1;
    global DATA_VI_05rs_dtds_out DATA_VI_05rs_dtds_out_gs1;
    global DATA_VI_05rs_dtdspw_out;
    global nyr nyr_gs;
    global nds mdi;
    global sos_map lgs_map;
    global mask_gs;
    
    datain=squeeze(DATA_VI_05rs_sm_out(:,:,:,:,:,mdi));
	[s1 s2 s3 s4 s5]=size(datain);
    datain=datain.*repmat(squeeze(mask_gs(:,:,mdi)),[1 1 s3 s4 s5]);
    
    % ---- Get growing season data based on the defined growing season
    if use_gs
        DATA_VI_05rs_out_gs=nan(s1,s2,12,nyr_gs,nds);

        adtmp_rs=reshape(datain,[s1 s2 s3*s4 s5]);

        % Growing season definition based on calculated SOS
        for dsid=1:nds
            for i=1:s1
                for j=1:s2

                    gssn_s=sos_map(i,j,dsid);

                    sosmi=2; % sos month index in array
                    gssn_s1=gssn_s-(sosmi-1);
                    if gssn_s1<1
                        gssn_s1=gssn_s1+12;
                    end

                    if ~isnan(gssn_s)
                        for y=1:nyr_gs
                            % % TOREMOVE comments: start from 1 month before the sos
                            m_s=(y-1)*12+gssn_s1;
                            m_e=m_s+11;

                            DATA_VI_05rs_out_gs(i,j,1:12,y,dsid)=squeeze(adtmp_rs(i,j,m_s:m_e,dsid));
                        end % year
                    end 
                end % lon/lat
            end % lon/lat
        end % data source
        
        % season 1 growing season dataset
        DATA_VI_05rs_out_gs1=nan(s1,s2,12,nyr_gs,nds);
        for dsid=1:nds
            for i=1:s1
                for j=1:s2

                    % Fill value based on the growing season length
                    m_lgs=lgs_map(i,j,dsid);
                    if ~isnan(m_lgs)

                        if m_lgs < 12
                            gs_prd=(1:m_lgs)+(sosmi-1);
                        else
                            gs_prd=1:12;
                        end

                        % mean of actual growing season (season 1)
                        DATA_VI_05rs_out_gs1(i,j,gs_prd,:,dsid)=DATA_VI_05rs_out_gs(i,j,gs_prd,:,dsid);
                    end
                end
            end
        end
        
        % MASKING
        DATA_VI_05rs_out_gs(DATA_VI_05rs_out_gs<=0)=nan;
        DATA_VI_05rs_out_gs1(DATA_VI_05rs_out_gs1<=0)=nan;
        
        ny=nyr_gs;
        dtmp=DATA_VI_05rs_out_gs1;
    else
        ny=nyr;
        dtmp=datain;
    end
    
   
    % --- Detrend and deseasonalize growing season dataset ---
    if strcmp(temp_int,'monthly_dsclim')

        % lon/lat lon/lat yr ds
        dtmp_movm=movmean(squeeze(nanmean(dtmp,3)),5, 3); % 5-year running mean of annual mean

        % Detrend the monthly data with the detrended annual mean
        atmp=repmat(dtmp_movm,[1 1 1 1 12]);
        dtmp=dtmp-permute(atmp,[1 2 5 3 4]);
        
        % Deseasonalize with the climatology seasonal cycle
        for dsid = 1:nds
            climssn=squeeze(nanmean(dtmp(:,:,:,:,dsid),4));
            dtmp(:,:,:,:,dsid)=dtmp(:,:,:,:,dsid)-repmat(climssn,[1 1 1 ny]);
            
            % Prewithening along the interannual direction
            nts=0;
            if pwi==1
                nhi=2;    % Highest AR order to consider
                k2=[1 2]; % set k2(1):1 to fit model from 1 to the defined nhi, 
                          % k2(2): 2 accept the lowest AIC
                          
                for i=1:s1
                    for j=1:s2
                        for m=1:12
                            ts_tmp=squeeze(dtmp(i,j,m,:,dsid));
                            
                            if sum(isnan(ts_tmp))==0
                                nts=nts+1;
                                fprintf(sprintf('Prewhitening DS: %d, %d time series data ... \n',dsid, nts));
                                [DATA_VI_05rs_dtdspw_out(i,j,m,:,dsid),dump1,dump2,dump3] = whit1(ts_tmp,nhi,k2);
                            end
                        end
                    end
                end
            end    

        end
        
        if use_gs
            DATA_VI_05rs_dtds_out_gs1=dtmp;
        else
            DATA_VI_05rs_dtds_out=dtmp;
        end
        clearvars dtmp climssn dtmp_movm atmp;
    end
    
end

