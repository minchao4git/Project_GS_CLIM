function []=month_to_tot_month_anom_calc_ratio(use_gs)
    % ------------------------------------------------------------------
    % The importance of monthly anomalies to the total monthly anomalies
    % ------------------------------------------------------------------

    global nds ds_s ds_e mdi;
    global window_wd nw nw_gs;
    global DATA_VI_05rs_out DATA_VI_05rs_sm_out;
    global DATA_VI_05rs_out_gs DATA_VI_05rs_out_gs1;
    global DATA_VI_05rs_dtds_out;
    global DATA_VI_05rs_dtds_out_gs DATA_VI_05rs_dtds_out_gs1;
    global DATA_VI_05rs_dtdspw_out;
    global mon_anm_ratio mon_anm_ratio_mean;
    global mw_ratio;
    global ratio_trnd ratio_trnd_p;
    global lgs_map;
    global nyr_gs;
    global mask_gs;
    
	data_astr=squeeze(DATA_VI_05rs_sm_out(:,:,:,:,:,mdi));
	[d1 d2 d3 d4 d5]=size(data_astr);
    data_astr=data_astr.*repmat(squeeze(mask_gs(:,:,mdi)),[1 1 d3 d4 d5]);
    
    % interannual variability of monthly anomalies ratio
    if use_gs
        dtmp=DATA_VI_05rs_dtds_out_gs1; % Growing season months
    else
        dtmp=DATA_VI_05rs_dtds_out; % Growing season months
    end

    [s1 s2 s3 s4 s5]=size(dtmp);
    % Definition of seasons (winter is actually starting from November, 1st month of the growing season)
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[]};
    nm=size(m_rng,2); % number of months/seasons
    
    % annual value array
    %   lat/lon, lat/lon, y, ds
    % Based on calculated growing season
	if use_gs
        tmp=squeeze(nanmean(DATA_VI_05rs_out_gs1,3));
    else
        tmp=squeeze(nanmean(data_astr,3));
	end

    clearvars tmp;

    mon_anm_ratio=nan(s1,s2,nm,s4,nds);
    mon_anm_ratio_mean=nan(s1,s2,nm,s4,nds);
    for dsid=ds_s:ds_e
        for m=1:nm
            fprintf(sprintf('DS_ID: %d, Calculating month %d ...\n',dsid,m));
            for i=1:s1
               for j=1:s2
                   
                   m_lgs=lgs_map(i,j,dsid);
                   [m_rng{13} m_rng{14}]=get_ssn12(m_lgs);
                    
                   mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,dsid),3));

                   if nansum(abs(mon_ts))>0 % skip the ocean
                                          
                      % ABS anomalies ratio
                      gs_anm_ts=squeeze(nansum(abs(dtmp(i,j,[m_rng{13} m_rng{14}],:,dsid)),3));
                      m_anm_ts=squeeze(nansum(abs(dtmp(i,j,m_rng{m},:,dsid)),3));
                      m_anm_ts_mean=squeeze(nanmean(abs(dtmp(i,j,m_rng{m},:,dsid)),3));
                      mon_anm_ratio(i,j,m,:,dsid)=m_anm_ts./gs_anm_ts;
                      mon_anm_ratio_mean(i,j,m,:,dsid)=m_anm_ts_mean./gs_anm_ts;
                   end
               end
            end
        end
    end

    return;
    
    % ------------------------------------------------------
    %  Monthly control (20 YEARS MOVING WINDOW MEAN RATIO)
    % ------------------------------------------------------
    if use_gs
        n=nw_gs;
    else
        n=nw;
    end
    mw_ratio=nan(s1,s2,nm,n,nds);
    
    for dsid=ds_s:ds_e
        for mwi=1:n % number of moving windows
            window4corcf=[mwi:(mwi+window_wd-1)];
            fprintf(sprintf('DS_ID: %d, === Window #%d ===\n',dsid,mwi));

            for m=1:nm
                fprintf(sprintf('Calculating moving windows for month %d ...\n',m));

                for i=1:s1
                   for j=1:s2
                       
                       mw_ratio(i,j,m,mwi,dsid)=squeeze(nanmean(mon_anm_ratio(i,j,m,window4corcf,dsid),4));

                   end
                end
            end
        end
    end

    % -------------------------------------------------
    % Calculate the trend (TREND OF 15 MOVING WINDOES)
    % -------------------------------------------------
    ratio_trnd=nan(s1,s2,nm,nds);
    ratio_trnd_p=nan(s1,s2,nm,nds);
    nE=0;
    for dsid=ds_s:ds_e
        for m=1:nm
            fprintf(sprintf('DS_ID: %d, Calculating trend for month %d ...\n',dsid,m));
            for i=1:s1
               for j=1:s2
                   a_tmp=squeeze(mw_ratio(i,j,m,:,dsid));

                    if nansum(abs(a_tmp))>0 % skip the ocean
                        % Linear trend
                        try
                            
                             mdl=fitlm((1:size(mw_ratio,4))',a_tmp);

                            ratio_trnd(i,j,m,dsid)=mdl.Coefficients.Estimate(2);
                            ratio_trnd_p(i,j,m,dsid)=mdl.coefTest;
                        catch
                            nE=nE+1;
                            fprintf(sprintf('DS_ID: %d, Error exists! Number of error: %d ...\n',dsid,nE));
                            ratio_trnd(i,j,m,dsid)=-999;
                            ratio_trnd_p(i,j,m,dsid)=-999;
                        end
                   end
               end
            end
        end
    end
end
