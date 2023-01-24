function []=month_to_tot_month_anom_calc_cv(use_gs)
    % ------------------------------------------------------------------
    % The importance of monthly anomalies to the total monthly anomalies
    % ------------------------------------------------------------------

    global nds ds_s ds_e;
    global window_wd nw nw_gs;
    global DATA_VI_05rs_out_gs DATA_VI_05rs_out_gs1;
    global DATA_VI_05rs_dtds_out;
    global DATA_VI_05rs_dtds_out_gs1;
    global mon_cv mon_std;
    global mw_cv mw_std;
    global cv_trnd cv_trnd_p std_trnd std_trnd_p;
    global lgs_map;
    
    % interannual variability of monthly anomalies ratio
    if use_gs
        dtmp=DATA_VI_05rs_dtds_out_gs1; % Growing season months
        utmp=DATA_VI_05rs_out_gs1;
    else
        dtmp=DATA_VI_05rs_dtds_out; 
        utmp=DATA_VI_05rs_out_gs;
    end

    [s1 s2 s3 s4 s5]=size(dtmp);
    % Definition of seasons (winter is actually starting from November, 1st month of the growing season)
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[],[]};
    nm=size(m_rng,2); % number of months/seasons
    
    mon_cv=nan(s1,s2,nm,nds); mon_std=nan(s1,s2,nm,nds);
    for dsid=ds_s:ds_e
        
        for m=1:nm
            fprintf(sprintf('DS_ID: %d, Calculating month %d ...\n',dsid,m));
            for i=1:s1
               for j=1:s2
                   
                   m_lgs=lgs_map(i,j,dsid);
                   [m_rng{13} m_rng{14}]=get_ssn12(m_lgs); % 1st and 2nd halfs
                   m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months
                   
                   mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,dsid),3));
                   if nansum(abs(mon_ts))>0 % skip the ocean
                                          
                      gs_u=squeeze(nanmean(nanmean(utmp(i,j,m_rng{15},:,dsid),3),4));
                      mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,dsid),3));
                      
                      if abs(gs_u)>0.1
                        mon_cv(i,j,m,dsid)=std(mon_ts,'omitnan')./abs(gs_u);
                        mon_std(i,j,m,dsid)=std(mon_ts,'omitnan');
                        
                        if m>=nm
                            mon_cv(i,j,m+1,dsid)=nanmean(mon_cv(i,j,m_rng{13},dsid),3);
                            mon_cv(i,j,m+2,dsid)=nanmean(mon_cv(i,j,m_rng{14},dsid),3);
                            mon_std(i,j,m+1,dsid)=nanmean(mon_std(i,j,m_rng{13},dsid),3);
                            mon_std(i,j,m+2,dsid)=nanmean(mon_std(i,j,m_rng{14},dsid),3);
                        end
                      end
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
    mw_cv=nan(s1,s2,nm,n,nds); mw_std=nan(s1,s2,nm,n,nds);
    
    for dsid=ds_s:ds_e
        for mwi=1:n % number of moving windows
            window4corcf=[mwi:(mwi+window_wd-1)];
            fprintf(sprintf('DS_ID: %d, === Window #%d ===\n',dsid,mwi));

            for m=1:nm
                fprintf(sprintf('Calculating moving windows for month %d ...\n',m));

                for i=1:s1
                   for j=1:s2
                       
                       m_lgs=lgs_map(i,j,dsid);
                       [m_rng{13} m_rng{14}]=get_ssn12(m_lgs); % 1st and 2nd halfs
                       m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months

                       mon_ts_w=squeeze(nanmean(dtmp(i,j,m_rng{m},window4corcf,dsid),3));
                       if nansum(abs(mon_ts_w))>0 % skip the ocean

                          % CV
                          gs_u_w=squeeze(nanmean(nanmean(utmp(i,j,m_rng{15},window4corcf,dsid),3),4));
                           % exclude the extreme small value before calculate
                          if gs_u_w>0.1
                            mw_cv(i,j,m,mwi,dsid)=std(mon_ts_w,'omitnan')./gs_u_w;
                            mw_std(i,j,m,mwi,dsid)=std(mon_ts_w,'omitnan');
                          end
                       end
                   end
                end
            end
        end
    end

    % ------------------------------------------------------
    % Calculate the trend of CV (TREND OF 15 MOVING WINDOES)
    % ------------------------------------------------------
    cv_trnd=nan(s1,s2,nm,nds);
    cv_trnd_p=nan(s1,s2,nm,nds);
    nE=0;
    for dsid=ds_s:ds_e

        for m=1:nm
            fprintf(sprintf('DS_ID: %d, Calculating trend of CV for month %d ...\n',dsid,m));
            for i=1:s1
               for j=1:s2
                   a_tmp=squeeze(mw_cv(i,j,m,:,dsid));

                    if nansum(abs(a_tmp))>0 % skip the ocean
                        % Linear trend
                        try
                            
                            mdl=fitlm((1:size(mw_cv,4))',a_tmp);
                            cv_trnd(i,j,m,dsid)=mdl.Coefficients.Estimate(2);
                            cv_trnd_p(i,j,m,dsid)=mdl.coefTest;

                        catch
                            nE=nE+1;
                            fprintf(sprintf('DS_ID: %d, Error exists! Number of error: %d ...\n',dsid,nE));
                            cv_trnd(i,j,m,dsid)=-999;
                            cv_trnd_p(i,j,m,dsid)=-999;
                        end
                   end
               end
            end
        end
    end
    
    % -------------------------------------------------------
    % Calculate the trend of STD (TREND OF 15 MOVING WINDOES)
    % -------------------------------------------------------
    std_trnd=nan(s1,s2,nm,nds);
    std_trnd_p=nan(s1,s2,nm,nds);
    nE=0;
    for dsid=ds_s:ds_e
        for m=1:nm
            fprintf(sprintf('DS_ID: %d, Calculating trend of STD for month %d ...\n',dsid,m));
            for i=1:s1
               for j=1:s2
                   a_tmp=squeeze(mw_std(i,j,m,:,dsid));

                    if nansum(abs(a_tmp))>0 % skip the ocean
                        % Linear trend
                        try

                            mdl=fitlm((1:size(mw_std,4))',a_tmp);
                            std_trnd(i,j,m,dsid)=mdl.Coefficients.Estimate(2);
                            std_trnd_p(i,j,m,dsid)=mdl.coefTest;
                        catch
                            nE=nE+1;
                            fprintf(sprintf('DS_ID: %d, Error exists! Number of error: %d ...\n',dsid,nE));
                            std_trnd(i,j,m,dsid)=-999;
                            std_trnd_p(i,j,m,dsid)=-999;
                        end
                   end
               end
            end
        end
    end
end
