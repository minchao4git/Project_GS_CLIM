function []=var_to_clminx_month_calc_r()
    % ------------------------------------------------------------------
    % Correlation coefficient between variable of interests and teleconnection climate indices
    % ------------------------------------------------------------------

    global DATA_VI_05rs_dtds_out_gs1 DATA_CLMINX_out_gs1 DATA_CLMINX_out_gs;
    global v_clminx_r v_clminx_p;
    global lgs_map;
    
    dtmp=DATA_VI_05rs_dtds_out_gs1; % Growing season months

    [s1 s2 s3 nyr s5]=size(dtmp);
    % Definition of seasons (winter is actually starting from November, 1st month of the growing season)
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[] [] []};
    nm=size(m_rng,2); % number of months/seasons
    nv=10;
    nlag=6;
    v_clminx_r=nan(s1,s2,nm,nlag,nv);
    v_clminx_p=nan(s1,s2,nm,nlag,nv);
    
	dsid=2; % growing season definition using EVI2
    for v=1:nv
        
        fprintf(sprintf('=== > Variable : %d \n',v));
        
        for lag=0:(nlag-1)
            for m=1:nm

                fprintf(sprintf(' * Calculating lag %d for month %d ...\n',lag, m));

                for i=1:s1
                   for j=1:s2

                       m_lgs=lgs_map(i,j,dsid);
                       [m_rng{13} m_rng{14}]=get_ssn12(m_lgs);
                       m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months
                       m_rng{16}=[m_rng{13} m_rng{14}];  % all growing season months

                       if m<15
                           v_mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,dsid),3));

                           clminx_mon_ts=squeeze(nanmean(DATA_CLMINX_out_gs1(i,j,m_rng{m},:,lag+1,v),3));
                           
                       elseif m==15 % Entire GS
                           endm_i=length(m_rng{m});
                           v_mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m}(1:endm_i),:,dsid),3));
                           

                           sos=m_rng{13}(1);
                           clminx_mon_ts=squeeze(nanmean(DATA_CLMINX_out_gs(i,j,sos:sos,:,lag+1,v),3));
                           
                       elseif m==16 % Parts of GS
                           endm_i=length(m_rng{m});
                           
                           % First several months of growing season
                           if endm_i>=3
                               endm_i=3;
                           end

                           v_mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m}(1:endm_i),:,dsid),3));
                           
                           sos=m_rng{13}(1);
                           clminx_mon_ts=squeeze(nanmean(DATA_CLMINX_out_gs(i,j,sos:sos,:,lag+1,v),3));
                           
                       end

                       if (nansum(abs(v_mon_ts))>0) && sum(isnan(clminx_mon_ts))<nyr/2 % skip the ocean and nan values

                          [r, p]=corrcoef(v_mon_ts,clminx_mon_ts);
                          v_clminx_r(i,j,m,lag+1,v)=r(1,2);
                          v_clminx_p(i,j,m,lag+1,v)=p(1,2);

                       end
                   end % j
                end % i
            end % month
        end % lag
    end % v
    
end
