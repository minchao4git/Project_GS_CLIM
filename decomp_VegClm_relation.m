function decomp_VegClm_relation()

    global DATA_VI_05rs_dtds_out_gs1 DATA_CLM_05rs_dtds_out_gs1;
    global lgs_map;
    global mon_vegclm_par_r mon_vegclm_par_p;
    
    [s1 s2 s3 s4 s5]=size(DATA_CLM_05rs_dtds_out_gs1);
   
    dtmp_v=DATA_VI_05rs_dtds_out_gs1;  % Growing season months
    dtmp_c=DATA_CLM_05rs_dtds_out_gs1;

    % Definition of seasons
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[],[]};
    nm=size(m_rng,2); % number of months/seasons
    
    nv=s5;
    mon_vegclm_par_r=nan(s1,s2,nm,nv);
    mon_vegclm_par_p=nan(s1,s2,nm,nv);
    
    dsid=2; % Based on EVI2 growing season
    v_ctl_i={[2 3],[1 3], [1 2]}; % controlling variables index
    
    for v=1:nv
        for m=1:nm
            fprintf(sprintf('Climate variable: %d, Calculating month %d ...\n',v,m));
            for i=1:s1
               for j=1:s2
                   
                   m_lgs=lgs_map(i,j,dsid);
                   [m_rng{13} m_rng{14}]=get_ssn12(m_lgs);
                   m_rng{15}=[m_rng{13} m_rng{14}];
                   
                   % Partial correlation for the current month (no lag)
                   mon_v_ts=squeeze(nanmean(dtmp_v(i,j,m_rng{m},:,dsid),3));
                   mon_c_ts=squeeze(nanmean(dtmp_c(i,j,m_rng{m},:,:),3)); % position 5: climate variabiles

                   if nansum(abs(mon_v_ts))>0 % skip the ocean

                      x=mon_c_ts(:,v); % 1 variable
                      y=mon_v_ts;      % 1 variable
                      z=mon_c_ts(:,v_ctl_i{v}); % 2 controlling variables
                      [r, p] = partialcorr(x,y,z);
                      
                      mon_vegclm_par_r(i,j,m,v)=r(1,1);
                      mon_vegclm_par_p(i,j,m,v)=p(1,1);
                      
                   end
               end
            end
        end
    end
    
end