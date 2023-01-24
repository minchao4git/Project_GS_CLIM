function [] = preproc_ClmInx(use_gs)
    % Preprocessing teleconnection climate indices
    global data_src;
    global chosen_prd nyr;
    filename=sprintf('%s/Climate_indices/Teleconnection_indices_1950-2020.csv', data_src);

    DataStartLine = 2;
    NumVariables = 13;
    VariableNames  = {'yyyy','mm','NAO','EA','WP','EP','PNA','EA','SCA','TNH','POL','PT','ExplVar'};
    VariableWidths = [   4  ,  4 ,   6 ,  6 ,  6 ,  6 ,   6 ,  6 ,   6 ,   6 ,   6 ,   6,   10];
    DataType       = {'int16','int16','double','double','double','double','double','double','double','double','double','double','double'};
    opts = fixedWidthImportOptions('NumVariables',NumVariables,...
                                   'DataLines',DataStartLine,...
                                   'VariableNames',VariableNames,...
                                   'VariableWidths',VariableWidths,...
                                   'VariableTypes',DataType);

    clim_indx_tab = readtable(filename,opts);

    global DATA_VI_05rs_out;
    global DATA_CLMINX_out DATA_CLMINX_out_gs DATA_CLMINX_out_gs1;
    global nyr nyr_gs;
    global sos_map lgs_map;
    
    % Climate indices data with lag
    %                   n lags x n year x 12 variable (year, month, indx1, indx2,
    %                   indx3, indx4, indx5, indx6, indx7, indx8, indx9, indx10)
    nlag=6;
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    data_sy=1950;
    
	[s1 s2 s3 s4 s5]=size(DATA_VI_05rs_out);
    
    % ---- Get growing season data based on the defined growing season
	%                   lon/lat,lon/lat,month*year, nvar
    ni=10;
    DATA_CLMINX_out=nan(s1,s2,s3*s4,ni);
    DATA_CLMINX_out_gs=nan(s1,s2,s3,nyr_gs,nlag,ni);
    
    dtmp=table2array(clim_indx_tab(((yr_s-data_sy)*12+1):((yr_e-data_sy+1)*12),3:12));
    DATA_CLMINX_out(:,:,1:s3*s4,1:ni)=permute(repmat(dtmp,[1 1 s1 s2]),[3 4 1 2]);
    
    if use_gs

        % Growing season definition based on calculated SOS
        dsid=2; % growing season definition using EVI2
        for v=1:ni
            
            fprintf(sprintf('--> Processing climate index : %d\n',v));
            
            for i=1:s1
                for j=1:s2

                    gssn_s=sos_map(i,j,dsid);
                    sosmi=2; % sos month index in array
                    
                    for lag=0:(nlag-1)
                        gssn_s1=gssn_s-(sosmi-1)-lag;
                        
                        if gssn_s1<1
                            gssn_s1=gssn_s1+12;
                        end

                        if ~isnan(gssn_s)

                            for y=1:nyr_gs
                                % start from 1 month before the sos
                                m_s=(y-1)*12+gssn_s1;
                                m_e=m_s+11;

                                % adjusted 12-month growing season datasets on the shift of the SOS
                                % Different from DATA_CLMINX_out_gs, here all 12 months have values
                                DATA_CLMINX_out_gs(i,j,1:12,y,lag+1,v)=squeeze(DATA_CLMINX_out(i,j,m_s:m_e,v));
                            end % year
                        end
                    end % lag
                    
                end % lon/lat
            end % lon/lat
        end % variables
        
        % season 1 growing season dataset
        DATA_CLMINX_out_gs1=nan(s1,s2,12,nyr_gs,nlag, ni);
        for v=1:ni
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
                        DATA_CLMINX_out_gs1(i,j,gs_prd,:,:,v)=DATA_CLMINX_out_gs(i,j,gs_prd,:,:,v); 
                    end
                end
            end
        end
        
        ny=nyr_gs;
        dtmp=DATA_CLMINX_out_gs1;
    else
        ny=nyr;
        dtmp=DATA_CLMINX_out;
    end
   
end

