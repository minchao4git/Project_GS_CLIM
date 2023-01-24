function iav_by_landclass()
    global DATA_LC_out DATA_KG_CLMZ_05rs_out;    
    global lc_dom_grp clmzone_grp;
    global ax_show;
    %% Climate zone and landclass grouping
    % ====> IGBP Land cover
    % > The calculation is based on the regridded dataset from MODIS IGBP land
    % cover type
    % > Two appraoches are used to achieve the dominant landcover type, and
    % both yield similar results
    % > The landcover type based on 0.5 grid and 0.05 grid are quite similar,
    % and also consistent with landclass from IGBP: 
    % https://e4ftl01.cr.usgs.gov/MOTA/MCD12C1.006/2008.01.01/BROWSE.MCD12C1.A2008001.006.2018053184623.1.jpg

    % ----------------------------
    %  Eight groups of land class 
    % ----------------------------
    % Skip igbp 1: water boday

    % Method I: Freqency approach: get the stable land cover
    [s1 s2 s3 s4]=size(DATA_LC_out);
    nyr_lc=s3;
    for i=1:s1
        for j=1:s2
            for y=1:nyr_lc
                [dump lc_dom_ts(i,j,y)]=max(DATA_LC_out(i,j,y,:));
            end
        end
    end

    lc_dom=mode(lc_dom_ts,3); % get the most freqent value
    lc_dom_stable=nan(s1,s2);
    for i=1:s1
        for j=1:s2

            if sum((find(lc_dom_ts(i,j,:)==lc_dom(i,j))>0))==nyr_lc
                lc_dom_stable(i,j)=lc_dom(i,j);
            end
        end
    end

    % grouping
    lc_dom_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch lc_dom(i,j)
                % evergreen
                case {2 3}
                    lc_grp=1;
                % decidous
                case {4 5}
                    lc_grp=2;
                % mixed forest
                case {6}
                    lc_grp=3;
                % closed shrubland, open shrubland and woody savannas
                case {7,8,9}
                    lc_grp=4;
                % grasslands, savannas
                case {10, 11}
                    lc_grp=5;
                % permanent wetland
                case {12}
                    lc_grp=6;
                % crop land, natural vegetation mosaic (we skip 14: Urban here)
                case {13,15}
                    lc_grp=7;
                % snow and ice, barren sparsely veg.
                case {16,17}
                    lc_grp=8;
                otherwise
                    lc_grp=nan;
            end

            lc_dom_grp(i,j)=lc_grp;
        end
    end

    color_lc_grp=[[0 129 27]; 
                  [124 255 100]; 
                  [142 204 51]; 
                  [214 236 163]; 
                  [244 181 120]; 
                  [0 104 150]; 
                  [255 236 131]; 
                  [170 255 255]
    ]/255;
    %%
    % -----------------------------------
    %     Three groups of climate zones 
    % -----------------------------------
    % grouping
    clmzone_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch DATA_KG_CLMZ_05rs_out(i,j)
                % Dfb, Dfc, ET (snow, fully humid, warm and cool summer, polar tundra)
                case {42, 43, 62}
                    cl_grp=1;
                % Cfa, Cfb, Cfc (warm temperate, fully humid, hot-warm-cool summer)
                case {31, 32, 33}
                    cl_grp=2;
                % Csa, Csb, Csc (warm temperate, summer dry, hot-warm-cool summer)
                % BSk (arid steppe cold arid)
                case {34,   35,  36, 26}
                    cl_grp=3;
                otherwise
                    cl_grp=nan;
            end

            clmzone_grp(i,j)=cl_grp;
        end
    end
    % You can skip the group with number of grid point < 10
    color_clm_grp=[[226 68 211]
                  [0 121 24];
                  [255 254 118];
    ]/255;
    %%
    figure('color','white','Position',[555   451   613   422]);
    gap_h=0.005; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(1,2,gap,marg_h,marg_w);
    yscal=0.92; xscal=0.03;
    cb_xR=1.0;
    cb_yR=1.6;
    cb_wR=1.0;
	ax_show.frame=0;

    % --- Land class ---
	axes(ha(1));
    geoplot(ax_show, lc_dom_grp');
    
    colormap(gca,color_lc_grp);
    caxis([0.5 8.5]);
    cbh1=colorbar('southoutside');
    cbh1.Ticks = [1:8];

    cbh1.TickLabels={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
    title('Landclass','FontSize',11);
    resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, '',0.32,30,9);

    % --- Climate zone ---
	axes(ha(2));
    bg_show=clmzone_grp;
    bg_show(isnan(bg_show))=3; % nan in North Africa, set it as Medit.
    
    geoplot(ax_show, bg_show');
    
    colormap(gca,color_clm_grp);
    caxis([0.5 3.5]);
    cbh1=colorbar('southoutside');
    cbh1.Ticks = [1:3];
    cbh1.TickLabels={'Boreal','Temperate', 'Medit.'};
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');
    title('Climate zone','FontSize',11);
    resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, '',0.32,30,9);
end
