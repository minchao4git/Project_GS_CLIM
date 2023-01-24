function [map_vegclm_r_rgb]=decomp_VegClm_relation_plot()

    global mon_vegclm_par_r mon_vegclm_par_p;
	global ax_show;
    global clm_cl_x;
    
    [s1 s2 nm nv]=size(mon_vegclm_par_r);

    map_vegclm_r_rgb=nan(s1,s2,nm,3);
    
    for m=1:nm
        for i=1:s1
            for j=1:s2

                wR=abs(mon_vegclm_par_r(i,j,m,1)); % temp.
                wG=abs(mon_vegclm_par_r(i,j,m,2)); % radiation
                wB=abs(mon_vegclm_par_r(i,j,m,3)); % soil moisture
                
                if ~isnan(wR)
                    map_vegclm_r_rgb(i,j,m,1:3)=WeightToRGB(wR, wG, wB);
                end

            end
        end
    end

    [lons,lats] = meshgrid(ax_show.south_s:0.5:ax_show.north_s, ax_show.west_s:0.5:ax_show.east_s);
    lons1=(lons(1:(end-1),1:(end-1))-0.25);
    lons2=(lons(1:(end-1),1:(end-1))+0.25);
    lats1=(lats(1:(end-1),1:(end-1))-0.25);
    lats2=(lats(1:(end-1),1:(end-1))+0.25);
if 1==2
    % ============== MAP ============== 
    figure('color','w','Position',[646   373   776   408]);
    gap=[0.06 0.06]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(1,2,gap,marg_h,marg_w);
    yscal=0.96; xscal=0.015;
    cb_xR=1.0;
    cb_yR=1.6;
    cb_wR=1.0;
	ax_show.frame=0;
    
	axes(ha(1));
    bg_show=squeeze(map_vegclm_r_rgb(:,:,13,:));
    mapH=axesm('MapProjection', 'miller', ...
      'MapLatLimit', [ax_show.south_s ax_show.north_s], ...
      'MapLonLimit', [ax_show.west_s ax_show.east_s]);
    axis off;
    framem;
    tightmap;
    
    hold on;
    % use only one big patchm
    lon_p=nan;lat_p=nan;cdata_ps=nan;
    n=0;
    for i=1:s1
        for j=1:s2
            cdata=squeeze(bg_show(i,j,:));
            if ~isnan(cdata)
                n=n+1;
                lon_ps(1,n)=lons1(i,j); 
                lon_ps(2,n)=lons2(i,j); 
                lon_ps(3,n)=lons2(i,j);
                lon_ps(4,n)=lons1(i,j);
                
                lat_ps(1,n)=lats2(i,j); 
                lat_ps(2,n)=lats2(i,j);
                lat_ps(3,n)=lats1(i,j);
                lat_ps(4,n)=lats1(i,j);
                
                cdata_ps(n,1:3)=cdata;
            end
        end
    end
	patchm(lon_ps,lat_ps,'FaceVertexCData',cdata_ps,'FaceColor', 'flat','EdgeColor','none');
    
    load coast;
    geoshow(mapH, flipud(lat), flipud(long));
    geoshow(mapH, flipud(lat), flipud(long), 'DisplayType', 'polygon', ...
        'FaceColor', 'white','LineStyle','none');
    plotm(flipud(lat), flipud(long),'Color','black');
    
    title(gca,'Early GS','FontSize',12);
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',12,'FontName','Arial','FontWeight','bold');
	set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');

	axes(ha(2));
    bg_show=squeeze(map_vegclm_r_rgb(:,:,14,:));
    mapH=axesm('MapProjection', 'miller', ...
      'MapLatLimit', [ax_show.south_s ax_show.north_s], ...
      'MapLonLimit', [ax_show.west_s ax_show.east_s]);
    axis off;
    framem;
	tightmap;
    
    hold on;
    % use only one big patchm
    lon_p=nan;lat_p=nan;cdata_ps=nan;
    n=0;
    for i=1:s1
        for j=1:s2
            cdata=squeeze(bg_show(i,j,:));
            if ~isnan(cdata)
                n=n+1;
                lon_ps(1,n)=lons1(i,j); 
                lon_ps(2,n)=lons2(i,j); 
                lon_ps(3,n)=lons2(i,j);
                lon_ps(4,n)=lons1(i,j);
                
                lat_ps(1,n)=lats2(i,j); 
                lat_ps(2,n)=lats2(i,j);
                lat_ps(3,n)=lats1(i,j);
                lat_ps(4,n)=lats1(i,j);
                
                cdata_ps(n,1:3)=cdata;
            end
        end
    end
	patchm(lon_ps,lat_ps,'FaceVertexCData',cdata_ps,'FaceColor', 'flat','EdgeColor','none');
    
    load coast;
    geoshow(mapH, flipud(lat), flipud(long));
    geoshow(mapH, flipud(lat), flipud(long), 'DisplayType', 'polygon', ...
        'FaceColor', 'white','LineStyle','none');
    plotm(flipud(lat), flipud(long),'Color','black');
    
    title(gca,'Late GS','FontSize',12);
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',12,'FontName','Arial','FontWeight','bold');
	set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
    
    export_fig 'decomp_VegClm_relation_map.png' -png -r600;
end

    % ============== Partial correlation by month ==============
    % ---- Plot by land class and climate zone ----
    global lc_dom_grp clmzone_grp;
    %                    nmon, ndriver, nlc, nclm
	mon_r_clmlc_summary=nan(12,3,8,3);
	legend_plot=false;
    xtmp=clm_cl_x;
    labels={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};

    lc_names={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
    clmz_name={'Boreal','Temperate', 'Medit.'};
    
    figure('Position',[144 222 1471 577],'Color','w');
    yscal=0.9; xscal=0.05;
    
    t = tiledlayout(3,6,'TileSpacing','Compact');
    n=0;
    for clm_grp=1:3
        col=0;
        
        for lc=xtmp(clm_grp,:)
            
            if clm_grp==3 && isnan(lc) && legend_plot==false
                nexttile;
                hold on;
                plot(1, 1,'-r','LineWidth',2); % temp.
                plot(1, 1,'-g','LineWidth',2); % rad.
                plot(1, 1,'-b','LineWidth',2); % sm.
                legend('Temp.','Rad.','Soil Moist.');
                set(gca,'XColor','none','YColor','none','TickDir','out');
                legend boxoff;
                set(gca,'Fontsize',12);
                hold off;
                legend_plot=true;
            end
            
            if isnan(lc)
                continue;
            end
            
            col=col+1;
            n=n+1;
            
            % spatial points
            amask=clmzone_grp; amask(amask~=clm_grp)=nan; amask(amask==clm_grp)=1;
            bmask=lc_dom_grp; bmask(bmask~=lc)=nan; bmask(bmask==lc)=1;
            p_thr=1;
            cmask=mon_vegclm_par_p; cmask(cmask>=p_thr)=nan;cmask(cmask<p_thr)=1;
            dtmp=mon_vegclm_par_r.*cmask.*repmat(amask,[1 1 nm 3]).*repmat(bmask,[1 1 nm 3]);
            
            for m=1:nm
                for v=1:3
                    if sum(sum(~isnan(dtmp(:,:,m,v)),1),2)<5 % too few number point
                        dtmp(:,:,m,v)=nan;
                    end
                end
            end

            mon_r_clmlc_summary(1:12,1:3,lc,clm_grp)=squeeze(nanmean(nanmean(dtmp(:,:,1:12,1:3),1),2));
            
            nexttile;
            hold on;
            plot(1:12, squeeze(nanmean(nanmean(dtmp(:,:,1:12,1),1),2)),'-r','LineWidth',2); % temp.
            plot(1:12, squeeze(nanmean(nanmean(dtmp(:,:,1:12,2),1),2)),'-g','LineWidth',2); % rad.
            plot(1:12, squeeze(nanmean(nanmean(dtmp(:,:,1:12,3),1),2)),'-b','LineWidth',2); % sm.

            plot([-100 100],[0 0],'--k');

            hold off
            box on;

            xticks([1:1:11.5]);
            xticklabels({'','SOS','+1','+2','+3','+4','+5','+6','+7','+8','+9'});
            set(gca,'Fontsize',8);
            title(lc_names{lc});

            if clm_grp==3
                xlabel('GS months');
            end
            
            if col==1
                ylabel(sprintf('%s\n%s',clmz_name{clm_grp},'Corr. changes (-)'),'Fontsize',11);
            end
            xlim([1.5 11.5]);
            ylim([-0.4 0.7]);
            
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{n},'FontSize',11,'FontName','Arial');
        end
    end
    export_fig 'decomp_VegClm_relation_plot_by_ClimLC.png' -png -r400;

if 1==1
    % ============== Partial correlation by climate zone ==============
    figure('Position',[1068         357         312         407],'Color','w');
    yscal=0.9; xscal=0.05;
    drv_color={'red','green','blue'};
	labels={'c)','d)','e)'};
    
%     t1 = tiledlayout(3,1,'TileSpacing','Compact');
    gap=[0.03 0.03]; marg_h=[0.08 0.05]; marg_w=[0.2 0.05];
    ha = tight_subplot(3,1,gap,marg_h,marg_w);
    
	for clm_grp=1:3
        axes(ha(clm_grp));
        pt_pos_low=nan(4,2);
        pt_pos_up=nan(4,2);

        for drv=1:3

            hold on;
            drv_mx=max(squeeze(mon_r_clmlc_summary(1:12,drv,:,clm_grp)), [], 2);
            drv_mean=nanmean(squeeze(mon_r_clmlc_summary(1:12,drv,:,clm_grp)), 2);
            drv_mn=min(squeeze(mon_r_clmlc_summary(1:12,drv,:,clm_grp)), [], 2);

            % align the length of the growing season
            for m=1:size(drv_mn,1)
                if (isnan(drv_mx(m)) || isnan(drv_mean(m)) || isnan(drv_mn(m)))
                    drv_mx(m)=nan;
                    drv_mean(m)=nan;
                    drv_mn(m)=nan;
                end
            end

            % get lower and upper bound of the polygon
            n=0;
            for m=1:size(drv_mn,1)
                if ~isnan(drv_mn(m))
                    n=n+1;
                    pt_pos_low(n,1)=m;
                    pt_pos_low(n,2)=drv_mn(m);
                end
            end

            n=0;
            for m=1:size(drv_mx,1)
                if ~isnan(drv_mx(m))
                    n=n+1;
                    pt_pos_up(n,1)=m;
                    pt_pos_up(n,2)=drv_mx(m);
                end
            end

            v1 = [pt_pos_low; flipud(pt_pos_up)];
            f1 = 1:size(pt_pos_low,1)*2;
            patch('Faces',f1,'Vertices',v1,'FaceColor',drv_color{drv},'FaceAlpha',.3, 'EdgeColor','none');
            plot((1:size(pt_pos_low,1))+1, drv_mean(~isnan(drv_mean)),drv_color{drv},'LineWidth',2);

            box on;
            xticks([1:1:11.5]);
            set(gca,'XTickLabel',{},'Box','on');
            
            set(gca,'Fontsize',8);
            xlim([1.5 11.5]);
            ylim([-0.2 0.6]);
            ylabel(sprintf('%s\n%s',clmz_name{clm_grp},'Corr. changes (-)'),'Fontsize',12);
            
            if clm_grp==3
                xticklabels({'','SOS','+1','+2','+3','+4','+5','+6','+7','+8','+9'});
                xlabel('GS months','FontSize',10);
            end
        end % drivers
        plot([0 100], [0 0],'--k');
        
        % Labels for sub-plots
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{clm_grp},'FontSize',14,'FontName','Arial');
    end
    export_fig 'decomp_VegClm_relation_plot_by_Clim.png' -png -r600;


    % ============== Partial correlation by drivers ==============
	figure('Position',[784   204   470   675],'Color','w');
    yscal=0.9; xscal=0.05;
    drv_color={'red','green','blue'};
    gap=[0.01 0.01]; marg_h=[0.08 0.05]; marg_w=[0.1 0.1];
    ha = tight_subplot(3,2,gap,marg_h,marg_w);
    
    drv_names={'Temp.','Rad.','Soil Moist.'};
    labels={'a)','b)','c)','d)','e)','f)'};
    
    n=0;
    for drv=1:3

        % ------ First columne ------
        n=n+1;
        axes(ha(n));
        
        bg_show=squeeze(mon_vegclm_par_r(:,:,13,drv))';
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
        if drv==1
            title('Early GS');
        end
        crng=[-0.7 0.7];
        colormap(gca,b2r(crng(1),crng(2),22));
        caxis([crng(1) crng(2)]);
        
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.4,drv_names{drv},'FontSize',11,'FontName','Arial','Rotation',90);
        
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{n},'FontSize',11,'FontName','Arial');

        % ------ First columne ------
        n=n+1;
        axes(ha(n));
        
        bg_show=squeeze(mon_vegclm_par_r(:,:,14,drv))';
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
        if drv==1
            title('Late GS');
        end
        
        crng=[-0.7 0.7];
        colormap(gca,b2r(crng(1),crng(2),22));
        caxis([crng(1) crng(2)]);
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{n},'FontSize',11,'FontName','Arial');
        
    end
    
    cbH=colorbar('SouthOutside');
    cb_x = cbH.Position(1);cb_y = cbH.Position(2);
    cb_w = cbH.Position(3);cb_h = cbH.Position(4);
    set(cbH, 'AxisLocationMode','manual');
    set(cbH, 'Position',[0.18 cb_y/2.5 cb_w*2.8 cb_h*0.8]);
end % 1==2
end