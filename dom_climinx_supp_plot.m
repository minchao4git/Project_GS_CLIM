function dom_climinx_supp_plot()

    global v_clminx_r v_clminx_p;
	global ax_show;
    global lons lats;
    
    [s1 s2 s3 s4 s5]=size(v_clminx_r);
    
    clearvars dtmp d_mxi d_mx;
    clmi_chosen=[1 2 3 5 6 7 9 ];
    clmi_names={'NAO','EA','WP','EP','PNA','EAWR','SCA','TNH','POL','PT'};
    nclmi=length(clmi_chosen);
    dtmp=v_clminx_r(:,:,:,:,clmi_chosen);
    
    % Significant points
    lons_sig=repmat(lons,[1 1 s3 s4 nclmi]);
    lats_sig=repmat(lats,[1 1 s3 s4 nclmi]);
    lons_sig(v_clminx_p(:,:,:,:,clmi_chosen)>0.05)=nan;
    lats_sig(v_clminx_p(:,:,:,:,clmi_chosen)>0.05)=nan;

    nmn=size(v_clminx_r,3);

    n=0;
	figure('color','w','Position',[415        -322         610        1282]);
    gap_h=0.005; gap_w=0.002;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(7,4,gap,marg_h,marg_w);
    title_s={'Entire GS', 'SOS','+2 mon.','+4 mon.'};
    ax_show.frame=0;
    
	global map_mxclm_lon_sig map_mxclm_lat_sig;
    map_mxclm_lon_sig=nan(s1,s2);
    map_mxclm_lat_sig=nan(s1,s2);
    
    for clmi=1:nclmi
        clearvars d_mxi d_mx map_mxclm_r;
        [d_mx d_mxi]=nanmax(abs(dtmp(:,:,:,:,clmi)),[],4); % get the highest values among lags
        d_mxi(isnan(d_mx))=0;

        d_mxi=squeeze(d_mxi);
        map_mxclm_r=nan(s1,s2,nmn);

        for i=1:s1
            for j=1:s2
                for m=1:nmn
                    mxi=d_mxi(i,j,m);
                    if mxi > 0
                        map_mxclm_r(i,j,m)=dtmp(i,j,m,mxi,clmi);
                        map_mxclm_lon_sig(i,j,m)=lons_sig(i,j,m,mxi,clmi); % to check if the largest values are significant, then get its coordinates
                        map_mxclm_lat_sig(i,j,m)=lats_sig(i,j,m,mxi,clmi);
                    end
                end
            end
        end

        nn=0;
        for m=[15 2 4 6]
            n=n+1;
            nn=nn+1;
            axes(ha(n));

            hold on;
            bg_show=squeeze(map_mxclm_r(:,:,m))';
            bg_show(isnan(bg_show))=0;
            geoplot(ax_show,bg_show);
            
            % plot the significant gridpoints
            a=squeeze(map_mxclm_lat_sig(1:2:end,1:2:end,m));
            b=squeeze(map_mxclm_lon_sig(1:2:end,1:2:end,m));
            plotm(b(:), a(:), 'Color',[80 80 80]/255,'LineStyle','none', 'Marker', '+', 'MarkerSize',2);
            
            colormap(gca, flipud(cbrewer('div','RdBu',21)));
            caxis([-0.8 0.8]);
            set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
            
            if nn==1
                a=get(gca);
                gca_w = (a.XLim(2)-a.XLim(1));
                gca_h = (a.YLim(2)-a.YLim(1));
                text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.35,sprintf('%s',clmi_names{clmi_chosen(clmi)}),'FontSize',9,'FontName','Arial','Rotation',90);
            end
            
            if clmi==1
                title(sprintf('%s',title_s{nn}));
            end
        end
    end % Climate index
    cbH = colorbar('southoutside');
    resizeCB(cbH, 0.31, 0.5, 7.5, 1.4, 'Max Corr.(-)',0.05,300,10);
end