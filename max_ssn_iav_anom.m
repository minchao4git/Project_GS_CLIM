function []=max_ssn_iav_anom()
    
    % Group it by season: mean of the whole period
    global mon_anm_ratio mon_anm_ratio_mean mon_cv;
    global nw_gs;
    global lgs_map;
    global ax_show;

    dss=1; dse=2;
    
    [s1 s2 s3 s4 s5]=size(mon_anm_ratio_mean);
    ssn_map=nan(s1,s2,4,3);
    dom_ssn_map=nan(s1,s2,3);
    dtmp1=mon_cv;
    dtmp1(isinf(dtmp1))=nan;
    for dsid=dss:dse
        for i=1:s1
            for j=1:s2
                 ts_tmp=squeeze(dtmp1(i,j,[13 14 15],dsid));
                if sum(ts_tmp)>0
                    dom_ssn_map(i,j,dsid)=((ts_tmp(1)-ts_tmp(2)))/ts_tmp(3);
                end
            end
        end
    end

    figure('color','white','Position',[432   462   995   391]);
    gap_h=0.005; gap_w=0.008;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(1,3,gap,marg_h,marg_w);
    yscal=0.92; xscal=0.03;
    cb_xR=1.0;
    cb_yR=1.1;
    cb_wR=1.0;
	ax_show.frame=1;
    ax_show.contourf=0;
        
    % ====== interannual variability ======
	axes(ha(1));
    stb_map=squeeze(nanmean(mon_cv(:,:,15,dss:dse),4))';
    geoplot(ax_show, stb_map);

    ctmp=colormap(parula(400));
    ctmp(1,:)=[1 1 1];
    colormap(gca,ctmp);
    cbh1=colorbar('southoutside');
    set(cbh1, 'AxisLocationMode','manual');
    caxis([0 0.08]);
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
    title('Interannual variability');
    resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, 'CV_G_S',0.32,30,9);

    % ====== Maximum seasonal variability by season ======
	axes(ha(2));
    set(gca,'Fontsize',8);

    use_gs=1;
    if use_gs
        n=nw_gs;
    else
        n=nw;
    end

    colors_gsl={[0 0 80]/255; 
                [0 0 140]/255; 
                [0 0 240]/255;
                [75 75 255]/255; 
                [236 236 0]/255;
                [255 198 105]/255;
                [232 155 0]/255;
                [142 94 0]/255;
                };
    % ----- Anomalies ------
    hold on;
    nlos=9;
    for lgs=3:nlos
        amask=permute(repmat(lgs_map(:,:,dss:dse),[1 1 1 s3 s4]),[1 2 4 5 3]);
        amask(amask~=lgs)=nan;amask(amask==lgs)=1;
        dtmp=mon_anm_ratio(:,:,:,:,dss:dse).*amask;
        dtmp(isinf(dtmp))=nan;

        x=1:12;
        y(1:12,1:(dse-dss+1))=squeeze(nanmean(nanmean(nanmean(dtmp(:,:,1:12,:,1:(dse-dss+1)),1),2),4))*100; % dimension: month and ds
        n=0;pt_pos_low=nan;pt_pos_up=nan;
        for m=1:12
            if sum(isnan(y(m,:)))==0 && sum((y(m,:)))>0

                n=n+1;
                pt_pos_low(n,1)=m;
                pt_pos_low(n,2)=min(y(m,:));
                
                pt_pos_up(n,1)=m;
                pt_pos_up(n,2)=max(y(m,:));

            end
        end
        
        if ~isnan(pt_pos_low)
            v1 = [pt_pos_low; flipud(pt_pos_up)];
    %         f1 = 1:(size(pt_pos_low,1))*2;
            f1 = 1:n*2;
            pat(lgs)=patch('Faces',f1,'Vertices',v1,'FaceColor',colors_gsl{lgs-1},'FaceAlpha',.3, 'EdgeColor','none');

            lines(lgs)=plot(x,nanmean(y,2),'color',colors_gsl{lgs-1},'LineWidth',2.5);
        end
    
    end
    
    % resize subplot
    gca_pos = get(gca,'Position');
    set(gca,'Position',gca_pos.*[1.05 1.35 0.9 0.96]);
    
    cbh=colorbar(gca,'southoutside');
    set(cbh,'visible','off');
    
    ax=gca;
    ax.YAxis(1).Color = 'k';
    ylabel('Anomalies percentage,\gamma_G_S(%)');
    ylim([5 45]);
    
    xlim([1 10.5]);
    xlabel('growing season months');
    xticks([1.5:1:9.5]);
    xticklabels({'          SOS','          +1','          +2','          +3','          +4','          +5','          +6','          +7','          +8'});

    box on;
    title('Anomalies percentage in GS','FontSize',11);

    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');

    legend(lines(3:nlos),{'LOS(3)','LOS(4)','LOS(5)','LOS(6)','LOS(7)','LOS(8)','LOS(9)'});
    legend('Location','NorthEast');
    legend boxoff;

    % ====== Monthly corr. with total anomalies and anomalies ratio ======
    % -------------------------
    %     monthly control
    % -------------------------
    % only include the significant point
	axes(ha(3));

    bg_show=nanmean(dom_ssn_map(:,:,dss:dse),3)';
    bg_show(isnan(bg_show))=0;
    bg_show(isinf(bg_show))=0;
    geoplot(ax_show, bg_show);
    
    cb_intv=0.05; cb_low=-0.6; cb_up=0.6;
    mycolor=(colormap(gca,cbrewer('div','PRGn',(cb_up-cb_low)/cb_intv+1)));
    colormap(gca,mycolor);
    caxis([cb_low-cb_intv/2 cb_up+cb_intv/2]);
    cbh2=colorbar('southoutside');
    set(cbh2, 'AxisLocationMode','manual');
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'c)','FontSize',11,'FontName','Arial');
    title('Early GS vs. late GS','FontSize',11);
    resizeCB(cbh2, cb_xR, cb_yR, cb_wR, 0.6, '\Delta\lambda',0,30,9);
end

