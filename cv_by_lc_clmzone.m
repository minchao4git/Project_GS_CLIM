function cv_by_lc_clmzone()

    % ===== CV by landclass and climate zone ===== 
    global lc_dom_grp clmzone_grp mon_cv;
    global ds_s ds_e nds;
    global mon_anm_ratio;
    global clm_cl_x;

    [s1 s2]=size(clmzone_grp);

    % Calculation
    data_calc=nan(s1,s2,5);
    clm_cl_x=nan(3,8);

    % number of gridpoint
    clm_cl_cnt=nan(3,8,5);

    % group in bar chart
    clm_cl_rpt=nan(3,8,5);
    clm_cl_stderr=nan(8,3,5);
    clm_lc_p=nan(8,3,5);

    dss=ds_s;dse=ds_e;
    
    data_calc(:,:,1,1:nds)=nanmean(mon_anm_ratio(:,:,13,:,dss:dse),4); % 1st half
    data_calc(:,:,2,1:nds)=nanmean(mon_anm_ratio(:,:,14,:,dss:dse),4); % 2nd half
    data_calc(:,:,3,1:nds)=nanmean(mon_cv(:,:,16,dss:dse),3); % 16 for early GS
    data_calc(:,:,4,1:nds)=nanmean(mon_cv(:,:,17,dss:dse),3); % 17 for late GS
    data_calc(:,:,5,1:nds)=nanmean(mon_cv(:,:,15,dss:dse),3); % 15 for annual value

    for clm=1:3
        for lc=1:8

            % spatial points
            amask=clmzone_grp; amask(amask~=clm)=nan;amask(amask==clm)=1;
            bmask=lc_dom_grp; bmask(bmask~=lc)=nan;bmask(bmask==lc)=1;
            dtmp=data_calc.*repmat(amask,[1 1 5 nds]).*repmat(bmask,[1 1 5 nds]);

            for v=1:5
                % dataset level
                for dsi=1:nds
                    atmp=squeeze(dtmp(:,:,v,dsi));
                    atmp=atmp(~isnan(atmp));

                    % Mean of stability across space
                    clm_cl_rpt(clm,lc,v,dsi)=nanmean(atmp);
                end
                
                % mean across dataset
                btmp=squeeze(nanmean(dtmp(:,:,v,:),4));
                btmp=btmp(~isnan(btmp));
                
                % number of gridpoints
                clm_cl_cnt(clm,lc,v)=length(btmp);
                % Standard error across space
                clm_cl_stderr(clm,lc,v)=std(btmp,0,1,'omitnan');

                % Confidence interval across space
                intv = tinv(0.975, length(btmp)-1); % Compute 95% Confidence Standard Deviation
                clm_lc_p(clm,lc,v)= clm_cl_stderr(clm,lc,v)*intv;
            end
            a=nanmean(dtmp(:,:,3,:),4);a=a(~isnan(a));
            b=nanmean(dtmp(:,:,4,:),4);b=b(~isnan(b));
            clm_cl_sig(clm,lc)=ttest2(a(:),b(:));
            
        end
    end

    clm_cl_sig
    
    % ===== plotting ===== 
    color_lc_grp=[[0 129 27];
                      [124 255 100]; 
                      [142 204 51];
                      [214 236 163]; 
                      [244 181 120]; 
                      [0 104 150]; 
                      [255 236 131]; 
                      [170 255 255]
        ]/255;

    lc_names={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
    clmz_name={'Boreal','Temperate', 'Medit.'};
    labels={'a)','b)','c)'};
    nv_show=1:4;
    figure('Position',[2475 247 467 499],'Color','w');

    t = tiledlayout(3,1,'TileSpacing','Compact');
    % --- land class ---
    for clm_grp=1:3

        [B sort_i(:,clm_grp)]=sort(nanmean(clm_cl_rpt(clm_grp,:,5,:),4),'descend','MissingPlacement','last');
        xtmp=sort_i(:,clm_grp);

        % --- climatezone A ---
        nexttile;

        n=0;
        x=nan;
        for i=1:length(xtmp)
          if ~isnan(nanmean(clm_cl_rpt(clm_grp,xtmp(i),5,:),4)) && clm_cl_cnt(clm_grp,xtmp(i),5) >= 10 % skim based on the CV values
             n=n+1;
             x(n)=xtmp(i);
             clm_cl_x(clm_grp,n)=xtmp(i); % output for decomp_VegClm_relation_plot
          end
        end

        hold on
        width=0.8;
        xoff=0;

        b1=bar((1:length(x))-xoff, squeeze(nanmean(clm_cl_rpt(clm_grp,x,[5 3 4],:),4)),width);
        b1(2).FaceColor = [22 170 116]/255;
        b1(3).FaceColor = [255 232 99]/255;
        
        % First error bar in the group
        x_bar=[b1(1,1).XEndPoints]';
        y_bar=[b1(1,1).YEndPoints; b1(1,2).YEndPoints; b1(1,3).YEndPoints]';
        er = errorbar(x_bar, squeeze(nanmean(clm_cl_rpt(clm_grp,x,[5],:),4)), squeeze(clm_lc_p(clm_grp,x,[5])));
        for v=1:1
            er(1,v).Color = [0 0 0];
            er(1,v).LineStyle = 'none'; % remove the plot line
        end
        
        % The second and the third bar in the group
        x_bar_2_3=[b1(1,2).XEndPoints; b1(1,3).XEndPoints]';

        for i=1:length(x)
            er_2_3 = errorbar(x_bar_2_3(i,:), squeeze(nanmean(clm_cl_rpt(clm_grp,x(i),[3 4],:),4)), squeeze(clm_lc_p(clm_grp,x(i),[3 4])));
            er_2_3.Color = [0 0 0];
            er_2_3.LineStyle = 'none'; % remove the plot line
            
            % solid line for significant values
            if clm_cl_sig(clm_grp,x(i))==0
                er_2_3.Bar.LineStyle = 'dashed';
            elseif clm_cl_sig(clm_grp,x(i))==1
                er_2_3.Bar.LineStyle = 'solid';
            end
        end
        
        % NDVI
        bar_off=0.22;
        sz=20;sz1=30;
        scatter((1:length(x))-xoff-bar_off,squeeze(clm_cl_rpt(clm_grp,x,[5],1)),sz,'d','MarkerEdgeColor',[0 0 0]);
        scatter((1:length(x))-xoff,squeeze(clm_cl_rpt(clm_grp,x,[3],1)),sz,'d','MarkerEdgeColor',[0 0 0]);
        fig_s1=scatter((1:length(x))-xoff+bar_off,squeeze(clm_cl_rpt(clm_grp,x,[4],1)),sz,'d','MarkerEdgeColor',[0 0 0]);
        
        % EVI2
        scatter((1:length(x))-xoff-bar_off,squeeze(clm_cl_rpt(clm_grp,x,[5],2)),sz1,'x','MarkerEdgeColor',[0 0 0]);
        scatter((1:length(x))-xoff,squeeze(clm_cl_rpt(clm_grp,x,[3],2)),sz1,'x','MarkerEdgeColor',[0 0 0]);
        fig_s2=scatter((1:length(x))-xoff+bar_off,squeeze(clm_cl_rpt(clm_grp,x,[4],2)),sz1,'x','MarkerEdgeColor',[0 0 0]);
        
        ax=gca;
        ax.YAxis(1).Color = 'k';
        ylabel(sprintf('%s \n anomalies (-)',clmz_name{clm_grp}));
        set(gca,'Fontsize',8);

        set(gca, 'XTick', 1:length(x));
        set(gca, 'XTickLabel',{lc_names{x}});
        xtickangle(0);

        plot([-100 100],[0 0],'-k');

        hold off
        box on;
        clmz_name={'Boreal','Temperate', 'Medit.'};
        ylabel(sprintf('%s (CV)',clmz_name{clm_grp}));

        if clm_grp==1
           title('CV by climate zone and land class','Fontsize',12);
        end

        if clm_grp==3
            xlabel('Landclass');
        end
        xlim([0.2 6.8]);
        if clm_grp==1
            ylim([0 0.18]);
        elseif clm_grp==2
            ylim([0 0.13]);
        else
            ylim([0 0.18]);
        end

        if clm_grp==1
            legend([fig_s1 fig_s2],'GIMMS3g-NDVI','VIP-EVI2','Location','Best');
            legend boxoff;
        end
        
        if clm_grp==3
            legend(b1,{'Entire GS', 'Early GS', 'Late GS'});
            legend boxoff;
        end
        % labels for sub-plots
        yscal=0.92; xscal=0.03;
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{clm_grp},'FontSize',11,'FontName','Arial');
        
    end
end