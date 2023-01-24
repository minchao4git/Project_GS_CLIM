function dom_clminx_plot()
%%
    global v_clminx_r v_clminx_p;
	global ax_show;

    % Pick the climate indices with highest correlation coefficient
    % 6 lags x 10 climate indices

    [s1 s2 s3 s4 s5]=size(v_clminx_r);

    % ===== Pre-calculation =====
    clearvars dtmp d_mx d_mxi dtmp_rs;
    % mask out the insignificant point
    clmi_chosen=[1 2 3 5 6 7 9 ];
    nclmi=length(clmi_chosen);
    dtmp=v_clminx_r(:,:,:,:,clmi_chosen);
    
    % Significant points
    nmn=size(v_clminx_r,3);

    dtmp_rs = reshape(dtmp,[s1 s2 s3 s4*nclmi]);
    [d_mx d_mxi]=nanmax(abs(dtmp_rs),[],4);
    
    d_mxi(isnan(d_mx))=0;

    d_mxi=squeeze(d_mxi);
    map_mxclm_r=nan(s1,s2,nmn);
    for i=1:s1
        for j=1:s2
            for m=1:nmn
                mxi=d_mxi(i,j,m);
                if mxi > 0
                    map_mxclm_r(i,j,m)=dtmp_rs(i,j,m,mxi);
                end
            end
        end
    end
    
%             1     2   3   4(x)   5     6      7    8(x)   9   10(x)
clmi_names={'NAO','EA','WP','EP','PNA','EAWR','SCA','TNH','POL','PT'};
mlabels={'SOS','+2 mon.','+4 mon.'};
%%
if 1==2
	% ===== Dominant map =====
    labels1={'a)','b)','c)'};
    labels2={'d)','e)','f)'};
    defcolormap;

    figure('color','w', 'Position',[135   272   614   467]);
    yscal=0.9; xscal=0.05;
    gap_h=0.01; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.1];
    ha = tight_subplot(2,3,gap,marg_h,marg_w);
    ax_show.frame=0;

    n=0;
    for m=[2 4 6] % IMPARTANT TO START FROM 2 here! because sosmi=2!!
        n=n+1;
        axes(ha(n));

        bg_show=squeeze(map_mxclm_r(:,:,m))';
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
            
        colormap(gca, flipud(cbrewer('div','RdBu',21)));
        caxis([-0.8 0.8]);
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        title(mlabels{n});
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels1{n},'FontSize',11,'FontName','Arial');
        
        if n==3
            cbh2 = colorbar('eastoutside');
            cb_x=cbh2.Position(1); cb_y=cbh2.Position(2);
            cb_w=cbh2.Position(3); cb_h=cbh2.Position(4);
            set(cbh2, 'AxisLocationMode','manual');
            set(cbh2, 'Position',[cb_x+0.04 cb_y cb_w*1.5 cb_h]);
        end
        
        if n==1
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.25,'Max. Corr.(-)','FontSize',12,'FontName','Arial','Rotation',90);
        end
    end
    
    n=0;

    for m=[2 4 6]
        n=n+1;
        axes(ha(3+n));
        
        bg_show=squeeze((squeeze(d_mxi(:,:,m))'+0.5));
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
        
        colormap(gca, mycolormap(1:(nclmi*6)+1,:)/255);
        caxis([0 nclmi*6+1]); % IMPORTANT TO INCLUDE 0 HERE, which corresponds to the white color
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels2{n},'FontSize',11,'FontName','Arial');
        
        if n==3
            cbh1 = colorbar('eastoutside');
            cbh1.Ticks = [0 linspace(1.5, nclmi*6+1.5, nclmi+1)];
            cbh1.TickLabels={'',clmi_names{clmi_chosen'},''};
            cb_x=cbh1.Position(1); cb_y=cbh1.Position(2);
            cb_w=cbh1.Position(3); cb_h=cbh1.Position(4);
            set(cbh1, 'AxisLocationMode','manual');
            set(cbh1, 'Position',[cb_x+0.04 cb_y cb_w*1.5 cb_h]);
        end
        
        if n==1
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.15,'Dominant climate index','FontSize',12,'FontName','Arial','Rotation',90);
        end
    end
end

%%
if 1==1
	% ===== Growing season evolution =====
    global clmzone_grp;
    defcolormap;
    % skip the white color
    mycolormap=mycolormap(2:end,:);
    nlag=6;
    clmz_name={'Boreal','Temperate', 'Medit.'};
	labels={'a)','b)','c)'};
    
	figure('color','w','Position',[-920    61   606   839]);
    yscal=0.90; xscal=0.92;
    gap_h=0.02; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.15 0.1];
    ha = tight_subplot(4,2,gap,marg_h,marg_w);

    n=0;
    thr_show=0.061;
    for clm_grp=1:3
        n=n+1;
        axes(ha(n));
        for clmi=fliplr(clmi_chosen)
            clim_ci = find(clmi_chosen==clmi);
            printed_label_str='';
            for lag=[0:(nlag-1)]
                % get the corresponding climate zone and landclass
                amask=clmzone_grp; amask(amask~=clm_grp)=nan; amask(amask==clm_grp)=1;
                ngrid=nansum(nansum(amask,1),2);
                dtmp=squeeze(v_clminx_r(:,:,16,lag+1,clmi));
                dmask=squeeze(v_clminx_p(:,:,16,lag+1,clmi));
                dtmp(dmask>=0.05)=nan;
                dtmp=dtmp.*amask;

                % Aggregate spatial values and standardized by number of
                % grid point
                sumR_gs=squeeze(sum(sum(dtmp,1,'omitnan'),2,'omitnan'))/ngrid; 
                meanR_gs=squeeze(nanmean(nanmean(dtmp,1),2));
                hold on;
                lincolor=mycolormap((clim_ci-1)*6+(lag+1),:)/255;

                plot([0.85 1.25],[sumR_gs sumR_gs],'color',lincolor, 'LineWidth',2.5);
                xlim([0 1.5]);
                
                % Label the line with large values
                lab_str=sprintf('%s_l_a_g_%d',clmi_names{clmi},lag);
                if clm_grp ==1  &&  (     (strcmp(clmi_names{clmi},'NAO')==true && (lag==0 || lag==1 || lag==2))... 
                                      ||  (strcmp(clmi_names{clmi},'SCA')==true && (lag==2)) ...
                                    )
                    if contains(printed_label_str,lab_str)==false
                        text(0.1, sumR_gs,lab_str,'FontSize',8,'FontName','Arial');
                        printed_label_str=[printed_label_str ' ' lab_str];

                    end
                elseif clm_grp ==2  &&  ( (strcmp(clmi_names{clmi},'NAO')==true && (lag==0 || lag==1))... 
                                      ||  (strcmp(clmi_names{clmi},'SCA')==true && (lag==2)) ...
                                      ||  (strcmp(clmi_names{clmi},'EA')==true && (lag==0 || lag==5)) ...
                                    )
                    if contains(printed_label_str,lab_str)==false
                        text(0.1, sumR_gs,lab_str,'FontSize',8,'FontName','Arial');
                        printed_label_str=[printed_label_str ' ' lab_str];

                    end
                elseif clm_grp ==3  &&  ( (strcmp(clmi_names{clmi},'NAO')==true && (lag==1 || lag==2))... 
                                      ||  (strcmp(clmi_names{clmi},'EA')==true && (lag==0 || lag==1)) ...
                                    )
                    if contains(printed_label_str,lab_str)==false
                        text(0.1, sumR_gs,lab_str,'FontSize',8,'FontName','Arial');
                        printed_label_str=[printed_label_str ' ' lab_str];

                    end
                end
                
                hold off;
                box on;
            end
        end
        if clm_grp==3
             ylim([-0.1 0.18]);
        else
             ylim([-0.1 0.18]);
        end
        ylabel(sprintf('%s \n Standardized Corr.(-)', clmz_name{clm_grp}));
        set(gca,'XTickLabel',{},'Box','on');
        
        gcaPos=get(gca,'Position');
        LenCutR=0.68;
        set(gca,'Position',gcaPos.*[1 1 (1-LenCutR) 1]);
        
        n=n+1;
        axes(ha(n));
        % show the low r lines in grey
        for clmi=fliplr(clmi_chosen)
            for lag=[0:(nlag-1)]
                % get the corresponding climate zone and landclass
                amask=clmzone_grp; amask(amask~=clm_grp)=nan; amask(amask==clm_grp)=1;
                ngrid=nansum(nansum(amask,1),2);
                dtmp=squeeze(v_clminx_r(:,:,1:12,lag+1,clmi));
                dmask=squeeze(v_clminx_p(:,:,1:12,lag+1,clmi));
                dtmp(dmask>=0.05)=nan;
                dtmp=dtmp.*repmat(amask,[1 1 12]);

                % Aggregate spatial values for 12 months and standardized by number of
                sumR_gs=squeeze(sum(sum(dtmp,1,'omitnan'),2,'omitnan'))/ngrid; 
                sumR_gs(1)=nan;
                hold on;
                if sum(abs(sumR_gs)>thr_show)==0
                    plot(1:12,sumR_gs,'color',[0.8 0.8 0.8], 'LineWidth',1);
                end
                hold off;
            end
        end
        
        for clmi=fliplr(clmi_chosen)
            clim_ci = find(clmi_chosen==clmi);
            printed_label_str='';
            for lag=[0:(nlag-1)]
                amask=clmzone_grp; amask(amask~=clm_grp)=nan; amask(amask==clm_grp)=1;
                ngrid=nansum(nansum(amask,1),2);
                dtmp=squeeze(v_clminx_r(:,:,1:12,lag+1,clmi));
                dmask=squeeze(v_clminx_p(:,:,1:12,lag+1,clmi));
                dtmp(dmask>=0.05)=nan;
                dtmp=dtmp.*repmat(amask,[1 1 12]);

                % Aggregate spatial values for 12 months and standardized by number of
                sumR_gs=squeeze(sum(sum(dtmp,1,'omitnan'),2,'omitnan'))/ngrid; 
                meanR_gs=squeeze(nanmean(nanmean(abs(dtmp),1),2));
                sumR_gs(1)=nan;sumR_gs(sumR_gs==0)=nan;meanR_gs(1)=nan;
                hold on;
                lincolor=mycolormap((clim_ci-1)*6+(lag+1),:)/255;

                last_i=1;
                for i=1:length(sumR_gs)
                    if abs(sumR_gs(i))>0
                        last_i=i;
                    end
                end
                if last_i>=11
                    last_i=11;
                end
                sumR_gs(last_i+1)=0;
                
                if sum(abs(sumR_gs)>thr_show)>0
                    plot(1:12,sumR_gs,'color',lincolor, 'LineWidth',2.5);
                    for s=2:12
                        if abs(meanR_gs(s))>0
                           [sz, mkr]=getSZ(meanR_gs(s));
                           scatter(s,sumR_gs(s),sz*1.5,lincolor,'filled',mkr);
                        end
                    end
                end
                
                % Label the line with large values
                for m=1:12
                    lab_str=sprintf('%s_l_a_g_%d',clmi_names{clmi},lag);
                    if clm_grp <=2  && (abs(sumR_gs(m)) > thr_show ... 
                                       || (lag==0 && ( strcmp(clmi_names{clmi},'NAO')==true ...
                                                    || strcmp(clmi_names{clmi},'SCA')==true ...
                                                    || strcmp(clmi_names{clmi},'EA')==true)))
                        if contains(printed_label_str,lab_str)==false
                            text(m+0.5, sumR_gs(m),lab_str,'FontSize',10,'FontName','Arial');
                            printed_label_str=[printed_label_str ' ' lab_str];

                        end
                    elseif clm_grp ==3  && abs(sumR_gs(m)) > 0.04
                        if contains(printed_label_str,lab_str)==false
                            text(m+0.5, sumR_gs(m),lab_str,'FontSize',10,'FontName','Arial');
                            printed_label_str=sprintf('%s %s',printed_label_str,lab_str);
                        end
                    end
                end
                
                hold off;
            end
        end
        xticks([1:1:11.5]);
        set(gca,'XTickLabel',{},'YTickLabel',{},'Box','on');
        
        if clm_grp==3
            xticklabels({'','SOS','+1','+2','+3','+4','+5','+6','+7','+8','+9'});
            xlabel('GS months');
        end
        
        if clm_grp==3
             ylim([-0.1 0.18]);
        else
             ylim([-0.1 0.18]);
        end
        xlim([1.8 12]);

        % Adjust axis size
        gcaPos=get(gca,'Position');
        set(gca,'Position',[gcaPos(1)-gcaPos(3)*(LenCutR*1.0) gcaPos(2) gcaPos(3)*(1+LenCutR) gcaPos(4)]);
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{clm_grp},'FontSize',11,'FontName','Arial');
        
        if clm_grp==3
            axes(ha(n+1));
            gca_img=imagescn(flipud(squeeze(d_mxi(:,:,7))'+0.5));
            colormap(gca, mycolormap(1:(nclmi*6),:)/255);
            cbh1 = colorbar('southoutside');
            caxis([0 nclmi*6]); 
            cbh1.Ticks = [linspace(0.5, nclmi*6+0.5, nclmi+1)];
            cbh1.TickLabels={clmi_names{clmi_chosen'},''};
            cb_x=cbh1.Position(1); cb_y=cbh1.Position(2);
            cb_w=cbh1.Position(3); cb_h=cbh1.Position(4);
            
            set(cbh1, 'Position',[cb_x+0.04 cb_y cb_w*1.5 cb_h]);
            set(gca_img,'Visible','off');
            box off;
            set(gca,'XColor', 'none','YColor','none')
            resizeCB(cbh1, 1, 1.7, 1.2, 1.6, '',nan,nan,1);
            
            gcaPos=get(gca,'Position');
            set(gca,'Position',[gcaPos(1) gcaPos(2) gcaPos(3)*0.1 gcaPos(4)*0.1]);

        end
    end
end
%%    
end

function [sz, mkr]=getSZ(meanR_gs)
    if abs(meanR_gs) > 0 && abs(meanR_gs) <= 0.2
        sz=5;
        mkr='o';
    elseif abs(meanR_gs) > 0.2 && abs(meanR_gs) <= 0.4
        sz=25;
        mkr='s';
    elseif abs(meanR_gs) > 0.4 && abs(meanR_gs) <= 0.6
        sz=35;
        mkr='d';
    elseif abs(meanR_gs) > 0.6 && abs(meanR_gs) <= 0.8
        sz=65;
        mkr='*';
    elseif abs(meanR_gs) > 0.8 && abs(meanR_gs) <= 1
        sz=85;
        mkr='x';
    else
        sz=1;
        mkr='xxxxxx';
    end
end

