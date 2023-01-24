function []=SOS_EOS_LGS()
    
    global nyr nds mdi;
    global DATA_VI_05rs_sm_out;
    global DATA_CRU_05rs_out;
    global sos_map eos_map lgs_map;
    global ax_show;
    global mask_gs;

    datain=DATA_VI_05rs_sm_out(:,:,:,:,:,mdi); %

    
    [s1 s2 s3 s4 s5]=size(datain);
    dtmp_map_sm=nan(s1,s2,12,nyr,nds);
	sos_map=nan(s1,s2,nds);eos_map=nan(s1,s2,nds);lgs_map=nan(s1,s2,nds);
    
    dtmp_mean=repmat(squeeze(nanmean(datain,3)),[1 1 1 1 12]); % all months average
    dtmp_mean=permute(dtmp_mean,[1 2 5 3 4]);
    dtmp_map_sm=datain-dtmp_mean;
    
    for dsid=1:nds
        
        fprintf(sprintf('==> DSID: %d, calculating SOS ...\n',dsid));
       
        clm_ssn=squeeze(nanmean(dtmp_map_sm(:,:,:,:,dsid),4)); % all years average
        
        dtmp_mean_ds=squeeze(nanmean(DATA_VI_05rs_sm_out(:,:,:,:,dsid,mdi),4));
        clm_ssn_max=squeeze(nanmax(nanmax(nanmax(dtmp_mean_ds))));
        fprintf(sprintf('--- MaxVI: %f, \n',clm_ssn_max));
        
        for i=1:s1
            for j=1:s2
                
                assn=squeeze(clm_ssn(i,j,:));
                
                if nansum(isnan(assn))==12 % skip the ocean
                    continue;
                end
                
                if nanmax(dtmp_mean_ds(i,j,:))<=clm_ssn_max*0.10 % skip bare land
                    continue;
                end
                
                atemp=squeeze(nanmean(DATA_CRU_05rs_out(i,j,:,:),4));
                
                [a_mx a_mxi]=max(assn);
                [a_mn a_mni]=min(assn);
                
                assn(isnan(assn))=a_mn;

                thr_sos=(a_mx-a_mn)*0.1+a_mn;
                thr_eos=(a_mx-a_mn)*0.1+a_mn;

                three_ssn=[assn assn assn];
                three_clm=[atemp atemp atemp];
                
                mark_3ssn=nan(36);
                frs_days=15;
                for m=1:35
                    if three_ssn(m)>=thr_sos && three_ssn(m)<three_ssn(m+1) && three_ssn(m)<=a_mx ... 
                       && three_clm(m)<=frs_days
                   
                       mark_3ssn(m)=1;

                    elseif three_ssn(m)<=a_mx && three_ssn(m)>three_ssn(m+1) && three_ssn(m)>=thr_eos ... 
                           && three_clm(m)<=frs_days
                       
                       mark_3ssn(m)=2;

                    end
                end
                
                % detect SOS
                found_sos=0;m_sos=nan;
                for m=2:36
                    if mark_3ssn(m)==1 && mark_3ssn(m-1)~=1
                        m_sos=m;
                        found_sos=1;
                        break;
                    end
                end
                if ~found_sos
                    fprintf(sprintf('SOS not found!! (i:%d j:%d)\n', i,j));
                    continue;
                end
                
                % detect EOS based on the found SOS month
                found_eos=0;m_eos=nan;
                for m=m_sos:35
                    if mark_3ssn(m)==2 && mark_3ssn(m+1)~=2
                        m_eos=m;
                        found_eos=1;
                        break;
                    end
                end
                if ~found_eos
                    fprintf(sprintf('EOS not found!! (i:%d j:%d)\n', i,j));
                    continue;
                end

                lgs_map(i,j,dsid)=m_eos-m_sos+1;
                
                sos_map(i,j,dsid)=mod(m_sos,12);
                eos_map(i,j,dsid)=mod(m_eos,12);

            end
        end
    end % dsid
    
    mask_tmp=repmat(squeeze(mask_gs(:,:,4)),[1 1 2]);
    
    sos_map(sos_map==0)=nan;sos_map=sos_map.*mask_tmp;
    eos_map(eos_map==0)=nan;eos_map=eos_map.*mask_tmp;
    lgs_map(lgs_map==0)=nan;lgs_map=lgs_map.*mask_tmp;
    
    ds_rng={[1:nds]};
    ax_show.frame=0;
    for dsid=1:length(ds_rng)
        yscal=0.95; xscal=0.02;
        cb_xR=1.0;
        cb_yR=1.6;
        cb_wR=1.0;

        figure('color','w','Position',[587   548   777   354]);
        gap_h=0.005; gap_w=0.002;
        gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
        ha = tight_subplot(1,3,gap,marg_h,marg_w);

        axes(ha(1));
        hold on;
        bg_show=nanmean(sos_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('Start of growing season');
        ctmp=colormap(gca,jet(24));
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb1=colorbar('Southoutside');
        cb1.Ticks=[0:12];
        cb1.TickLabels={'','Jan','','Mar','','May','','Jul','','Sep','','Nov',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
        resizeCB(cb1, cb_xR, cb_yR, cb_wR, 0.6, '',nan,nan,9);

        axes(ha(2));
        hold on;
        bg_show=nanmean(eos_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('End of growing season');
        ctmp=colormap(gca,jet(24));
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb2=colorbar('Southoutside');
        cb2.Ticks=[0:12];
        cb2.TickLabels={'','Jan','','Mar','','May','','Jul','','Sep','','Nov',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');
        resizeCB(cb2, cb_xR, cb_yR, cb_wR, 0.6, '',nan,nan,9);

        axes(ha(3));
        hold on;
        bg_show=nanmean(lgs_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('Length of growing season');
        a=jet(30);
        ctmp=a(1:24,:);
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb3=colorbar('Southoutside');
        cb3.Ticks=[0:12];
        cb3.TickLabels={'','1','','3','','5','','7','','9','','11',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'c)','FontSize',11,'FontName','Arial');
        resizeCB(cb3, cb_xR, cb_yR, cb_wR, 0.6, '(month)',19,10,8);
    end
end
