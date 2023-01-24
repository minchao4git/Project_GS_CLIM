function Mcolor(res)
   
    nTri_ASide=res;
    dTri_l=1/nTri_ASide;
    figure('color','w','Position',[740   174   399   324]);
    hold on;
    for i = 1:nTri_ASide
        for j=1:nTri_ASide-(i-1)

            % Draw the triangle matrix, define the x and y coordinates
            TP1=[(i-1)*dTri_l+0.5*(j-1)*dTri_l              (j-1)*(sqrt(3)/2)*dTri_l];
            TP2=[(i)*dTri_l+0.5*(j-1)*dTri_l                (j-1)*(sqrt(3)/2)*dTri_l];
            TP3=[(i-1)*dTri_l+(1/2)*dTri_l+0.5*(j-1)*dTri_l (j)*(sqrt(3)/2)*dTri_l];

            % Normal stand triangle
            v1 = [TP1; TP2; TP3];
            f1 = [1 2 3];

            wr=(nTri_ASide-(i-1) - (j-1));
            wg= (i-1);
            wb= (j-1);
            [RGB]=WeightToRGB(wr, wg, wb);

            patch('Faces',f1,'Vertices',v1,'FaceColor',RGB,'EdgeColor','none');

            % Headstand triangle
            if j < nTri_ASide-(i-1)
                TP4=[(i-1)*dTri_l+(1/2)*dTri_l+0.5*(j-1)*dTri_l+1*dTri_l (j)*(sqrt(3)/2)*dTri_l];
                v1 = [TP3; TP2; TP4];
                f1 = [1 2 3];

                wr1=wr-1;
                wg1=wg+0.5;
                wb1=wb+0.5;
                [RGB1]=WeightToRGB(wr1, wg1, wb1);

                patch('Faces',f1,'Vertices',v1,'FaceColor',RGB1,'EdgeColor','none');
                
            end
        end
    end
    xlim([-0.15 1.15]);
    ylim([-0.15 1.05]);
    box off;
    set(gca,'XColor','none','YColor','none','TickDir','out');

    text(-0.13, -0.1,'Temp.','FontSize',24,'FontName','Arial');
    text( 0.9, -0.1,'Rad.','FontSize',24,'FontName','Arial');
    text( 0.25,  0.93,'Soil Moist.','FontSize',24,'FontName','Arial');
    
    hold off;
end
