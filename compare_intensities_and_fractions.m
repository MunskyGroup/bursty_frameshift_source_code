function    [ fit_meanIntensity,fit_ratio,conditionFractions] = compare_intensities_and_fractions (intensityVector_0,intensityVector_1,plottingCondition,folderName)

backgroundDetectionLimit =24;
numberOfProbes =12;
%% Spot classification
I_y_0 = intensityVector_0(intensityVector_0 >= backgroundDetectionLimit & intensityVector_1 <= backgroundDetectionLimit);
I_p_m1 = intensityVector_1(intensityVector_1 >= backgroundDetectionLimit & intensityVector_0 <= backgroundDetectionLimit);
I_w_0 = intensityVector_0(intensityVector_1 >= backgroundDetectionLimit & intensityVector_0 >= backgroundDetectionLimit);
I_w_m1 = intensityVector_1(intensityVector_1 >= backgroundDetectionLimit & intensityVector_0 >= backgroundDetectionLimit);

n_I_y_0 = size(I_y_0,1);
n_I_w_0 = size(I_w_0,1);
n_I_p_m1 = size(I_p_m1,1);
n_I_w_m1 = size(I_w_m1,1);

%% Comparing experimental and simulation intensity mean values
means_mod = [mean(I_y_0), mean(I_w_0),mean(I_w_m1),mean(I_p_m1)]/numberOfProbes;
means_dat = [6.3, 4.3, 2, 3.1];
spot_numbers =[ 367,9,8,1197];
spots_BF= spot_numbers(2)+spot_numbers(3);
spot_intens_sem_data = [0.3842, 1.0754 , 0.5276, 0.6472];

sum_IntExpData = sum(means_dat);

spot_intens_sem_model= [std(I_y_0)/sqrt(n_I_y_0), std(I_w_0)/sqrt(n_I_w_0),std(I_w_m1)/sqrt(n_I_w_m1),std(I_p_m1)/sqrt(n_I_p_m1)];

fit_meanIntensity = abs(sum(means_dat-means_mod))/sum_IntExpData;
fit_meanIntensity(isnan(fit_meanIntensity)==1)= 1;

%% Comparing experimental and simulation intensity ratios
total_b2g_ratio_mod = (sum(I_w_m1)+sum(I_p_m1))/(sum(I_y_0)+sum(I_w_0));
total_b2g_ratio_dat = 4.1065e-02;
total_b2g_ratio_sd = 4.16e-16;

[mean_ratio_sim,sd_ratio_sim]  = bootstrp_Intensities (intensityVector_0,intensityVector_1);
fit_ratio= abs(mean_ratio_sim- total_b2g_ratio_dat)/total_b2g_ratio_dat;
if fit_ratio < 0.5 && fit_meanIntensity<0.4
    conditionFractions=1;
else
    conditionFractions=0;
end

if plottingCondition==1
    
    for p=1:2
        %% Plotting ratio
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        if p ==1
            %    fig1.PaperPosition = [0, 0, 2.2, 2.3]; % [left bottom width height]
            fig1.PaperPosition = [0, 0, 2,4.54];
        else
            fig1.PaperPosition = [0, 0,1.67, 2.79];
        end
        num = 1; %number of different subcategories
        c = 1:num;
        axes1 = axes;
        barColor1 = [0.1 0.1 0.1];     %black
        barColor2 = [0.8 0.8 0.8];     %grey
        hold on
        %bar([total_b2g_ratio_dat, total_b2g_ratio_mod])
        for i = 1:num
            bar(1,total_b2g_ratio_dat,0.8,'FaceColor',barColor1);
            bar(2,mean_ratio_sim,0.8,'FaceColor',barColor2);
        end
        errH1 = errorbar(1,total_b2g_ratio_dat,total_b2g_ratio_sd,'.','Color',barColor2);
        errH2 = errorbar(2,mean_ratio_sim,sd_ratio_sim,'.','Color',barColor1);
        errH1.LineWidth = 0.7;
        errH2.LineWidth = 0.7;
        errH1.Color = barColor2;
        errH2.Color = barColor1;
        max_Y = max([total_b2g_ratio_dat, mean_ratio_sim]);
        box on
        set(gca,'linewidth',2)
        set(gca,'xtick',[0.7 2.2],'xticklabels',{'Da',' Mo'},'ylim',[0 max_Y*1.2],'xlim',[0,3])
        if p ==1
            nameplot = horzcat('ratio');
            set (gca ,'FontSize',18); set(gca, 'FontName', 'Arial')
            ylabel('FS : nonFS intensity ratio','FontSize',20,'FontName', 'Arial')
        else
            nameplot = horzcat('ratio_SI');
            set (gca ,'FontSize',16); set(gca, 'FontName', 'Arial')
            ylabel('FS : nonFS intensity','FontSize',18,'FontName', 'Arial')
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
    
    
    %% Plotting Intensities
    
    for p=1:2
        figure('visible', 'off');
        fig2= gcf;
        fig2.PaperUnits = 'inches';
        if p ==1
            fig2.PaperPosition = [0, 0,6.5, 3.13]; % [left bottom width height]
        else
            fig2.PaperPosition = [0, 0,2.5, 2.79];
        end
        num = 4; %number of different subcategories
        c = 1:num;
        axes1 = axes;
        hold on
        box on
        set(gca,'linewidth',2)
        % Bar(s)
        barColor1 = [0.1 0.1 0.1];     %black
        barColor2 = [0.8 0.8 0.8];     %grey
        for i = 1:num
            b1= bar(c(i)-0.2,means_dat(i),0.3,'FaceColor',barColor1);
            b2= bar(c(i)+0.2,means_mod(i),0.3,'FaceColor',barColor2);
        end
        errH1 = errorbar(c-0.2,means_dat,spot_intens_sem_data,'.','Color',barColor2);
        errH2 = errorbar(c+0.2,means_mod,spot_intens_sem_model,'.','Color',barColor1);
        errH1.LineWidth = 0.7;
        errH2.LineWidth = 0.7;
        errH1.Color = barColor2;
        errH2.Color = barColor1;
        max_Y = max([means_mod+ spot_intens_sem_model, means_dat+ spot_intens_sem_data]);
        box on
        set(gca,'linewidth',2)
        ylim([0 max_Y*1.25]);
        set(axes1,'Xlim',[c(1)-0.5 c(end)+0.5]);
        set(axes1,'XTick',[c],'XTickLabel',...
            {'0F','BF_{0}','BF_{-1}','-1F'});
        %lgd=legend([b1,b2],'Data','Model');
        %set(lgd,'FontSize',10);
        %lgd.Location='NorthEast';
        NumTicks=3;
        yT = get(gca,'YLim');
        set(gca,'YTick',round(linspace(0,max_Y,NumTicks),0))
        if p ==1
            nameplot = horzcat('Intensities');
            set (gca ,'FontSize',20); set(gca, 'FontName', 'Arial');
            ylabel('Intensity  (UMP)','FontSize',22,'FontName', 'Arial');
        else
            nameplot = horzcat('Intensities_SI');
            set (gca ,'FontSize',16); set(gca, 'FontName', 'Arial');
            ylabel('Intensity (UMP)','FontSize',18,'FontName', 'Arial');
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
    
end


