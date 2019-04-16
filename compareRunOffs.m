function [fit_val_0F, fit_val_1F] = compareRunOffs(pre_namePlot, folderName, rawData_0F, rawData_1F, T_array,sim_0F, sim_1F, sim_Err_0F, sim_Err_1F,plottingCondition,independent_Repetitions)
%% Intensity Deffinition
if independent_Repetitions>1
    mean_sim_0F = mean(sim_0F,1);
    mean_sim_1F = mean(sim_1F,1);
else
    mean_sim_0F = sim_0F;
    mean_sim_1F = sim_1F;
end

std_sim_0F = sim_Err_0F(1,:);
std_sim_1F = sim_Err_1F(1,:);

data_0F = rawData_0F(:,2)';
data_1F = rawData_1F(:,2)';
expTime = rawData_1F(:,4)';

%% Plotting
if plottingCondition ==1
    
    %% PLOTTING. 2 frames square.  USED IN PAPER.
    gray =[0.7,0.7,0.7];
    dark_gray = [0.3,0.3,0.3];
    black = [ 0,0,0];
    green = [0 , 170, 0]./255;
    cyan = [0 , 255,255]./255;
    
    for p=1:2
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        if p==1
            fig1.PaperPosition = [0, 0, 5,4.54]; % [left bottom width height]
        else
            fig1.PaperPosition = [0, 0, 3.34,2.79];
        end
        hold on
        h1 = errorbar(expTime, data_0F, rawData_0F(:,3), '^', 'MarkerEdgeColor',black,'Color',green,'MarkerFaceColor',green,'MarkerSize',11,'LineStyle','none','LineWidth',1.5);
        h2 = errorbar(expTime, data_1F, rawData_1F(:,3), 'o', 'MarkerEdgeColor',black,'Color',cyan,'MarkerFaceColor',cyan,'MarkerSize',11,'LineStyle','none','LineWidth',1.5);
               % PLOTTING 0 FRAME
        lineProps.col= {gray}; lineProps.width = 3;
        A1 = mseb(expTime,mean_sim_0F,std_sim_0F,lineProps,1);
        % PLOTTING -1 FRAME
        lineProps.col = {black}; lineProps.width = 3;
        A2 = mseb(expTime,mean_sim_1F,std_sim_1F,lineProps,1);
        
        box on
        set(gca,'linewidth',2.5)
        % LABLES
        ylim([-0.05 1.4])
        xlim([0 1300]);
        % legend([h1, h2], {'0F','mean(0F,-1F)'})
        grid on
        set(gca,'XGrid','off')
        set(gca,'YGrid','on')
        yticks([0 0.25 0.50 0.75 1.00 1.25]);
        yticklabels({'0.00','0.25','0.50', '0.75', '1.0','1.25'});
        if p ==1
            nameplot = horzcat(pre_namePlot,'_OR');
            set (gca ,'FontSize',18);  set(gca, 'FontName', 'Arial');
            xlabel('Time (sec)','FontSize',22);
            ylabel('Total int. (a.u.)','FontSize',22);
        else
            nameplot = horzcat(pre_namePlot,'_SI_OR');
            set (gca ,'FontSize',18);  set(gca, 'FontName', 'Arial');
            xlabel('Time (sec)','FontSize',20);
            ylabel('Total int. (a.u.)','FontSize',20);
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
    
    
    
    %   save runOffs.mat pre_namePlot folderName expTime data_0F data_1F rawData_0F rawData_1F T_array mean_sim_0F sim_Err_0F sim_Err_1F sim_0F sim_1F repetitions sampling mean_sim_1F std_sim_1F std_sim_0F std_sim_1F
end
%% Comparing simulations and experiments.
%dataPoints = 0:60:1320; dataPoints(1)=1;
dataPoints =30;
%fit_val_0F = sum(abs(data_0F(1:24) - mean_sim_0F(1:24))) / length(dataPoints);
%fit_val_1F = sum(abs(data_1F(1:24) - mean_sim_1F(1:24))) / length(dataPoints);
fit_val_0F = sqrt(sum((data_0F(1:dataPoints) - mean_sim_0F(1:dataPoints)).^2)/ dataPoints);
fit_val_1F = sqrt(sum((data_1F(1:dataPoints) - mean_sim_1F(1:dataPoints)).^2)/ dataPoints);

end
