function  [fit_Percent,conditionPercentage] = percentagePerFrame(numberOfSpots,folderName,intensityVector_0,intensityVector_1,plottingCondition )
%trayectoryAnalysis. This function calculates distributions for  FS signals.  % FS % 0F % -1F % 0 and -1F.
%% Number of trajectories with 0F, -1F and 0F & -1F.
for i =1:numberOfSpots
    StrData(1,i).trajectory_1 = intensityVector_0(i,:);
    StrData(1,i).trajectory_2 = intensityVector_1(i,:);
    % calculating size of trajectories
    StrData(1,i).Length = length(StrData(1,i).trajectory_1);
end
for i =1:numberOfSpots
    if max( StrData(1,i).trajectory_1) >0
        StrData(1,i).trajectory_1 =  StrData(1,i).trajectory_1./max( StrData(1,i).trajectory_1);
    end
    if max( StrData(1,i).trajectory_2) >0
        StrData(1,i).trajectory_2 = StrData(1,i).trajectory_2 ./ max(StrData(1,i).trajectory_2);
    end
end
observationTime = 0; % in seconds
[sim_0F,sim_1F, sim_BF,err_sim_0F,err_sim_1F,err_sim_BF]  = simulationsToPercentage (numberOfSpots,StrData,observationTime);

%% Experimental Data
exp_0F = 95.6;
exp_1F = 2.34;
exp_BF = 2.08;

err_exp_0F = 1.3;
err_exp_1F = 0.5;
err_exp_BF = 1.1;

%% Evaluating objective function
fit_Percent = abs(log10(sim_0F/ exp_0F)) + abs(log10(sim_1F/exp_1F))+ abs(log10(sim_BF / exp_BF));

%% Section that reject parameters if those are not in a biological feasible range
if fit_Percent<1
    conditionPercentage =1;
else
    conditionPercentage=0;
end

%%
if plottingCondition ==1
    ys = [exp_0F, sim_0F ;exp_1F, sim_1F; exp_BF, sim_BF ];
    y_sem = [err_exp_0F,err_sim_0F; err_exp_1F,err_sim_1F; err_exp_BF,err_sim_BF ];
    
    for p=1:2
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        if p ==1
            fig1.PaperPosition = [0, 0,4, 4.54]; % [left bottom width height]
        else
            fig1.PaperPosition = [0, 0,2.5, 2.79];
        end
        num = 3; %number of different subcategories
        c = 1:num;
        axes1 = axes;
        hold on
        box on
        set(gca,'linewidth',2)
        % Bar(s)
        barColor1 = [0.1 0.1 0.1];     %black
        barColor2 = [0.8 0.8 0.8];     %gray
        for i = 1:num
            b1=bar(c(i)-0.2,ys(i,1),0.3,'FaceColor',barColor1);
            b2=bar(c(i)+0.2,ys(i,2),0.3,'FaceColor',barColor2);
        end
        errH1 = errorbar(c-0.2,ys(:,1),y_sem(:,1),'.','Color','k');
        errH2 = errorbar(c+0.2,ys(:,2),y_sem(:,2),'.','Color','k');
        errH1.LineWidth = 0.8;
        errH2.LineWidth = 0.8;
        errH1.Color = [0.5 0.5 0.5];
        errH2.Color = [0.5 0.5 0.5];
        box on
        yticks([1 50 100]);
        yticklabels({'1','50','100'});
        % Set x-ticks
        set(axes1,'Xlim',[0.5 3.7]);
        set(axes1,'XTick',[1 1.5 2 2.5 3 3.5],'XTickLabel',...
            {'0F ','  ',' -1F ','  ',' BF',' '});
        if p ==1
            nameplot = horzcat('Percent');
            set (gca ,'FontSize',20);
            ylabel('Translation %','FontSize',22, 'FontName', 'Arial');
            
            lgd=legend([b1,b2],'Data','Model');
            set(lgd,'FontSize',14);
            lgd.Location='NorthEast';
        else
            nameplot = horzcat('Percent_SI');
            set (gca ,'FontSize',18);
            ylabel('Translation %','FontSize',20, 'FontName', 'Arial');
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
end
end
