function [fit_ribOcup,conditionRibosomes] = numberOfRibosomes(intensityVector_0,intensityVector_1,folderName,plottingCondition )
% Function to calculate the ribosome occupancy
numberOfProbes =12;
numberOfRibosomes_0=intensityVector_0./numberOfProbes;
numberOfRibosomes_1 =intensityVector_1./numberOfProbes;
%detectionLimit = 1;
detectionLimit = 3*12; % Experimental limit that is detected by the microscope.

%% making zero the spots with intensity in both frames.
numberOfRibosomes_0(intensityVector_1>=detectionLimit|intensityVector_0<detectionLimit)=0;
numberOfRibosomes_1(intensityVector_0>=detectionLimit|intensityVector_1<detectionLimit)=0;

%% Experimental data values
expRibosome_0 = 5.32;
expRibosome_1 = 3.12;

err_expRibosome_0 = 0.32;
err_expRibosome_1 = 0.74;

exp_repetitions_0F =927;
exp_repetitions_1F =15;

%% Removing zero values from the calculation of the number of ribosomes
sim_0F = mean (nonzeros(numberOfRibosomes_0));
sim_1F = mean (nonzeros(numberOfRibosomes_1));
sim_0F(isnan(sim_0F)==1)=0;
sim_1F(isnan(sim_1F)==1)=0;
% BF_0= mean(nonzeros(intensityVector_1(intensityVector_1>0)/12));
% BF_1= mean(nonzeros(intensityVector_1(intensityVector_1>0)/12));

err_sim_0F = std (nonzeros(numberOfRibosomes_0))/sqrt(exp_repetitions_0F);
err_sim_1F = std (nonzeros(numberOfRibosomes_1))/sqrt(exp_repetitions_1F);

if (sim_0F>0&& sim_0F<=6.5)  && (sim_1F>0 && sim_1F<=4)
    conditionRibosomes =1;
else
    conditionRibosomes=0;
end

%% Caclulating the objective function
fit_ribOcup = (abs(expRibosome_0-sim_0F ) + abs(expRibosome_1-sim_1F) )/ (expRibosome_0+expRibosome_1);
%fit_ribOcup = (( (expRibosome_0-sim_0F)^2)/err_sim_0F^2 + ((expRibosome_1-sim_1F)^2)/err_sim_1F^2 );

%% Plotting
if plottingCondition ==1
    
    %% Calculating values for plotting
    ys= [expRibosome_0, sim_0F ;expRibosome_1, sim_1F ];
    y_sem= [err_expRibosome_0,err_sim_0F; err_expRibosome_1,err_sim_1F];
    nr = [expRibosome_0+ err_expRibosome_0, sim_0F+err_sim_0F ,expRibosome_1+err_expRibosome_1, sim_1F +err_sim_1F];
    
    
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2.2, 2.3]; % [left bottom width height]
    num = 2; %number of different subcategories
    c = 1:num;
    axes1 = axes;
    % xlabel('Frames','FontSize',20, 'FontName', 'Arial');
    ylabel('Ribosomes','FontSize',20, 'FontName', 'Arial');
    hold on
    % Bar(s)
    barColor1 = [0.1 0.1 0.1];     %black
    barColor2 = [0.8 0.8 0.8];     %grey
    for i = 1:num
        bar(c(i)-0.2,ys(i,1),0.3,'FaceColor',barColor1);
        bar(c(i)+0.2,ys(i,2),0.3,'FaceColor',barColor2);
    end
    box on
    set(gca,'linewidth',2)
    yticks([1 5 10]);
    yticklabels({'1','5','10'});
    errH1 = errorbar(c-0.2,ys(:,1),y_sem(:,1),'.','Color','k');
    errH2 = errorbar(c+0.2,ys(:,2),y_sem(:,2),'.','Color','r');
    errH1.LineWidth = 0.7;
    errH2.LineWidth = 0.7;
    errH1.Color = [0.5 0.5 0.5];
    errH2.Color = [0.5 0.5 0.5];
    % Set x-ticks
    set(axes1,'Ylim',[0 max(nr)+1]);
    set(axes1,'Xlim',[0.5 2.5]);
    set(axes1,'XTick',[1 1.5 2 2.5],'XTickLabel',...
        {'0F',' ','-1F',' '});
    set (gca ,'FontSize',18); set(gca, 'FontName', 'Arial')
    nameplot = horzcat('No_Ribosomes');
    print('-dpng','-r600',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end
end
