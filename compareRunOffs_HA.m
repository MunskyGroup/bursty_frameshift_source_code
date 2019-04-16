function [fit_HA_FS,fit_HA_NonFs,fit_Flag] = compareRunOffs_HA(pre_namePlot, folderName,rawData_SunTag_inFS,rawData_HA_inFS,rawData_HA_inNonFS,rawData_FLAG_Only, T_array, sim_HA_FS, sim_1F_FS, sim_FLAG,sim_HA_NonFS, sim_Err_HA_FS, sim_Err_1F_FS, sim_Err_HA_NonFS,sim_Err_Flag, plottingCondition,independent_Repetitions)
%% Intensity vector deffinition
numberSpots = size(sim_HA_FS,2);

if independent_Repetitions>1
    sim_HA_FS = sim_HA_FS(any(sim_HA_FS,2),:);
    sim_1F_FS = sim_1F_FS(any(sim_1F_FS,2),:);
    sim_FLAG = sim_FLAG(any(sim_FLAG,2),:);
    sim_HA_NonFS = sim_HA_NonFS(any(sim_HA_NonFS,2),:);
end

if isempty (sim_HA_FS)==1
    sim_HA_FS = zeros(1, numberSpots);
end
if isempty (sim_1F_FS)==1
    sim_1F_FS = zeros(1, numberSpots);
end
if isempty (sim_FLAG)==1
    sim_FLAG = zeros(1, numberSpots);
end
if isempty (sim_HA_NonFS)==1
    sim_HA_NonFS = zeros(1, numberSpots);
end

if independent_Repetitions>1
    if size(sim_HA_FS,1)>1
        mean_sim_HA_FS = mean(sim_HA_FS,1);
    else
        mean_sim_HA_FS = sim_HA_FS;
    end
    
    if size(sim_1F_FS,1)>1
        mean_sim_1F_FS = mean(sim_1F_FS,1); % FLAG ONLY
    else
        mean_sim_1F_FS = sim_1F_FS;
    end
    
    if size(sim_FLAG,1)>1
        mean_sim_Flag =mean(sim_FLAG,1);
    else
        mean_sim_Flag = sim_FLAG; % FLAG ONLY
    end
    
    if size(sim_HA_NonFS,1)>1
        mean_sim_HA_NonFS =mean(sim_HA_NonFS,1);
    else
        mean_sim_HA_NonFS =sim_HA_NonFS;
    end
    
end

std_sim_HA_FS = sim_Err_HA_FS(1,:);
std_sim_1F_FS = sim_Err_1F_FS(1,:);
std_sim_HA_NonFS = sim_Err_HA_NonFS(1,:); % FLAG ONLY
std_sim_FLAG =  sim_Err_Flag(1,:);

%% Loading exprimental data
data_HA_FS = rawData_HA_inFS(:,2)';
data_1F = rawData_SunTag_inFS(:,2)';
data_0F_NonFS = rawData_HA_inNonFS(:,2)';
data_Flag = rawData_FLAG_Only(:,2)';
expTime = rawData_HA_inFS(:,4)';

%% Plotting
if plottingCondition ==1
    
    %% PLOTTING. Run-Off Time 3 Frames Square
    dark_gray =[0.6,0.6,0.6];
    gray = [0.2,0.2,0.2];
    black = [ 0,0,0];
    orange = [255, 170, 85]./255;
    for p=1:2
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        if p==1
            fig1.PaperPosition = [0, 0, 5,4.54]; % [left bottom width height]
        else
            fig1.PaperPosition = [0, 0, 3.34,2.79];
        end
        % PLOTTING
        hold on
        h1=  errorbar(expTime, data_HA_FS, rawData_HA_inFS(:,3), 'o','Color',orange,'MarkerEdgeColor',black,'MarkerFaceColor',orange,'MarkerSize',12,'LineStyle','none','LineWidth',1.5);
        h2 = errorbar(expTime, data_0F_NonFS, rawData_HA_inNonFS(:,3), '^','Color',orange,'MarkerEdgeColor',black,'MarkerFaceColor',orange,'MarkerSize',12,'LineStyle','none','LineWidth',2);
        lineProps.col= {gray}; lineProps.width = 3; % black
        A1 = mseb(expTime,mean_sim_HA_FS,std_sim_HA_FS,lineProps,1);
        lineProps.col= {dark_gray}; lineProps.width = 3; % gray
        A2 = mseb(expTime,mean_sim_HA_NonFS,std_sim_HA_NonFS,lineProps,1);
        box on
        ylim([-0.05 1.4]);
        xlim([0 3480]);
        grid on
        set(gca,'linewidth',2.5)
        set(gca,'XGrid','off')
        set(gca,'YGrid','on')
        if p==1
            yticks([0 0.25 0.50 0.75 1.00 1.25 ]);
            set(gca,'Yticklabel',[])
            nameplot = horzcat(pre_namePlot,'_HA');
            set (gca ,'FontSize',18);  set(gca, 'FontName', 'Arial');
            xlabel('Time (sec)','FontSize',22);
        else
            nameplot = horzcat(pre_namePlot,'_SI_HA');
            yticks([0 0.25 0.50 0.75 1.00 1.25 ]);
            yticklabels({'0.00','0.25','0.50', '0.75', '1.0','1.25' });
            set (gca ,'FontSize',18);  set(gca, 'FontName', 'Arial');
            ylabel('Total int. (a.u.)','FontSize',20);
            xlabel('Time (sec)','FontSize',20);
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
    
end

%% Comparing Simulations and Experiemntal data
%dataPoints = 0:120:3480; dataPoints(1)=1;
dataPoints =30;
%fit_HA_FS = sum ( abs(data_HA_FS(1:30) - mean_sim_HA_FS(1:30))) / length(dataPoints);
%fit_HA_NonFs = sum(abs(data_0F_NonFS(1:30) - mean_sim_HA_NonFS(1:30))) / length(dataPoints);
%fit_Flag = sum(abs(data_Flag(1:30) - mean_sim_Flag(1:30))) / length(dataPoints);
%fit_1F = sum ( abs(data_1F(1:30) - mean_sim_1F_FS(1:30))) / length(dataPoints);

fit_HA_FS = sqrt(sum((data_HA_FS(1:dataPoints) - mean_sim_HA_FS(1:dataPoints)).^2)/ dataPoints);
fit_HA_NonFs = sqrt(sum((data_0F_NonFS(1:dataPoints) - mean_sim_HA_NonFS(1:dataPoints)).^2)/ dataPoints);
fit_Flag = sqrt(sum((data_Flag(1:dataPoints) - mean_sim_Flag(1:dataPoints)).^2)/ dataPoints);
fit_1F = sqrt(sum((data_1F(1:dataPoints) - mean_sim_1F_FS(1:dataPoints)).^2)/ dataPoints);


end
