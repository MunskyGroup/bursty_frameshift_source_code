function ribosomeOccupancy(pre_ribosomeDensity_cell, geneLength_1,geneLength_2,folderName,position_FS,plottingCondition )
% Function to calculate the ribosome occupancy
nSpots = size(pre_ribosomeDensity_cell,2);
geneLength = size(pre_ribosomeDensity_cell{1,1},2);
pre_ribosomeDensity = zeros(nSpots,geneLength);

for i =1:nSpots
    pre_ribosomeDensity (i,:)= pre_ribosomeDensity_cell{1,i};
end

ribosomeFootprint =9;
ribosomeDensity = mean(pre_ribosomeDensity);
ribosomeDensity_0 = ribosomeDensity(1:geneLength_1);
ribosomeDensity_1 = [zeros(1,position_FS),ribosomeDensity([geneLength_1+1: end])];
minSizeDensities = min( length(ribosomeDensity_0), length(ribosomeDensity_1));
ribosomeDensity_BF = ribosomeDensity_0(1:minSizeDensities) + ribosomeDensity_1(1:minSizeDensities);
ribosomeDensity_Hist_BF = ribosomeDensity_BF./sum(ribosomeDensity_BF);
bins_geneSize = linspace(1,length(ribosomeDensity_Hist_BF),length(ribosomeDensity_Hist_BF));
smothed_ribosomeDensity = smooth(ribosomeDensity_Hist_BF,ribosomeFootprint)';

blue = [0 0.6 1];
green =[0.2 1 0.2];
cyan = [0 1 1];
red = [1 0 0];
gray = [0.5 0.5 0.5];
orange = [255, 170, 85]./255;
black =[0 0 0];

probePosition_0F = [54    82   110   138   166   194   222   250   278   306   334   362   400   428   456   484   512   540   568   596   624   652   680   708];
probePosition_1F = [62    90   118   146   174   202   230   258   286   314   342   370   408   436   464   492   520   548   576   604   632   660   688   716];
probePosition_0F_HA = [ 427   455   483   511   539   567   595   623   651   679   707];
probePosition_1F_HA = [407   435   463   491   519   547   575   603   631   659   687   715];
probePosition_HA_HA = [2    12    22   207   218   231   242   321   331   341];

if plottingCondition ==1
    if position_FS ==25
        geneType ='OR';
    else
        geneType ='HA';
    end
    
    %% Plot Ribosome Occupancy
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 5, 3.13];
    hold on
    if  position_FS ==25
        for i = 2:length(bins_geneSize)
            h1=bar(i,smothed_ribosomeDensity(i));
            if (bins_geneSize(i) == position_FS)==1 || ( bins_geneSize(i)-1 == position_FS )==1
                set(h1,'FaceColor',red, 'EdgeColor',red);
            elseif any(bins_geneSize(i) == probePosition_0F )==1 || any(bins_geneSize(i)+2 == probePosition_0F )==1
                set(h1,'FaceColor',green, 'EdgeColor',green);
            elseif any(bins_geneSize(i) == probePosition_1F)==1 || any(bins_geneSize(i)+2 == probePosition_1F )==1
                set(h1,'FaceColor',blue, 'EdgeColor',blue);
            else
                set(h1,'FaceColor',gray, 'EdgeColor','none');
            end
        end
        hold off
    else
        for i = 2:length(bins_geneSize)
            h1=bar(i,smothed_ribosomeDensity(i));
            if (bins_geneSize(i) == position_FS ) ==1 || ( bins_geneSize(i)-1 == position_FS )==1
                set(h1,'FaceColor',red, 'EdgeColor',red);
            elseif any(bins_geneSize(i) == probePosition_0F_HA) ==1 || any(bins_geneSize(i)+2 == probePosition_0F_HA )==1
                set(h1,'FaceColor',green, 'EdgeColor',green);
            elseif any(bins_geneSize(i) == probePosition_1F_HA) ==1 || any(bins_geneSize(i)+2 == probePosition_1F_HA )==1
                set(h1,'FaceColor',blue, 'EdgeColor',blue);
            elseif any(bins_geneSize(i) == probePosition_HA_HA) ==1 || any(bins_geneSize(i)+2 == probePosition_HA_HA )==1
                set(h1,'FaceColor',orange, 'EdgeColor',orange);
            else
                set(h1,'FaceColor',gray, 'EdgeColor','none');
            end
        end
        hold off
    end
    xlim([-10 geneLength_1])
    ylim([0 8e-3])
    xticks([200 400 600 800 1000 1200 1400])
    box on
    %axis tight
    set(gca,'linewidth',2)
    set(gca,'fontsize',18 , 'FontName', 'Arial')
    xlabel('Codon Position','FontSize',20);
    ylabel('Probability','FontSize',20);
    nameplot = horzcat('Rib_Dens_',geneType);
    print('-dpng','-r600',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end
end
