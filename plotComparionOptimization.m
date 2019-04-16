function plotComparionOptimization(folderName,mean_fitValue,std_fitValue,sig_fitValue,pval)
%% Plot all data
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0,9.5, 3];
num = 6; %number of different subcategories
c = 1:num;
axes1 = axes;
ylabel('\it J(\theta)','Interpreter','tex','FontSize',22);
hold on
box on
set(gca,'linewidth',2)
% Bar(s)
barColor1 = [ 0.1 0.1 0.1];     %back
barColor2 = [1 1 1];     %gray
for i = 1:num
    b1= bar(c(i)-0.2,mean_fitValue(1,i),0.3,'FaceColor',barColor1);
    b2= bar(c(i)+0.2,mean_fitValue(2,i),0.3,'FaceColor',barColor2);
    b1.LineWidth = 1;
b2.LineWidth = 1;
end

errH1 = errorbar(c-0.2,mean_fitValue(1,:),std_fitValue(1,:),'.','Color','k');
errH2 = errorbar(c+0.2,mean_fitValue(2,:),std_fitValue(2,:),'.','Color','k');
errH1.LineWidth = 1.2;
errH2.LineWidth = 1.2;
errH1.Color = [0.3,0.3,0.3];
errH2.Color = [0 0 0];
max_Y = max([mean_fitValue(1,1)+ std_fitValue(1,1), mean_fitValue(2,1)+ std_fitValue(2,1)]);
%% ploting the significance asterisk
for j =1:num
    
    if pval(j) > 0.05
        plot([c(j)-0.2, c(j)+0.2], [1 1]*max_Y*1.1, '-k')
        text(c(j)-0.1,max_Y*1.15, 'ns','FontSize',12)
        
    elseif pval(j) < 0.05 && pval(j) > 0.01
        plot([c(j)-0.2, c(j)+0.2], [1 1]*max_Y*1.1, '-k')
        text(c(j)-0.1,max_Y*1.15, '*','FontSize',12)
        
    elseif pval(j) <= 0.01 && pval(j) > 0.001
        plot([c(j)-0.2, c(j)+0.2], [1 1]*max_Y*1.1, '-k')
        text(c(j)-0.1,max_Y*1.15, '**','FontSize',12)
        
    elseif pval(j) <= 0.001
        plot([c(j)-0.2, c(j)+0.2], [1 1]*max_Y*1.1, '-k')
        text(c(j)-0.1,max_Y*1.15, '***','FontSize',12)
    end
end
box on
set(gca,'linewidth',2)
ylim([0 max_Y*1.25]);
set(axes1,'Xlim',[c(1)-0.5 c(end)+0.5]);
set(axes1,'XTick',[c],'XTickLabel',...
    {'Total','Frac','Int','Ratio','HT','HT_{HA}'});
lgd=legend([b1,b2],'Bursting','Constitutive');
set(lgd,'FontSize',14);
lgd.Location='northeastoutside';
NumTicks=3;
yT = get(gca,'YLim');
set(gca,'YTick',round(linspace(0,max_Y,NumTicks),2))
set (gca ,'FontSize',18); set(gca, 'FontName', 'Arial');
print('-dpng','-r300','plots_comparingOptimizations')
movefile(horzcat('plots_comparingOptimizations', '.png'),horzcat(folderName),'f');
close
end