function [] = make_uncertainty_plot(data,varargin)
%%%%%%%%%%%%%%%%%%%
% Make a plot of uncertainties for a given covariance
% data is a matrix with n_samples rows and n_parameters columns
% log_trans is a boolean; if true, take the log of the data, otherwise false.
% scatter_plot is a boolean; if true, add a scatter plot
% VARARGIN:
%
%%%%%%%%%%%%%%%%%%%
close all

options = struct('nameplot',{{}},'log_trans',0,'scatter_plot',1,'contour_plot',0,'ellipse',1,'parameter_names',{{}},'crlb',[],'kde',0,'true_parameters',[]);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %
    inpName = lower(pair{1}); %# make case insensitive
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end


% add parameter names
npars = size(data,2);
if length(options.parameter_names)==0
    for i=1:npars
        options.parameter_names{i} = ['\lambda_' num2str(i)];
    end
end

if options.log_trans
    data = log(data);
    options.true_parameters = log(options.true_parameters);
end

mu = mean(data);
CV = cov(data);
c = parula;

%% Plotting
fig1 = figure('visible', 'on');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 8, 8];
for i = 1:npars
    for j=1:i
        subplot(npars,npars,(i-1)*npars+j)
        if i==j
            if options.kde
                % do some kernel density estimation
                [f,xi] = ksdensity(data(:,i));
                %                 fill([xi fliplr(xi)],[f zeros(1,length(f))],c(32,:));
                f= f./sum(f);
                fill([xi fliplr(xi)],[f zeros(1,length(f))],[0.5,0.5,0.5]);
                
                if length(options.true_parameters)>0
                    hold on
                    plot([options.true_parameters(i) options.true_parameters(i)],[0,max(f)],'w','LineWidth',2)
                    plot([options.true_parameters(i) options.true_parameters(i)],[0,max(f)],'r','LineWidth',1.5)
                    
                    
                    NumTicks = 2;
                xT = get(gca,'XLim');
                set(gca,'XTick',linspace(xT(1),xT(2),NumTicks))
                yT = get(gca,'YLim');
                set(gca,'YTick',linspace(yT(1),yT(2),NumTicks)) % set(gca,'xticklabel',num2str(tix,'%.1f'))
                
                if xT(2)>1;
                    xtickformat('%.1f')
                else
                    xtickformat('%.2f')
                end
                if yT(2)>1;
                    ytickformat('%.1f')
                else
                    ytickformat('%.2f')
                end
                    
                    
                end
            else
                % or a histogram as default.
                histogram(data(:,i))
                if length(options.true_parameters)>0
                    hold on
                    N = histcounts(data(:,i));
                    plot([options.true_parameters(i) options.true_parameters(i)],[0,max(N)],'r')
                end
                
                
                
                
            end
            
        else
            hold on
            if options.ellipse
                error_ellipse(CV([j i],[j i]),[mu(j),mu(i)]);
                if length(options.crlb>0)
                    hold on
                    error_ellipse(options.crlb([j i],[j i]),[mu(j),mu(i)],'conf',0.9);
                end
                
            end
            if options.scatter_plot
                hold on
                scatter(data(:,j),data(:,i),1,'k','filled')
                alpha(0.5)
            end
            if options.contour_plot
                hold on
                [counts,C] = hist3(data(:,[j i]),[10,10]);
                contour(C{1},C{2},counts')
                %fix axis ticks
                NumTicks = 2;
                xT = get(gca,'XLim');
                set(gca,'XTick',linspace(xT(1),xT(2),NumTicks))
                yT = get(gca,'YLim');
                set(gca,'YTick',linspace(yT(1),yT(2),NumTicks)) % set(gca,'xticklabel',num2str(tix,'%.1f'))
                
                if xT(2)>1;
                    xtickformat('%.1f')
                else
                    xtickformat('%.2f')
                end
                if yT(2)>1;
                    ytickformat('%.1f')
                else
                    ytickformat('%.2f')
                end
                
            end
            if length(options.true_parameters)>0
                hold on
                % add a marker
                scatter(options.true_parameters(j),options.true_parameters(i),'filled','r')
                
                NumTicks = 2;
                xT = get(gca,'XLim');
                set(gca,'XTick',linspace(xT(1),xT(2),NumTicks))
                yT = get(gca,'YLim');
                set(gca,'YTick',linspace(yT(1),yT(2),NumTicks)) % set(gca,'xticklabel',num2str(tix,'%.1f'))
                
                if xT(2)>1;
                    xtickformat('%.1f')
                else
                    xtickformat('%.2f')
                end
                if yT(2)>1;
                    ytickformat('%.1f')
                else
                    ytickformat('%.2f')
                end
                
            end
        end
        if i==npars
            xlabel(options.parameter_names{j})
        end
        if j==1
            ylabel(options.parameter_names{i})
        end
    end
end

% set (gca ,'FontSize',12,'FontName', 'Arial');

% % add parameter names
% npars = size(data,2);
% if length(options.parameter_names)==0
%     for i=1:npars
%         options.parameter_names{i} = ['\lambda_' num2str(i)];
%     end
% end
%
% if options.log_trans
%     data = log10(data);
% end
%
% mu = mean(data);
% CV = cov(data);
% c = parula;
%
% %% Plotting
% fig1 = figure('visible', 'on');
% fig1= gcf;
% fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0, 0, 8, 6];
%
%
% for i = 1:npars
%     for j=1:i
%         subplot(npars,npars,(i-1)*npars+j)
%         if i==j
%             if options.kde
%                 % do some kernel density estimation
%                 [f,xi] = ksdensity(data(:,i));
%                 fill([xi fliplr(xi)],[f zeros(1,length(f))],c(32,:));
%             else
%                 % or a histogram as default.
%                 histogram(data(:,i))
%             end
%
%         else
%             hold on
%             if options.ellipse
%                 error_ellipse(CV([j i],[j i]),[mu(j),mu(i)]);
%                 if length(options.crlb>0)
%                     hold on
%                     error_ellipse(options.crlb([j i],[j i]),[mu(j),mu(i)]);
%                 end
%
%             end
%             if options.scatter_plot
%                 scatter(data(:,j),data(:,i),1,'k','filled')
%                 alpha(0.5)
%             end
%             if options.contour_plot
%                 [counts,C] = hist3(data(:,[j i]),[10,10]);
%                 contour(C{1},C{2},counts')
%             end
%         end
%        if i==npars
%           xlabel(options.parameter_names{j})
%        end
%        if j==1
%            ylabel(options.parameter_names{i})
%        end
%     end
% end
print('-dpng','-r600',options.nameplot)