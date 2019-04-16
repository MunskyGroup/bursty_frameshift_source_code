clear all; close all; clc;
UsingBurstingModel = 1; % 1 for bursting model; 0 for constitutive model.
if UsingBurstingModel==1
    folderName = horzcat('Results_Optimization_B'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    initialMatrix = [0.0244, 3, 0.0234,0.0139 ,9.6e-5]; % Best parameters for a Bursty model
    folderName_2 = horzcat('Results_Optimization_C'); if exist (folderName_2, 'dir') ~= 7; mkdir(folderName_2); end % second folder for the comparison
else
    folderName = horzcat('Results_Optimization_C'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    initialMatrix = [0.0223    9.8295    0.0800    0.0021 , 0]; %  9.8295
    folderName_2 = horzcat('Results_Optimization_B'); if exist (folderName_2, 'dir') ~= 7; mkdir(folderName_2); end % second folder for the comparison
end
%% Using fast SSA. Hybrid stochastic before the FSS and deterministic after FSS.
elongationFast=0; % 0 For fully stochastic system. 1 for hybrid system.
%% deffining time and repetitions
nonConsiderTime = 10000;
nonConsiderTime_short =2000;
nRepetitions = 10000;
nRepetitions_Harringtonine = 192;
runningParameterScan = 0;
independent_Repetitions = 6;
%% Optimization implementation
k_off =0.0013; %
ki_min = 0.01; ki_max = 0.1;
ke_min = 1; ke_max =20;
k_fss_min = 0.0001;  k_fss_max =0.1;
k_s_fss_min = 0.0001;  k_s_fss_max =0.1;
kon_min = 9.6e-5*0.5;  kon_max =9.6e-5*1.5;

%% deffining parameter ranges
lb = [ ki_min, ke_min, k_fss_min, k_s_fss_min, kon_min];
ub = [ki_max, ke_max, k_fss_max, k_s_fss_max, kon_max];
nvars= length(lb);

%% Parameter matrix
fitFnc = @(x)fitnessFunction(x, nonConsiderTime,nonConsiderTime_short, nRepetitions, nRepetitions_Harringtonine,UsingBurstingModel,runningParameterScan,independent_Repetitions,k_off,elongationFast);

%% Running the Pattern Search Algorithm
options = optimoptions('patternsearch','MaxIterations',1000,'MaxTime', 3600*10);
x = patternsearch(fitFnc,initialMatrix,[],[],[],[],lb,ub,options);

%% Reporting Results
plottingCondition=0; runningParameterScan=1;
[fitValue, fit_meanIntensity, fit_ratio, fit_HT_Org, fit_HT_HA,fit_Percent] = reportingResults_Optimization( x, nonConsiderTime,nonConsiderTime_short, nRepetitions, nRepetitions_Harringtonine ,folderName, UsingBurstingModel,k_off,runningParameterScan,independent_Repetitions,elongationFast,plottingCondition);

%% Save best values
cd (folderName)
S = dir('sel_param*.*');
try
    numberOf_Optimization_Repetitions =  length (S);
catch
    numberOf_Optimization_Repetitions = 0;
end
fname=[pwd,'/sel_param_',num2str(numberOf_Optimization_Repetitions+1),'.mat'];
fname2=[pwd,'/sel_fval_',num2str(numberOf_Optimization_Repetitions+1),'.mat'];
save(fname,'x')
save(fname2,'fitValue','fit_meanIntensity','fit_ratio','fit_HT_Org','fit_HT_HA','fit_Percent');
cd ..

%% Section that generates Bar plots alwas takes all information and plots the result with the maximum number of repetitions
folderName_Cell ={'Results_Optimization_B','Results_Optimization_C'};
for k=1:2 % loading data for bursting and constitutive models
    cd (folderName_Cell{k})
    % load for all existing repetitions of the Optimization
    S_load = dir('sel_fval_*.*');
    numberOf_Optimization_Repetitions(k) =  length (S_load);
    if numberOf_Optimization_Repetitions(k)>0
        if k==1
            vector_fitValue = nan(2,numberOf_Optimization_Repetitions(k));
            vector_fit_Percent = nan(2,numberOf_Optimization_Repetitions(k));
            vector_fit_meanIntensity = nan(2,numberOf_Optimization_Repetitions(k));
            vector_fit_ratio = nan(2,numberOf_Optimization_Repetitions(k));
            vector_fit_HT = nan(2,numberOf_Optimization_Repetitions(k));
            vector_fit_HT_HA = nan(2,numberOf_Optimization_Repetitions(k));
        end
        
        for i = 1: numberOf_Optimization_Repetitions(k)
            fileNames =  ['sel_fval_',num2str(i),'.mat'];
            load (fileNames)
            vector_fitValue(k,i) = fitValue;
            vector_fit_Percent(k,i) = fit_Percent;
            vector_fit_meanIntensity(k,i) = fit_meanIntensity;
            vector_fit_ratio(k,i) = fit_ratio;
            vector_fit_HT(k,i) = fit_HT_Org;
            vector_fit_HT_HA(k,i) = fit_HT_HA;
        end
        % Calculating mean
        % matrix mean_fitValue and std_fitValue rows for bursing and constitutive models.
        % columsn represent the Objective Function total, no rib, percentage, HT.
        mean_fitValue(k,1) = nanmean(vector_fitValue(k,:));
        mean_fitValue (k,2)= nanmean(vector_fit_Percent(k,:));
        mean_fitValue(k,3) = nanmean(vector_fit_meanIntensity(k,:));
        mean_fitValue(k,4) = nanmean(vector_fit_ratio(k,:));
        mean_fitValue (k,5)= nanmean(vector_fit_HT(k,:));
        mean_fitValue (k,6)= nanmean(vector_fit_HT_HA(k,:));
        
        % Calculating std
        std_fitValue (k,1)= nanstd(vector_fitValue(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        std_fitValue(k,2) = nanstd(vector_fit_Percent(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        std_fitValue(k,3) = nanstd(vector_fit_meanIntensity(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        std_fitValue(k,4) = nanstd(vector_fit_ratio(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        std_fitValue(k,5) = nanstd(vector_fit_HT(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        std_fitValue(k,6) = nanstd(vector_fit_HT_HA(k,:))./sqrt(numberOf_Optimization_Repetitions(k));
        
    end
    cd ..
end
% check significance
if min(numberOf_Optimization_Repetitions)>1
    sig_fitValue = zeros (1,6);
    pval = zeros (1,6);
    [pval(1),sig_fitValue(1)] = ranksum(vector_fitValue(1,:),vector_fitValue(2,:));
    [pval(2),sig_fitValue(2)] = ranksum(vector_fit_Percent(1,:), vector_fit_Percent(2,:));
    [pval(3),sig_fitValue(3)] = ranksum(vector_fit_meanIntensity(1,:), vector_fit_meanIntensity(2,:));
    [pval(4),sig_fitValue(4)] = ranksum(vector_fit_ratio(1,:), vector_fit_ratio(2,:));
    [pval(5),sig_fitValue(5)] = ranksum(vector_fit_HT(1,:), vector_fit_HT(2,:));
    [pval(6),sig_fitValue(6)] = ranksum(vector_fit_HT_HA(1,:), vector_fit_HT_HA(2,:));
    % Function that plots the data
    plotComparionOptimization('Results_Optimization_B',mean_fitValue,std_fitValue,sig_fitValue,pval)
end