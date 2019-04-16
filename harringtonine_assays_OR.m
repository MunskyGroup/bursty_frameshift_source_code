function [fit_HT_Org ] = harringtonine_assays_OR(folderName,nonConsiderTime,nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,independent_Repetitions,elongationFast)
timePerturbationApplication = nonConsiderTime+ 0; evaluatingFRAP = 0; evaluatingInhibitor = 1;
timePerturbationApplication_short = nonConsiderTime_short+ 0;

%% %%%%%%%%%%%%% Original sequence %%%%%%%%%%%%%%%%%%%%%%%%%
totalSimulationTime = 1740;
% loading the experimental data
fileName_0F = '0F.xls';
fileName_1F = '0and-1F.xls';
[rawData_0F,~,~]=xlsread(fileName_0F);
[rawData_1F,~,~]=xlsread(fileName_1F);
expTime = rawData_1F(:,4)';

% Running the SSA
for j= 1: independent_Repetitions
    [~,intensityVector_0, intensityVector_1,~, T_array, ~, ~, ~,~,~] = masterFunction(nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, totalSimulationTime, timePerturbationApplication_short, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss, 'user', expTime,0);
    [sim_0F(j,:),sim_1F(j,:),sim_Err_0F(j,:),sim_Err_1F(j,:)] = runOff_Intesities_OR(intensityVector_0, intensityVector_1);
end
% Plotting
namePlot = 'HT_';
[fit_HT_0F,fit_HT_1F ] = compareRunOffs(namePlot, folderName, rawData_0F, rawData_1F, T_array,sim_0F, sim_1F,sim_Err_0F,sim_Err_1F,plottingCondition,independent_Repetitions);

fit_HT_0F (isnan(fit_HT_0F))=10;
fit_HT_1F  (isnan(fit_HT_1F))=10;
fit_HT_Org = (fit_HT_0F+fit_HT_1F)/2;

end
