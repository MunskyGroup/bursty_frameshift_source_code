function [fit_HT_Org, fit_HT_HA ] = harringtonine_assays(folderName,nonConsiderTime,nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,independent_Repetitions,elongationFast)
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

fit_HT_Org = (fit_HT_0F+fit_HT_1F)/2;

if fit_HT_Org < 0.19
    condition_Harringtonine = 1;
else
    condition_Harringtonine = 0;
end

%% %%%%%%%%%%%%% HA TAG sequence %%%%%%%%%%%%%%%%%%%%%%%%%
if condition_Harringtonine ==1 || plottingCondition ==1 || runningParameterScan ==  1
    totalSimulationTime = 3481;
    fileName_SunTag_inFS = 'FS_smHA_SunTag.xls'; % SunTAg in FS spots
    fileName_HA_inFS = 'FS_smHA.xls';        % HA in FS spots
    fileName_HA_inNonFS = 'nFS_smHA.xls'; % HA in NON-FS spots
    fileName_FLAG_Only = 'FS_smHA_0F_only.xls';  % Data FLAG only
    
    [rawData_SunTag_inFS,~,~]=xlsread(fileName_SunTag_inFS);
    [rawData_HA_inFS,~,~]=xlsread(fileName_HA_inFS);
    [rawData_HA_inNonFS,~,~]=xlsread(fileName_HA_inNonFS);
    [rawData_FLAG_Only ,~,~]=xlsread(fileName_FLAG_Only);
    
    geneFileName_0F_HA = '0F_HA.txt';
    geneFileName_1F_HA = '1F_HA.txt';
    position_FS_HA = 368;
    expTime = rawData_SunTag_inFS(:,4)';
    
    % Running the SSA
    for j= 1: independent_Repetitions
        [~,intensityVector_0_HA, intensityVector_1_HA,IntensityVectors_HA_HA, T_array, ~, ~, ~,~,~] = masterFunction(nonConsiderTime, geneFileName_0F_HA,geneFileName_1F_HA, ke, position_FS_HA, ki, nRepetitions, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'user', expTime,elongationFast);
        %         [~,intensityVector_0_HA, intensityVector_1_HA,IntensityVectors_HA_HA, T_array, ~, ~, ~,~,~] = masterFunction(nonConsiderTime, geneFileName_0F_HA,geneFileName_1F_HA, ke, position_FS_HA, ki, nRepetitions, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'user', expTime,elongationFast);
        [sim_HA_FS(j,:),sim_1F_FS(j,:),sim_FLAG(j,:),sim_HA_NonFS(j,:), sim_Err_HA_FS(j,:),sim_Err_1F_FS(j,:),sim_Err_HA_NonFS(j,:),sim_Err_Flag(j,:)] = runOff_Intesities_HA(intensityVector_0_HA, intensityVector_1_HA,IntensityVectors_HA_HA);
    end
    % Plotting
    namePlot = 'HT_';
    [fit_HA_FS,fit_HA_NonFs,fit_Flag] = compareRunOffs_HA(namePlot, folderName, rawData_SunTag_inFS,rawData_HA_inFS,rawData_HA_inNonFS,rawData_FLAG_Only, T_array, sim_HA_FS, sim_1F_FS, sim_FLAG,sim_HA_NonFS, sim_Err_HA_FS, sim_Err_1F_FS, sim_Err_HA_NonFS,sim_Err_Flag, plottingCondition,independent_Repetitions);
else
    fit_HA_FS = 1;
    fit_HA_NonFs = 1;
    fit_Flag = 1;
end
% set NANs to 1. The maximum value for the OF.
fit_HT_0F (isnan(fit_HT_0F))=1;
fit_HT_1F  (isnan(fit_HT_1F))=1;
fit_HA_FS  (isnan(fit_HA_FS))=1;
fit_HA_NonFs (isnan(fit_HA_NonFs))=1;

% sum of all harringtonine data sets.
fit_HT_HA = (fit_HA_FS + fit_HA_NonFs)/2;
end
