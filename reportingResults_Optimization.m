function [fitValue,fit_meanIntensity,fit_ratio, fit_HT_Org, fit_HT_HA,fit_Percent] = reportingResults_Optimization( x, nonConsiderTime, nonConsiderTime_short, nRepetitions, nRepetitions_Harringtonine,folderName ,UsingBurstingModel,k_off,runningParameterScan,independent_Repetitions,elongationFast,plottingCondition)
%% OPTIMIZED Parameters
ki = x(1);
ke = x(2);
k_fss = x(3);
k_s_fss = x(4);
k_on = x(5);
totalSimulationTime =121;

%% 2.- Loading the experimental data
geneFileName_0F = 'Frame_0.txt';
geneFileName_1F = 'Frame_1.txt';
geneFileName_0F_HT = '2X_FS_0F.txt';
geneFileName_1F_HT = '2X_FS_1F.txt';

%% 3.- Running stochastic simulations for the original sequence
position_FS =25;
timePerturbationApplication = 0; evaluatingFRAP = 0; evaluatingInhibitor = 0; 
[~,intensityVector_0, intensityVector_1,IntensityVectors_HA, T_array, geneLength_0, geneLength_1, pre_ribosomeDensity,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = masterFunction(nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'finonly',[],elongationFast);
% Calculating intensities and ratios
[ fit_meanIntensity,fit_ratio,~] = compare_intensities_and_fractions (intensityVector_0,intensityVector_1,plottingCondition,folderName);
%% 4.- Running the simulation for the Original-sequence
% Calculating percentage of spots per frame
[fit_Percent,~] = percentagePerFrame(nRepetitions, folderName, intensityVector_0, intensityVector_1,plottingCondition);
%% 5.- Running the simulation for the HA-sequence
[fit_HT_Org, fit_HT_HA] = harringtonine_assays(folderName,nonConsiderTime, nonConsiderTime_short,geneFileName_0F_HT,geneFileName_1F_HT, ke,position_FS, ki, nRepetitions_Harringtonine, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,independent_Repetitions,elongationFast);
%% Final fit value
fitValue = fit_meanIntensity+fit_ratio+ fit_HT_Org+ fit_HT_HA+fit_Percent;
end
