%% fitnes function
function fitValue = fitnessFunction (x, nonConsiderTime,nonConsiderTime_short, nRepetitions, nRepetitions_Harringtonine,UsingBurstingModel,runningParameterScan,independent_Repetitions,k_off,elongationFast)
%% PARAMETERS TO OPTIMIZE
ki = x(1);
ke = x(2);
k_fss = x(3);
k_s_fss = x(4);
k_on = x(5);
%% Get Gene info
geneFileName_0F_HA = '2X_FS_0F.txt';
geneFileName_1F_HA = '2X_FS_1F.txt';
geneFileName_0F = 'Frame_0.txt';
geneFileName_1F = 'Frame_1.txt';
totalSimulationTime = 121;
% parameters for the FS
position_FS = 25;
timePerturbationApplication = 0; evaluatingFRAP = 0; evaluatingInhibitor = 0; plottingCondition =0;
%% Running stochastic simulations
[~,intensityVector_0, intensityVector_1,~, ~, ~, ~, ~,~,~] = masterFunction(nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'finonly',[],elongationFast);
%% Plotting Section
% Calculating percentage of spots per frame
[fit_Percent,conditionPercentage] = percentagePerFrame(nRepetitions, '', intensityVector_0, intensityVector_1,plottingCondition);
% Calculating intensities and ratios
[ fit_meanIntensity,fit_ratio,conditionFractions] = compare_intensities_and_fractions (intensityVector_0,intensityVector_1,plottingCondition,'');
% Performing the Harringtonine Assays
if (conditionFractions ==1 && conditionPercentage==1) || runningParameterScan ==1
    [fit_HT_Org, fit_HT_HA]= harringtonine_assays('',nonConsiderTime,nonConsiderTime_short, geneFileName_0F_HA,geneFileName_1F_HA, ke,position_FS, ki, nRepetitions_Harringtonine, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,independent_Repetitions,elongationFast);
else
    fit_HT_Org =1;
    fit_HT_HA =1;
end
%% Final fit value
fitValue = fit_meanIntensity + fit_ratio+ fit_HT_Org + fit_HT_HA + fit_Percent;
fitValue(isnan(fitValue)==1)=inf;

end