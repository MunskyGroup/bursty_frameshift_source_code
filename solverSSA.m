function [t_out,ribosomesPerRNA_BeforeFSS, pre_ribosomeDensity,IntensityVectors_0F,IntensityVectors_1F,IntensityVectors_HA,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = solverSSA( nRepetitions, UsingBurstingModel, totalSimulationTime, k_elongation_0F, k_elongation_1F, gene_total, position_FS, k_on, k_off, k_fss, k_s_fss, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, nonConsiderTime,probePosition_0F,probePosition_1F,probePosition_HA,itims,def_times,elongationFast)
switch itims
    case 'all'
        TimeVectorFixedSize = linspace (0,totalSimulationTime,totalSimulationTime);
    case 'allpos'
        TimeVectorFixedSize = [0,linspace(nonConsiderTime,totalSimulationTime,totalSimulationTime-nonConsiderTime+1)];
    case 'sparsepos'
        NT = 50;
        TimeVectorFixedSize = [0,linspace(nonConsiderTime,totalSimulationTime,NT+1)];
    case 'finonly'
        TimeVectorFixedSize = [0,nonConsiderTime,totalSimulationTime];
    case 'user'
        TimeVectorFixedSize = [0,nonConsiderTime,nonConsiderTime+def_times];
end
nonConsiderTime = 2;
% We can save a little time by ignoring the early time points.
%% Running the SSA
maximum_Number_Ribosomes = 100;
if position_FS ==25
    X_SteadyState = zeros(maximum_Number_Ribosomes,1);
else
    X_SteadyState = zeros(maximum_Number_Ribosomes,1);
end

parfor k = 1: nRepetitions
    try
        [RibsomePositions{k}] = SSA_BEM_mex(UsingBurstingModel, position_FS, gene_total,maximum_Number_Ribosomes, X_SteadyState, k_elongation_0F, k_elongation_1F, TimeVectorFixedSize, k_on, k_off, k_fss, k_s_fss,timePerturbationApplication, evaluatingInhibitor, nonConsiderTime,elongationFast);
    catch
        [RibsomePositions{k}] = SSA_BEM(UsingBurstingModel, position_FS, gene_total,maximum_Number_Ribosomes, X_SteadyState, k_elongation_0F, k_elongation_1F, TimeVectorFixedSize, k_on, k_off, k_fss, k_s_fss,timePerturbationApplication, evaluatingInhibitor, nonConsiderTime,elongationFast);
    end
end

%% Saving ribosome position in time
parfor k = 1: nRepetitions
    [ribosomesPerRNA_BeforeFSS{k},pre_ribosomeDensity{k}, pre_numOfRibosomes_0{k},pre_numOfRibosomes_1{k},IntensityVectors_0F(k,:),IntensityVectors_1F(k,:),IntensityVectors_HA(k,:)] = SSAtoIntensity(nonConsiderTime,TimeVectorFixedSize,gene_total,RibsomePositions{k},probePosition_0F,probePosition_1F,probePosition_HA,k_elongation_0F,position_FS );
end
t_out = TimeVectorFixedSize(3:end)-TimeVectorFixedSize(2);
