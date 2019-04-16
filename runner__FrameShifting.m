close all; clc; clear all;
%% This code is intended to simulate frameshifting dynamics.

UsingBurstingModel = 1; % 1 for bursting model; 0 for constitutive model.
if UsingBurstingModel==1
    folderName = horzcat('Results_B'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    % Parameters for the bursting model
    ki  = 0.0244 ; % Initiation rate
    ke = 3;  % % Elongation rate.
    k_fss = 0.0234; % Pause in the 0 frame  %0.03
    k_s_fss =  0.0139; % Pause in the -1 frame
    k_on = 9.6e-5; % Swithing to -1 frame. 30 min
    k_off = 0.0013; % Swithcing to 0 frame.
else
    folderName = horzcat('Results_C'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    % Parameters for the constitutive model
    k_on = 0; % Switching rates are not used on the constitutive model
    k_off = 0; % Switching rates are not used on the constitutive model
    ki  = 0.0223 ; % Initiation rate
    ke = 15.83 ;  % % Elongation rate.
    k_fss = 0.078; % Pause in the 0 frame
    k_s_fss =  0.0021; % Pause in the -1 frame
end
%% Using fast SSA. Hybrid stochastic before the FSS and deterministic after FSS.
elongationFast=0; % 0 For fully stochastic system. 1 dor hybrid system.
%% 1.- Deffining simulation times and number of repetitions
nonConsiderTime = 10000;
nonConsiderTime_short = 10000;
totalSimulationTime = 5000;
nRepetitions_TimeCourses = 10000;
nRepetitions_Harringtonine = 1000;
noIndependent_Repetitions_Harringtonine = 8;
position_FS =25;
position_FS_HA = 368;
%% 2.- Loading the experimental data
geneFileName_0F = 'Frame_0.txt'; % (0F) FSS MF tag
geneFileName_1F = 'Frame_1.txt'; % (-1F) FSS MF tag
geneFileName_0F_HT = '2X_FS_0F.txt'; % (0F) FSS 2xMF tag
geneFileName_1F_HT = '2X_FS_1F.txt'; % (-1F) FSS 2xMF tag
geneFileName_0F_HA = '0F_HA.txt'; % (0F) HA MF tag
geneFileName_1F_HA = '1F_HA.txt'; % (-1F) HA MF tag

%% 3.- Running stochastic simulations for the original sequence
timePerturbationApplication = 0; evaluatingFRAP = 0; evaluatingInhibitor = 0; plottingCondition =1; runningParameterScan=0;
[ribosomesPerRNA_BeforeFSS,intensityVector_0, intensityVector_1,IntensityVectors_HA, T_array, geneLength_0F,geneLength_1F, pre_ribosomeDensity,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = masterFunction(nonConsiderTime, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions_TimeCourses, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'finonly',[],0);
% Calculating percentage of spots per frame
[fit_Percent,conditionPercentage] = percentagePerFrame(nRepetitions_TimeCourses, folderName, intensityVector_0, intensityVector_1,plottingCondition);
% Calculating ratio and intensities
[ fit_meanIntensity,fit_ratio,conditionFractions] = compare_intensities_and_fractions (intensityVector_0,intensityVector_1,plottingCondition,folderName);
% Calculating Ribosome Distributions.
ribosomeOccupancy(pre_ribosomeDensity, geneLength_0F, geneLength_1F, folderName, position_FS, plottingCondition);
% Running simulations for the HA tag
[ribosomesPerRNA_BeforeFSS_HA,intensityVector_0_HA, intensityVector_1_HA,IntensityVectors_HA_HA, T_array, geneLength_0_HA, geneLength_1_HA, pre_ribosomeDensity_HA,pre_numOfRibosomes_0_HA,pre_numOfRibosomes_1_HA] = masterFunction(nonConsiderTime, geneFileName_0F_HA,geneFileName_1F_HA, ke,position_FS_HA, ki, nRepetitions_TimeCourses, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'allpos',[],1);
% Calculating Ribosome Distributions.
ribosomeOccupancy( pre_ribosomeDensity_HA, geneLength_0_HA, geneLength_1_HA, folderName, position_FS_HA, plottingCondition);
%% 6.- Performing the Harringtonine Assays
[fit_HT_Org, fit_HT_HA]= harringtonine_assays(folderName,nonConsiderTime,nonConsiderTime_short, geneFileName_0F_HT,geneFileName_1F_HT, ke,position_FS, ki, nRepetitions_Harringtonine, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,noIndependent_Repetitions_Harringtonine,elongationFast);
