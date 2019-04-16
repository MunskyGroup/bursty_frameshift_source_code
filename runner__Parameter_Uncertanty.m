close all; clc; clear all;
UsingBurstingModel = 1; % 1 for bursting model; 0 for constitutive model.
k_off = 0.0013; % Swithcing to 0 frame.

if UsingBurstingModel==1
    folderName = horzcat('Results_ParameterUncertanty_B'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    pars_best = [0.0244, 3, 0.0234,0.0139, 9.6e-5];
    pars_real = pars_best;
else
    folderName = horzcat('Results_ParameterUncertanty_C'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    pars_best = [0.0223 ,   9.8295 ,   0.0800  ,  0.0021, 0];
    pars_real = pars_best;
end
numberTotalEvaluations=1000;

%% Using fast SSA. Hybrid stochastic before the FSS and deterministic after FSS.
elongationFast=0; % 0 For fully stochastic system. 1 dor hybrid system.
finvalmin=inf;
ftvals_best = [10 10 10 10 10];
pars = pars_best;
pars_sv = [];
inorout = [];
J_sv = [];
Ftvals = [1 1 1 1 1];
for ip=1:numberTotalEvaluations
    
    ip
    if ip==1
        delt = 0;
    else
        delt = 0.1;
    end
    %% Deffing parameters.
    %%%%%%%%%%%%%%%%%%
    pars_new = pars.*(1+delt*randn(size(pars)));
    pars_sv = [pars_sv;pars_new];
    
    ki  = pars_new(1); % Initiation rate
    ke = pars_new(2);  % % Elongation rate.
    k_fss =pars_new(3); % Pause in the 0 frame
    k_s_fss = pars_new(4); % Pause in the -1 frame
    k_on = pars_new(5); % Swithing to -1 frame. 30 min
    % k_off =pars_new(6); % Swithcing to 0 frame.
    %% 1.- Deffining simulation times and number of repetitions
    nonConsiderTime = 10000;
    nonConsiderTime_short = 2000;
    totalSimulationTime = 200;
    nRepetitions = 4000;
    totalSimulationTime_short = 121;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nRepetitions_Harringtonine = 192;
    independent_Repetitions = 6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    position_FS =25;
    position_FS_HA = 368;
    
    %% 2.- Loading the experimental data
    %     folderName = horzcat('Results'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end
    geneFileName_0F = 'Frame_0.txt';
    geneFileName_1F = 'Frame_1.txt';
    geneFileName_0F_HT = '2X_FS_0F.txt';
    geneFileName_1F_HT = '2X_FS_1F.txt';
    geneFileName_0F_HA = '0F_HA.txt';
    geneFileName_1F_HA = '1F_HA.txt';
    % deffining additional parameters
    %% 3.- Running stochastic simulations for the original sequence
    timePerturbationApplication = 0; evaluatingFRAP = 0; evaluatingInhibitor = 0; plottingCondition =0; runningParameterScan=1;
    [ribosomesPerRNA_BeforeFSS,intensityVector_0, intensityVector_1,IntensityVectors_HA, T_array, ~, ~, ~,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = masterFunction(nonConsiderTime_short, geneFileName_0F,geneFileName_1F, ke,position_FS, ki, nRepetitions, totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,'finonly',[],elongationFast);
    % Calculating intensities and ratios
    [ fit_meanIntensity,fit_ratio,~] = compare_intensities_and_fractions (intensityVector_0,intensityVector_1,plottingCondition,folderName);
    [fit_Percent,~] = percentagePerFrame(nRepetitions, folderName, intensityVector_0, intensityVector_1,plottingCondition);
    
    %% 6.- Performing the Harringtonine Assays
    if ip==1||(fit_meanIntensity+fit_ratio+fit_Percent + 0.90*min(Ftvals(~isnan(Ftvals(:,4)),4))) < 1.1*sum(ftvals_best(1:4))
        
        [fit_HT_Org, fit_HT_HA]= harringtonine_assays(folderName,nonConsiderTime,nonConsiderTime_short, geneFileName_0F_HT,geneFileName_1F_HT, ke,position_FS, ki, nRepetitions_Harringtonine, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss,plottingCondition,runningParameterScan,independent_Repetitions,elongationFast);
        fit_HT= (fit_HT_Org+ fit_HT_HA)/2;
        finval = fit_meanIntensity + fit_ratio + fit_HT+ fit_Percent;
        temp_finval = [fit_meanIntensity, fit_ratio, fit_Percent, fit_HT];
        Ftvals = [Ftvals;fit_meanIntensity, fit_ratio,fit_Percent,fit_HT, finval];
        
        % Updates best parametes if it finds something better than given best parameters
        % if finval<finvalmin
        if  all (temp_finval < ftvals_best(1:4))==1
            pars_best = pars_new;
            ftvals_best =  [fit_meanIntensity, fit_ratio, fit_Percent,fit_HT, finval];
            finvalmin = finval;
        end
        
        % Classify parameters if they are above a threshold. After 20
        % Random search
        %         if finval<=1.05*finvalmin && ip>20
        if all (temp_finval<=1.1*ftvals_best(1:4))==1 && ip>10
            pars = pars_new;
            inorout = [inorout,1];
            J_sv = [J_sv;Ftvals];
        else
            inorout = [inorout,0];
        end
        
    else
        Ftvals = [Ftvals;fit_meanIntensity, fit_ratio,fit_Percent, NaN, NaN];
        % disp('Skipping this parameter set based upon initial evaluations.')
        inorout = [inorout,0];
    end
    
end
Parameter_AboveThreshold = pars_sv(inorout==1,:);

%% Save best values
cd (folderName)
S = dir('randomSearchData*.*');
try
    numberOf_RS_Repetitions =  length (S);
catch
    numberOf_RS_Repetitions = 0;
end
fname=[pwd,'/randomSearchData_',num2str(numberOf_RS_Repetitions+1),'.mat'];
save(fname,'Parameter_AboveThreshold','pars_best','J_sv')

%% Load previous list of parameters that fullfil the selection criterion
sel_param_list = [0,0,0,0,0];
% load all repetitions to Plot data.
S_load = dir('randomSearchData*.*');
numberOf_RS_Repetitions =  length (S_load);
for i = 1: numberOf_RS_Repetitions
    fileNames =  ['randomSearchData_',num2str(i),'.mat'];
    load (fileNames)
    sel_param_list= [sel_param_list;Parameter_AboveThreshold];
end
% removing rows with only zeros.
sel_param_list = sel_param_list(any(sel_param_list,2),:);
cd ..

%% Plotting the Covariance matrix.
if UsingBurstingModel==1
    namePlot = ['Param_UC_Burst_',num2str(numberOf_RS_Repetitions)];
    make_uncertainty_plot(sel_param_list(:,[1:end]),'nameplot',namePlot,'contour_plot',1,'ellipse',0,'parameter_names',{'k_i','k_e' 'k_{fss}','k_{sfss}','k_{on}'},'kde',1,'true_parameters',pars_real(1:end))
    movefile(horzcat(namePlot, '.png'),horzcat(folderName),'f');
else
    namePlot = ['Param_UC_Const_',num2str(numberOf_RS_Repetitions)];
    make_uncertainty_plot(sel_param_list(:,[1:end-2]),'nameplot',namePlot,'contour_plot',1,'ellipse',0,'parameter_names',{'k_i','k_e' 'k_{fss}','k_{sfss}'},'kde',1,'true_parameters',pars_real(1:end-2))
    movefile(horzcat(namePlot, '.png'),horzcat(folderName),'f');
end