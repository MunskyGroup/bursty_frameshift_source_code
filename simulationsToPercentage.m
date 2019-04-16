function [sim_0F,sim_1F,sim_BF,err_sim_0F,err_sim_1F,err_sim_BF]  = simulationsToPercentage (numberOfSpots,StrData,observationTime)

%% Bootstraping to calculate error in number of spots per frame
parfor j =1: 1000
    numberOfTrajectoriesInZeroFrame = 0;
    numberOfTrajectoriesInMinusOneFrame = 0;
    numberOfTrajectoriesInZeroAndMinusOneFrame = 0;
    Sampling_numberOfSpots =1000; % experimental number of spots
    for i =1:Sampling_numberOfSpots
        ii= randi(numberOfSpots);
        if max(StrData(1,ii).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,ii).trajectory_2(end-observationTime:end)) == 0
            numberOfTrajectoriesInZeroFrame = numberOfTrajectoriesInZeroFrame + 1;
        end
        if max(StrData(1,ii).trajectory_2(end-observationTime:end)) > 0  && max(StrData(1,ii).trajectory_1(end-observationTime:end)) == 0
            numberOfTrajectoriesInMinusOneFrame = numberOfTrajectoriesInMinusOneFrame + 1;
        end
        if max(StrData(1,ii).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,ii).trajectory_2(end-observationTime:end)) > 0
            numberOfTrajectoriesInZeroAndMinusOneFrame = numberOfTrajectoriesInZeroAndMinusOneFrame + 1;
        end
    end
    numberOfTrajectoriesInZeroFrame = numberOfTrajectoriesInZeroFrame + 1/2*(numberOfTrajectoriesInZeroAndMinusOneFrame);
    numberOfTrajectoriesInMinusOneFrame = numberOfTrajectoriesInMinusOneFrame + 1/2*(numberOfTrajectoriesInZeroAndMinusOneFrame);
    vec_sim_0F(j) = (numberOfTrajectoriesInZeroFrame./Sampling_numberOfSpots)*100;
    vec_sim_1F(j) = (numberOfTrajectoriesInMinusOneFrame./Sampling_numberOfSpots)*100;
    vec_sim_BF(j) = (numberOfTrajectoriesInZeroAndMinusOneFrame./Sampling_numberOfSpots)*100;
end

%% Error Propagation. 1/2 of spots in Both frames is added to the -1 and 0 spots.
err_sim_0F = std(vec_sim_0F);
err_sim_1F = std(vec_sim_1F);
err_sim_BF = std(vec_sim_BF);


% err_sim_0F_withPropagation =sqrt (err_sim_0F.^2 + err_sim_1F.^2);
% err_sim_1F_withPropagation =sqrt (err_sim_0F.^2 + err_sim_1F.^2);


%% Calculating the mean number of spots per frame
numberOfTrajectoriesInZeroFrame_m = 0;
numberOfTrajectoriesInMinusOneFrame_m = 0;
numberOfTrajectoriesInZeroAndMinusOneFrame_m = 0;
numberOfTrajectoriesNonTranslating_m = 0;
for i =1:numberOfSpots
    if max(StrData(1,i).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,i).trajectory_2(end-observationTime:end)) == 0
        numberOfTrajectoriesInZeroFrame_m = numberOfTrajectoriesInZeroFrame_m + 1;
    end
    if max(StrData(1,i).trajectory_2(end-observationTime:end)) > 0  && max(StrData(1,i).trajectory_1(end-observationTime:end)) == 0
        numberOfTrajectoriesInMinusOneFrame_m = numberOfTrajectoriesInMinusOneFrame_m + 1;
    end
    if max(StrData(1,i).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,i).trajectory_2(end-observationTime:end)) > 0
        numberOfTrajectoriesInZeroAndMinusOneFrame_m = numberOfTrajectoriesInZeroAndMinusOneFrame_m + 1;
    end
    if max(StrData(1,i).trajectory_1(end-observationTime:end)) == 0 && max(StrData(1,i).trajectory_2(end-observationTime:end)) == 0
        numberOfTrajectoriesNonTranslating_m = numberOfTrajectoriesNonTranslating_m + 1;
    end
end
numberOfSpots = numberOfSpots - numberOfTrajectoriesNonTranslating_m;
sim_0F = (numberOfTrajectoriesInZeroFrame_m./numberOfSpots)*100;
sim_1F = (numberOfTrajectoriesInMinusOneFrame_m./numberOfSpots)*100;
sim_BF = (numberOfTrajectoriesInZeroAndMinusOneFrame_m./numberOfSpots)*100;

end