function [sim_0F,sim_avg_0F_FS_1F,sim_Err_0F,sim_error_0F_FS_1F] = runOff_Intesities_OR(intensityVector_0, intensityVector_1)
%% function that separates intensity vectors on the following cases:
% % % sim_0F % INTENSITY FOR 0 FRAME NON FRAMESHIFTING SPOTS
% % % sim_1F % INTENSITY FOR -1 FRAME FRAMESHIFTING SPOTS

%% Loading experimental data
fileName_0F = '0F.xls';
fileName_1F = '0and-1F.xls';
[rawData_0F,~,~]=xlsread(fileName_0F);
[rawData_1F,~,~]=xlsread(fileName_1F);

%% Deffining the experimental time points
time_Normalization = 3;

numberOfTimePoints = size(intensityVector_0,2);
numberOfRibosomePositions = size(intensityVector_0,1);

%% Prealocating memory

intensity_0F_ONLY = zeros (numberOfRibosomePositions,numberOfTimePoints);
intensity_0F_FS =  zeros (numberOfRibosomePositions,numberOfTimePoints);
intensity_1F_FS = zeros (numberOfRibosomePositions,numberOfTimePoints);

%% BACKGROUND INTENSITIES ARE DEFFINED AS THE LAST 3 TIME POINTS OF THE EXPERIMENT.
BG_rawData_0F = mean (rawData_0F(end-3:end,2));
BG_rawData_1F = mean (rawData_1F(end-3:end,2));
BG_rawData_0F(BG_rawData_0F<0)=0;
BG_rawData_1F(BG_rawData_1F<0)=0;

%% Intensity classification
for i =1:numberOfRibosomePositions
    FS = sum(intensityVector_1(i,:));
    % 0 frame only
    for tp = 1:size(intensityVector_0,2)
        FS_untilTP = any (intensityVector_1(i,1:tp));
        if FS ==0
            intensity_0F_ONLY(i,tp) = intensityVector_0(i,tp);
        end
        % 0F FS
        if FS ~=0 &&  FS_untilTP ~=0
            intensity_0F_FS(i,tp) = intensityVector_0 (i,tp);
        end
        % -1F frameshifting
        if FS ~=0 &&  FS_untilTP ~=0
            intensity_1F_FS(i,tp) = intensityVector_1 (i,tp);
        end
    end
end

%% removing lines with empty intensities.
rm_zeros_intensity_0F_ONLY = intensity_0F_ONLY(any(intensity_0F_ONLY,2),:);
rm_zeros_intensity_0F_FS = intensity_0F_FS(any(intensity_0F_FS,2),:);
rm_zeros_intensity_1F_FS = intensity_1F_FS(any(intensity_1F_FS,2),:);

%% Calculating total mean intensity

% INTENSITY FOR 0 FRAME in FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_0F_FS)==0
    if size(rm_zeros_intensity_0F_FS,1)>1
        rm_zeros_intensity_0F_FS = rm_zeros_intensity_0F_FS./max(max(rm_zeros_intensity_0F_FS)); % First normalization with respect to the maximum value in all repetitions
        sim_0F_FS = mean (rm_zeros_intensity_0F_FS ); % Calculating mean intensity
        sim_0F_FS = sim_0F_FS./mean(sim_0F_FS(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_0F_FS = sim_0F_FS + BG_rawData_0F; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_0F_FS = sim_0F_FS./mean(sim_0F_FS(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        
    else
        sim_0F_FS = zeros(1, numberOfTimePoints);
    end
else
    sim_0F_FS = zeros(1, numberOfTimePoints);
end

% INTENSITY FOR 0 FRAME in NON FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_0F_ONLY)==0
    if size(rm_zeros_intensity_0F_ONLY,1)>1
        rm_zeros_intensity_0F_ONLY = rm_zeros_intensity_0F_ONLY./max(max(rm_zeros_intensity_0F_ONLY)); % First normalization with respect to the maximum value in all repetitions
        sim_0F = mean (rm_zeros_intensity_0F_ONLY ); % Calculating mean intensity
        sim_0F = sim_0F./mean(sim_0F(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_0F = sim_0F + BG_rawData_0F; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_0F = sim_0F./mean(sim_0F(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_0F = nanstd (rm_zeros_intensity_0F_ONLY(:,:))./sqrt(rawData_0F(:,1)');
    else
        sim_0F = zeros(1, numberOfTimePoints);
        sim_Err_0F = zeros(1, numberOfTimePoints);
    end
else
    sim_0F = zeros(1, numberOfTimePoints);
    sim_Err_0F = zeros(1, numberOfTimePoints);
end

% INTENSITY FOR -1 FRAME in NON-FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_1F_FS)==0
    if size(rm_zeros_intensity_1F_FS,1)>1
        rm_zeros_intensity_1F_FS = rm_zeros_intensity_1F_FS./max(max(rm_zeros_intensity_1F_FS)); % First normalization with respect to the maximum value in all repetitions
        sim_1F = mean (rm_zeros_intensity_1F_FS ); % Calculating mean intensity
        sim_1F = sim_1F + BG_rawData_1F; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_1F = sim_1F./mean(sim_1F(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_1F = nanstd (rm_zeros_intensity_1F_FS(:,:))./sqrt(rawData_1F(:,1)');
    else
        sim_1F = zeros(1, numberOfTimePoints);
        sim_Err_1F = zeros(1, numberOfTimePoints);
    end
else
    sim_1F = zeros(1, numberOfTimePoints);
    sim_Err_1F = zeros(1, numberOfTimePoints);
end

%% Calculating Reported value as the average on -1 frame and 0 frame FS spots
sim_avg_0F_FS_1F = mean([sim_0F_FS;sim_1F ]);

%% Error propagation for the reported error on -1 frame and 0 frame FS spots
sim_error_0F_FS_1F = sqrt (sim_Err_0F.^2 + sim_Err_1F.^2);

end
