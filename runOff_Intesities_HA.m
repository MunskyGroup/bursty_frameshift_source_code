function [sim_HA_FS,sim_1F_FS,sim_FLAG,sim_HA_NonFS,sim_Err_HA_FS,sim_Err_1F_FS,sim_Err_HA_NonFS,sim_Err_Flag] = runOff_Intesities_HA(intensityVector_0, intensityVector_1,intensityVector_HA)

%% function that separates intensity vectors on the following cases:
% % % sim_HA_FS % INTENSITY FOR 0 FRAME FRAMESHIFTING SPOTS
% % % sim_1F_FS % INTENSITY FOR -1 FRAME FRAMESHIFTING SPOTS
% % % sim_FLAG % INTENSITY FOR 0 FRAME ONLY NO FRAMESHIFTING SPOTS
% % % sim_HA_NonFS % INTENSITY FOR 0 FRAME FRAMESHIFTING AND NON FRAMESHIFTING SPOTS

%% Loading experimental data

fileName_SunTag_inFS = 'FS_smHA_SunTag.xls'; % SunTAg in FS spots
fileName_HA_inFS = 'FS_smHA.xls';        % HA in FS spots
fileName_HA_inNonFS = 'nFS_smHA.xls'; % HA in NON-FS spots
fileName_FLAG_Only = 'FS_smHA_0F_only.xls';  % Data FLAG only

[rawData_SunTag_inFS,~,~]=xlsread(fileName_SunTag_inFS);
[rawData_HA_inFS,~,~]=xlsread(fileName_HA_inFS);
[rawData_HA_inNonFS,~,~]=xlsread(fileName_HA_inNonFS);
[rawData_FLAG_Only ,~,~]=xlsread(fileName_FLAG_Only);

time_Normalization = 3;
numberOfTimePoints = size(intensityVector_0,2);
numberOfRibosomePositions = size(intensityVector_0,1);


%% Prealocating memory
intensity_sim_HA_FS =zeros (numberOfRibosomePositions,numberOfTimePoints);
intensity_sim_1F_FS = zeros (numberOfRibosomePositions,numberOfTimePoints);
intensity_sim_FLAG =  zeros (numberOfRibosomePositions,numberOfTimePoints);
intensity_sim_HA_NonFS = zeros (numberOfRibosomePositions,numberOfTimePoints);

%% BACKGROUND INTENSITIES ARE DEFFINED AS THE LAST 3 TIME POINTS OF THE EXPERIMENT.
BG_rawData_SunTag_inFS= mean (rawData_SunTag_inFS(end-3:end,2));
BG_rawData_HA_inFS = mean (rawData_HA_inFS(end-3:end,2));
BG_rawData_HA_inNonFS = mean (rawData_HA_inNonFS(end-3:end,2));
BG_rawData_FLAG_Only = mean (rawData_FLAG_Only(end-3:end,2));

% Deffining backgraund intensities as positive values.
BG_rawData_SunTag_inFS(BG_rawData_SunTag_inFS<0)=0;
BG_rawData_HA_inFS(BG_rawData_HA_inFS<0)=0;
BG_rawData_HA_inNonFS(BG_rawData_HA_inNonFS<0)=0;
BG_rawData_FLAG_Only(BG_rawData_FLAG_Only<0)=0;

%% Intensity classification
for i =1:size(intensityVector_0,1)
    FS = sum(intensityVector_1(i,:));
    for tp = 1:size(intensityVector_0,2)
        FS_untilTP = any (intensityVector_1(i,1:tp));
        % % sim_HA_FS        % INTENSITY FOR HA in FS SPOTS
        if  FS ~=0 && FS_untilTP ~= 0
            intensity_sim_HA_FS (i,tp) = intensityVector_HA (i,tp);
        end
        % % % sim_1F_FS    % INTENSITY FOR -1 FRAME, FRAMESHIFTING SPOTS & HA
        if FS ~=0 && FS_untilTP ~= 0
            intensity_sim_1F_FS (i,tp) = intensityVector_1 (i,tp);
        end
        % % % sim_FLAG   % INTENSITY FOR 0 FRAME ONLY, NO FRAMESHIFTING SPOTS & HA
        if FS ==0 &&  FS_untilTP == 0
            intensity_sim_FLAG (i,tp) = intensityVector_0 (i,tp);
        end
        % % % sim_HA_NonFS        % INTENSITY FOR HA in 0F
        if FS ==0 && FS_untilTP == 0
            intensity_sim_HA_NonFS (i,tp) = intensityVector_HA (i,tp);
        end
    end
end

%% removing lines with empty intensities.
rm_zeros_intensity_sim_HA_FS = intensity_sim_HA_FS(any(intensity_sim_HA_FS,2),:);
rm_zeros_intensity_sim_1F_FS = intensity_sim_1F_FS(any(intensity_sim_1F_FS,2),:);
rm_zeros_intensity_sim_FLAG = intensity_sim_FLAG(any(intensity_sim_FLAG,2),:);
rm_zeros_intensity_sim_HA_NonFS = intensity_sim_HA_NonFS(any(intensity_sim_HA_NonFS,2),:);

%% Calculating total mean intensity

% INTENSITY FOR HA in FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_sim_HA_FS)==0
    if size(rm_zeros_intensity_sim_HA_FS,1)>1
        rm_zeros_intensity_sim_HA_FS = rm_zeros_intensity_sim_HA_FS./max(max(rm_zeros_intensity_sim_HA_FS)); % First normalization with respect to the maximum value in all repetitions
        sim_HA_FS = mean (rm_zeros_intensity_sim_HA_FS ); % Calculating mean intensity
        sim_HA_FS = sim_HA_FS./mean(sim_HA_FS(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_HA_FS = sim_HA_FS + BG_rawData_HA_inFS; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_HA_FS = sim_HA_FS./mean(sim_HA_FS(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_HA_FS   = std(rm_zeros_intensity_sim_HA_FS(:,:))./sqrt(rawData_HA_inFS(:,1)');
    else
        sim_HA_FS = zeros(1, numberOfTimePoints);
        sim_Err_HA_FS   = zeros(1, numberOfTimePoints);
    end
else
    sim_HA_FS = zeros(1, numberOfTimePoints);
    sim_Err_HA_FS   = zeros(1, numberOfTimePoints);
end

% INTENSITY FOR HA in NON-FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_sim_HA_NonFS)==0
    if size(rm_zeros_intensity_sim_HA_NonFS,1)>1
        rm_zeros_intensity_sim_HA_NonFS = rm_zeros_intensity_sim_HA_NonFS./max(max(rm_zeros_intensity_sim_HA_NonFS)); % First normalization with respect to the maximum value in all repetitions
        sim_HA_NonFS = mean (rm_zeros_intensity_sim_HA_NonFS ); % Calculating mean intensity
        sim_HA_NonFS = sim_HA_NonFS./mean(sim_HA_NonFS(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_HA_NonFS = sim_HA_NonFS + BG_rawData_HA_inNonFS; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_HA_NonFS = sim_HA_NonFS./mean(sim_HA_NonFS(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_HA_NonFS = std(rm_zeros_intensity_sim_HA_NonFS(:,:))./sqrt(rawData_HA_inNonFS(:,1)');
    else
        sim_HA_NonFS = zeros(1, numberOfTimePoints);
        sim_Err_HA_NonFS = zeros(1, numberOfTimePoints);
    end
else
    sim_HA_NonFS = zeros(1, numberOfTimePoints);
    sim_Err_HA_NonFS = zeros(1, numberOfTimePoints);
end

% INTENSITY FOR -1 in FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_sim_1F_FS)==0
    if size(rm_zeros_intensity_sim_1F_FS,1)>1
        rm_zeros_intensity_sim_1F_FS = rm_zeros_intensity_sim_1F_FS./max(max(rm_zeros_intensity_sim_1F_FS)); % First normalization with respect to the maximum value in all repetitions
        sim_1F_FS = mean (rm_zeros_intensity_sim_1F_FS ); % Calculating mean intensity
        sim_1F_FS = sim_1F_FS./mean(sim_1F_FS(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_1F_FS = sim_1F_FS + BG_rawData_SunTag_inFS; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_1F_FS = sim_1F_FS./mean(sim_1F_FS(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_1F_FS   = std(rm_zeros_intensity_sim_1F_FS(:,:))./sqrt(rawData_SunTag_inFS(:,1)');
    else
        sim_1F_FS = zeros(1, numberOfTimePoints);
        sim_Err_1F_FS   = zeros(1, numberOfTimePoints);
    end
else
    sim_1F_FS = zeros(1, numberOfTimePoints);
    sim_Err_1F_FS   = zeros(1, numberOfTimePoints);
end

% INTENSITY FOR FLAG NON-FRAMESHIFTING SPOTS
if isempty(rm_zeros_intensity_sim_FLAG)==0
    if size(rm_zeros_intensity_sim_FLAG,1)>1
        rm_zeros_intensity_sim_FLAG = rm_zeros_intensity_sim_FLAG./max(max(rm_zeros_intensity_sim_FLAG)); % First normalization with respect to the maximum value in all repetitions
        sim_FLAG = mean (rm_zeros_intensity_sim_FLAG ); % Calculating mean intensity
        sim_FLAG = sim_FLAG./mean(sim_FLAG(1:time_Normalization)); % First normalization with respect the the first experimental time points
        sim_FLAG = sim_FLAG + BG_rawData_FLAG_Only; % ADDING BACKGROUND INTENSITY FROM EXPERIMENTAL DATA.
        sim_FLAG = sim_FLAG./mean(sim_FLAG(1:time_Normalization)); % Second normalization with respect the the first experimental time points
        sim_Err_Flag = std(rm_zeros_intensity_sim_FLAG(:,:))./sqrt(rawData_FLAG_Only(:,1)');
    else
        sim_FLAG = zeros(1, numberOfTimePoints);
        sim_Err_Flag =  zeros(1, numberOfTimePoints);
    end
else
    sim_FLAG = zeros(1, numberOfTimePoints);
    sim_Err_Flag =  zeros(1, numberOfTimePoints);
end




end
