function [mean_ratio,sd_ratio]  = bootstrp_Intensities (intensityVector_0,intensityVector_1)
backgroundDetectionLimit =24;
numberOFSpots = size(intensityVector_0,1);
parfor j =1: 1000
sampIndex =randi(numberOFSpots,[1,numberOFSpots]); % Sampling vector.

intensityVector_0_sample = intensityVector_0(sampIndex,1); % selecting elements from the intensity vectors
intensityVector_1_sample = intensityVector_1(sampIndex,1); % selecting elements from the intensity vectors
%% Spot classification
I_y_0 = intensityVector_0_sample(intensityVector_0_sample >= backgroundDetectionLimit & intensityVector_1_sample <= backgroundDetectionLimit);
I_p_m1 = intensityVector_1_sample(intensityVector_1_sample >= backgroundDetectionLimit & intensityVector_0_sample <= backgroundDetectionLimit);
I_w_0 = intensityVector_0_sample(intensityVector_1_sample >= backgroundDetectionLimit & intensityVector_0_sample >= backgroundDetectionLimit);
I_w_m1 = intensityVector_1_sample(intensityVector_1_sample >= backgroundDetectionLimit & intensityVector_0_sample >= backgroundDetectionLimit);

%% Calculaitng the ratio 
b2g_ratio_mod(j) = (sum(I_w_m1)+sum(I_p_m1))/(sum(I_y_0)+sum(I_w_0));
end
%% Calculating the error
sd_ratio = std(b2g_ratio_mod);
mean_ratio = mean(b2g_ratio_mod);
end

