function [ribosomesPerRNA_BeforeFSS,pre_ribosomeDensity, pre_numOfRibosomes_0,pre_numOfRibosomes_1,IntensityVectors_0F,IntensityVectors_1F,IntensityVectors_HA] = SSAtoIntensity(nonConsiderTime,t_array,gene_total,RibsomePositions,probePosition_0F,probePosition_1F,probePosition_HA,k_elongation_0F,position_FS )
gene_Length_0F = length(k_elongation_0F(2:end-1));
numberOftimePoints_withoutBurningTime = length(t_array)- (nonConsiderTime);
edges = linspace(1,gene_total,gene_total);
%% Number of ribosomes per mRNA
pre_ribosomesPerRNA_BeforeFSS =  RibsomePositions(:,end);
pre_ribosomesPerRNA_BeforeFSS (pre_ribosomesPerRNA_BeforeFSS>position_FS) = 0;
ribosomesPerRNA_BeforeFSS = sum( pre_ribosomesPerRNA_BeforeFSS(:,:)>1);
% matrix
pre_ribosomeDensity = histcounts(RibsomePositions(:,:), edges);
X_output_0F = RibsomePositions;
X_output_0F (X_output_0F>gene_Length_0F) = 0;
X_output_1F =RibsomePositions;
X_output_1F (X_output_1F<gene_Length_0F) =0;
X_output_1F = (X_output_1F+ position_FS)- gene_Length_0F;

%% Distribution of ribosomes in the mRNA
pre_numOfRibosomes_0 = sum( X_output_0F(:,:)>1);
pre_numOfRibosomes_1 = sum( X_output_1F(:,:)>1);

%% Saving the intensity vectors in a cell array.
IntensityVectors_0F = 0;
IntensityVectors_1F = 0;
IntensityVectors_HA = 0;
for ip = 1:length(probePosition_0F)
    IntensityVectors_0F = IntensityVectors_0F+sum(X_output_0F(:,1:numberOftimePoints_withoutBurningTime)>=probePosition_0F(ip));
end
for ip = 1:length(probePosition_1F)
    IntensityVectors_1F = IntensityVectors_1F+sum(X_output_1F(:,1:numberOftimePoints_withoutBurningTime)>=probePosition_1F(ip));
end
for ip = 1:length(probePosition_HA)
    IntensityVectors_HA = IntensityVectors_HA+sum(X_output_0F(:,1:numberOftimePoints_withoutBurningTime)>=probePosition_HA(ip));
end

%% Maxing Intensity vectors a row vector of zeros in case the vector is zero
if numberOftimePoints_withoutBurningTime>1
    if length(IntensityVectors_0F)==1
        IntensityVectors_0F = zeros(1,numberOftimePoints_withoutBurningTime);
    end
    if  length(IntensityVectors_1F)==1
        IntensityVectors_1F = zeros(1,numberOftimePoints_withoutBurningTime);
    end
    if  length(IntensityVectors_HA)==1
        IntensityVectors_HA = zeros(1,numberOftimePoints_withoutBurningTime);
    end
    
    if length(pre_ribosomeDensity) ==1
        pre_ribosomeDensity = zeros(1,edges);
    end
    
    if length(ribosomesPerRNA_BeforeFSS) ==1
        ribosomesPerRNA_BeforeFSS = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
    
    if length(pre_numOfRibosomes_0) ==1
        pre_numOfRibosomes_0 = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
    
    if length(pre_numOfRibosomes_1) ==1
        pre_numOfRibosomes_1 = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
end


if numberOftimePoints_withoutBurningTime==1
    if isempty(pre_ribosomeDensity)==1
        pre_ribosomeDensity = zeros(1,edges);
    end
    
    if isempty(ribosomesPerRNA_BeforeFSS)==1
        ribosomesPerRNA_BeforeFSS = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
    
    if  isempty(pre_numOfRibosomes_0)==1
        pre_numOfRibosomes_0 = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
    
    if isempty(pre_numOfRibosomes_1)==1
        pre_numOfRibosomes_1 = zeros(1,numberOftimePoints_withoutBurningTime+1);
    end
end

end
