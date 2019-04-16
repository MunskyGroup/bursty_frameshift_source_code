function [ribosomesPerRNA_BeforeFSS, IntensityVectors_0F,IntensityVectors_1F,IntensityVectors_HA,time, geneLength_0F,geneLength_1F,pre_ribosomeDensity,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = masterFunction(nonConsiderTime,geneFileName_0F, geneFileName_1F, k_elongationMean,position_FS, k_initiation,nRepetitions,totalSimulationTime, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, UsingBurstingModel,k_on, k_off, k_fss, k_s_fss, itms, def_times,elongationFast)
% zero frame
temp_geneSequence_Fasta_0F = fastaread(geneFileName_0F);
str_geneSequence_0F = temp_geneSequence_Fasta_0F;
[geneSequence_0F, ~,probePosition_0F,probePosition_HA ] = sequenceAnalyzer(str_geneSequence_0F);
% minus one frame
temp_geneSequence_Fasta_1F = fastaread(geneFileName_1F);
str_geneSequence_1F = temp_geneSequence_Fasta_1F;
[geneSequence_1F, ~,probePosition_1F,~] = sequenceAnalyzer(str_geneSequence_1F);

%%  Separating sequence in codons for 0 frame.
codons = geneSequence_0F.Sequence;
geneLength_0F = length(codons)/3;
counter = 1;
for i =1 : geneLength_0F
    separated_codons(i,:) = codons(counter : counter + 2);
    counter = counter + 3;
end
%% Elongation constant.
k_elongation=  zeros (1,geneLength_0F );
genomicCopyNumber;
tRNA_copyNumber = zeros (1,geneLength_0F );
for i = 1 : geneLength_0F
    field = separated_codons(i,:);
    tRNA_copyNumber (i) = strGenCopy.(field) ;
end
% calculating the mean values in the structure
mean_tRNACopyNumber = mean (mean(reshape(struct2array(strGenCopy),numel(fieldnames(strGenCopy)),[]),2));
for i = 1 : geneLength_0F
    k_elongation(i)= (tRNA_copyNumber (i)/ mean_tRNACopyNumber)* k_elongationMean;
    %    k_elongation(i)=  k_elongationMean;
end
k_elongation_0F = [k_initiation, k_elongation, 10];
%%  Separating sequence in codons for -1 frame.
codons = geneSequence_1F.Sequence;
geneLength_1F = length(codons)/3;
counter = 1;
for i =1 : geneLength_1F
    separated_codons(i,:) = codons(counter : counter + 2);
    counter = counter + 3;
end
%% Elongation constant.
k_elongation_1=  zeros (1,geneLength_1F );
genomicCopyNumber;
tRNA_copyNumber_1 = zeros (1,geneLength_1F );
for i = 1 : geneLength_1F
    field = separated_codons(i,:);
    tRNA_copyNumber_1 (i) = strGenCopy.(field) ;
end
% calculating the mean values in the structure
mean_tRNACopyNumber = mean (mean(reshape(struct2array(strGenCopy),numel(fieldnames(strGenCopy)),[]),2));
for i = 1 : geneLength_1F
    k_elongation_1(i)= (tRNA_copyNumber_1 (i)/ mean_tRNACopyNumber)* k_elongationMean;
end
k_elongation_1F = [k_initiation,k_elongation_1,10];
% time = linspace(0,totalSimulationTime,totalSimulationTime);
gene_total= geneLength_0F + (geneLength_1F - position_FS);
% Running the methods
[time,ribosomesPerRNA_BeforeFSS, pre_ribosomeDensity,IntensityVectors_0F,IntensityVectors_1F,IntensityVectors_HA,pre_numOfRibosomes_0,pre_numOfRibosomes_1] = solverSSA(nRepetitions, UsingBurstingModel, nonConsiderTime+totalSimulationTime, k_elongation_0F, k_elongation_1F,gene_total, position_FS, k_on, k_off, k_fss, k_s_fss, timePerturbationApplication, evaluatingFRAP, evaluatingInhibitor, nonConsiderTime , probePosition_0F,probePosition_1F,probePosition_HA,itms,def_times,elongationFast);
end
