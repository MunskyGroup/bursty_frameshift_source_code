function [generated_geneSequence, typeOfTag,tagPositions, tagPositions_HA] = sequenceAnalyzer(str_geneSequence)

%% Defining the gene sequences
whole_gene = str_geneSequence;
typeOfTag = [];
tagPositions = [];
geneSequence = [];
tagPositions_HA = [1,1];

%% Defining the distinct tag sequences
T_SunTag = 'EELLSKNYHLENEVARLKK';
T_Flag = 'DYKDDDDK';
T_Hemagglutinin = 'YPYDVPDYA';

try
%% Detecting all the possible ORFs in the gene
orfPositions = seqshoworfs(whole_gene.Sequence,'MinimumLength',100, 'nodisplay', 'true','Frames',1);

%% Read each possible ORF and trasnlate those sequences from DNA to protein
counter = 1;
for i =1: length (orfPositions)
for j =1: length (orfPositions(i).Start)
try
temporal_ORF{counter} = whole_gene.Sequence (orfPositions(i).Start(j):orfPositions(i).Stop(j)-1);
temporal_Proteins {counter} =   nt2aa(temporal_ORF{counter});
counter = counter+1;
catch
end
end
end

%% Find if one of those proteins contains the tags.

position_SunTag = strfind(temporal_Proteins,T_SunTag);
position_Flag = strfind(temporal_Proteins,T_Flag);
position_Hemagglutinin = strfind(temporal_Proteins,T_Hemagglutinin);

%# find empty cells in the cell array
emptyCells_SunTag = cellfun(@isempty,position_SunTag);
emptyCells_Flag = cellfun(@isempty,position_Flag);
emptyCells_Hemagglutinin = cellfun(@isempty,position_Hemagglutinin);

% selecting the Number of Protein that contains the Tag region
selected_counter_SunTag = find(~emptyCells_SunTag,1);
selected_counter_Flag = find(~emptyCells_Flag,1);
selected_counter_Hemagglutinin = find(~emptyCells_Hemagglutinin,1);

%# remove empty cells
position_SunTag(emptyCells_SunTag) = [];
position_Flag(emptyCells_Flag) = [];
position_Hemagglutinin(emptyCells_Hemagglutinin) = [];


%% Finding the type of Tags
if isempty(position_SunTag) ==0
typeOfTag = 'SunTag';
tagPositions = position_SunTag{1};
geneSequence = upper(temporal_ORF{selected_counter_SunTag});
end

if isempty(position_Flag) == 0
typeOfTag = 'Flag';
tagPositions = position_Flag{1};
geneSequence = upper(temporal_ORF{selected_counter_Flag});
end

if isempty(position_Hemagglutinin) ==0
typeOfTag = 'HA';
tagPositions = position_Hemagglutinin{1};
geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
end

if isempty(position_Hemagglutinin) ==0 && isempty(position_Flag) == 0
typeOfTag = 'HA_Flag';
tagPositions = position_Flag{1};
tagPositions_HA = [position_Hemagglutinin{1}];
geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
end

if isempty(position_Hemagglutinin) ==0 && isempty(position_SunTag) ==0
typeOfTag = 'HA_SunTag';
tagPositions = [position_SunTag{1}];
tagPositions_HA = [position_Hemagglutinin{1}];
geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
end

generated_geneSequence.Sequence= geneSequence;
catch
end


end



