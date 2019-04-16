function [RibsomePositions] = SSA_BEM(UsingBurstingModel, position_FS, gene_total,maximum_Number_Ribosomes, X_State, k_elongation_0F, k_elongation_1F, t_array, k_on, k_off, k_fss, k_s_fss, timePerturbationApplication, evaluatingInhibitor, nonConsiderTime,elongationFast)
% Function SSA_BEM is intended to stochastically simulate translation with frameshifting.
%% Deffining parameters
exclusion = 9;  % exclusion volume
k_elongation = k_elongation_0F(2:end-1); % elongation constants
k_bind = k_elongation_0F(1);  % Initiation rate
k_termination = k_elongation_0F(end);    % Termination rate.
gene_Length_0F = length(k_elongation);  % number of codons.
%% Steady state for the HA system
t = t_array(1); % time array
%% Parameters for FS
if UsingBurstingModel ==0
    z = 1;
    k_on = 0;
    k_off = 0;
else
    z=0;
end
%% temporal solution for the elongations in the 2nd frame.
k_elongation_2F = k_elongation_1F(position_FS:end);
k_elongation_total = [k_elongation, k_elongation_2F];

%% Prealocating memory
number_TimePoints = length(t_array);
t_final = t_array(number_TimePoints);
X_output = zeros(maximum_Number_Ribosomes,number_TimePoints);
x_terms0 = zeros(1,1e4);ixt0=0;
x_terms1 = zeros(1,1E4);ixt1=0;
fss_to_end = gene_Length_0F-position_FS;  % length of the gene after the frameshift site (portion to be copied for -1 frame)

%% Run SSA
iteration = 1;
while t < t_final
    if t >= timePerturbationApplication && evaluatingInhibitor ==1
        Inhibitor_Condition = 0; % Set to 0 during its application
    else
        Inhibitor_Condition = 1; % Set to 1 before its application
    end
    %% Compute propensity functions
    Nribosomes = sum(X_State > 0);
    X_State(Nribosomes + 1,1) = 0;
    temp_Stoichiometry = eye(Nribosomes + 1);
    temp_Propensities = zeros(Nribosomes + 4,1);% N elongations/termination, 1 initiation, 2-z shifts, one FS.
    temp_Propensities(X_State > 0) = k_elongation_total(X_State(X_State > 0));
    
    if elongationFast ==1
        Lcut =  max(1,fss_to_end-1);
        %% BRIAN -- Determn. Termination
        % Change to make the termination reaction the one that passes the FSS.
        if (X_State(1)==(gene_Length_0F-Lcut))||(X_State(1)==(gene_total-Lcut))
            temp_Stoichiometry(:,1) = [X_State(2:Nribosomes+1);0] - X_State(1:Nribosomes+1);
            temp_Propensities(1) = k_termination;
        end
    else
        Lcut =  max(1,fss_to_end-400);
        %  termination in the normal case
        if X_State(1) == gene_Length_0F || X_State(1)>= gene_total-1
            temp_Stoichiometry(:,1) = [X_State(2:Nribosomes+1);0] - X_State(1:Nribosomes+1);
            temp_Propensities(1) = k_termination;
        end
    end
    
    %%
    % initiation
    if Nribosomes == 0 || X_State(Nribosomes) > exclusion % >=
        temp_Propensities(Nribosomes+1) = k_bind * Inhibitor_Condition;
    end
    
    %   Elongation
    if elongationFast ==1     %         case 'fast'
        qw = X_State;
        qw(qw>gene_Length_0F)=qw(qw>gene_Length_0F)-fss_to_end;
        elongation_Condition = qw(1:Nribosomes-1)~=(qw(2:Nribosomes)+exclusion);
    else         %        case 'orig'
        elongation_Condition = ...
            (X_State(2:Nribosomes)<gene_Length_0F)&(X_State(1:Nribosomes-1)<gene_Length_0F)&...
            ((X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1))|...   % both in first frame
            (X_State(2:Nribosomes)>gene_Length_0F)&(X_State(1:Nribosomes-1)>gene_Length_0F)&...
            ((X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1))|...   % both in second frame
            (X_State(2:Nribosomes)<gene_Length_0F)&(X_State(1:Nribosomes-1)>gene_Length_0F)&...
            ((X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1)-fss_to_end)|...  % first in second frame, follower in first frame
            (X_State(2:Nribosomes)>gene_Length_0F)&(X_State(1:Nribosomes-1)<gene_Length_0F)&...
            ((X_State(2:Nribosomes) + exclusion - fss_to_end) < X_State(1:Nribosomes-1));  % first in first frame, follower in second frame
    end
    if length(elongation_Condition)>=1
        if ~min(elongation_Condition)
            temp_Propensities(2:Nribosomes) = temp_Propensities(2:Nribosomes) .* elongation_Condition;
        end
    end
    % FS reactions for the zero frame
    FS_Condition = (X_State(1:Nribosomes)) ==  position_FS;  % X(i)=1 if that ribosome is in the FSS.  Otherwise zero.
    % FS reactions for the zero frame using the constitutive or bursting model
    if UsingBurstingModel ==1
        temp_Propensities(1:Nribosomes) = temp_Propensities(1:Nribosomes).* ~FS_Condition + (FS_Condition .*( k_fss .* (1-z)));
    else
        temp_Propensities(1:Nribosomes) = temp_Propensities(1:Nribosomes).* ~FS_Condition + (FS_Condition .*( k_fss .* (z)));
    end
    %% FS propensities
    temp_Propensities(Nribosomes+2:Nribosomes+4) = ...
        [any(FS_Condition)* k_s_fss* z,... % propensity to change to -1 frame
        k_on * (1-z),... % propensity to make z= 1;
        k_off * z]; % propensity to make z =0;
    %% Update time
    sum_Propensities = sum(temp_Propensities);
    t = t - log(rand) / sum_Propensities;
    %% Generate output
    while iteration <= number_TimePoints && t > t_array(iteration)
        %BRIAN        X_output(1:size(X_State,1),iteration) = X_State;
        X_output(1:maximum_Number_Ribosomes,iteration) = X_State;
        iteration = iteration + 1;
    end
    %% Update state
    if t < t_final
        %% Select reaction
        r2 = sum_Propensities * rand;
        i = 1;
        tmp = temp_Propensities(i);
        while tmp < r2
            i = i + 1;
            tmp = tmp + temp_Propensities(i);
        end
        %% update state
        if (i==1)&&((X_State(1)==gene_Length_0F-Lcut)||(X_State(1)==(gene_total-Lcut)))
            % BRIAN -- Termination Event.
            if X_State(1)<gene_Length_0F
                ixt0=ixt0+1;
                x_terms0(ixt0)=t;  % Record the time that the Ribosome passes the FSS region
            else
                ixt1=ixt1+1;
                x_terms1(ixt1)=t;  % Record the time that the Ribosome passes the FSS region
            end
            X_State(1:Nribosomes+1) = X_State(1:Nribosomes+1) + temp_Stoichiometry(:,i);
        elseif i<= length(temp_Propensities)- 3
            X_State(1:Nribosomes+1) = X_State(1:Nribosomes+1) + temp_Stoichiometry(:,i);
        elseif i == length(temp_Propensities)-2
            X_State (X_State==position_FS) = gene_Length_0F+1;
        elseif i == length(temp_Propensities)-1
            z = z+1;
        elseif i == length(temp_Propensities)
            z = z-1;
        end
        if UsingBurstingModel ==0
            z =1;
        end
    end
end

%% Saving ribosome position in time
RibsomePositions = X_output(:,nonConsiderTime:end);
RibsomePositions = RibsomePositions(sum(X_output,2)~=0,:);  % Truncate to remove zeros.

%% Reconstructing states after the FSS in the Fast elongation model
if elongationFast ==1
    Lcut =  max(1,fss_to_end-1);
    %%
    dt0 = cumsum(1./k_elongation(end-Lcut:end)); % Times of motion from codon to codon.
    ps = [gene_Length_0F-Lcut:gene_Length_0F];   % Positions of the codons from FSS to end.
    x_terms0=x_terms0(x_terms0~=0);    % Remove the zeros form the termination times.
    nte = length(x_terms0); % Number of ribosomes that terminated.
    RibPosn = zeros(nte,number_TimePoints);  % Allocate Rib position matrix to handle terminated Ribs.
    for i=1:length(x_terms0)  % Step through the terminated Ribs.
        ts = x_terms0(i)+dt0;  % Find times of codon transitions
        if ts(end)>=t_array(nonConsiderTime)
            v = (interp1(ts,ps,t_array,'previous')); % Find the codon posn at the printing time.
            v(isnan(v))=0;  % remove non-sensical values.
            RibPosn(i,:) = v;  % record positions.
        end
    end
    RibPosn = RibPosn(:,nonConsiderTime:end);   % remove the non-sense times.
    RibsomePositions = [RibsomePositions;RibPosn];% append the terminated Ribs.
    
    %% do same thing for the -1 Ribs.
    dt1 = cumsum(1./k_elongation_1F(end-Lcut:end));
    ps = [gene_total-Lcut:gene_total];
    x_terms1=x_terms1(x_terms1~=0);
    nte = length(x_terms1);
    RibPosn = zeros(nte,number_TimePoints);
    for i=1:length(x_terms1)
        ts = x_terms1(i)+dt1;
        if ts(end)<t_array(nonConsiderTime)
        else
            v = (interp1(ts,ps,t_array,'previous'));
            v(isnan(v))=0;
            RibPosn(i,:) = v;
        end
    end
    RibPosn = RibPosn(:,nonConsiderTime:end);
    RibsomePositions = [RibsomePositions;RibPosn];
end