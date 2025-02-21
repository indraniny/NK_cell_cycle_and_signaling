%function get_init_condition_Indrani_ODE()
all clear;
% Author: Darren Wethington
% Developer's note: Check for %%% ___ %%% for spots where edits can be made
% for individual changes

% This is designed to run a gillespie simulation for a very long time in
% order to get an "initial condition".  This is "before" stimulation occurs
% which would perturb the cells from a steady state condition.

% If a parameter scan is desired, you can use this counter to loop.  Make
% sure parameters are randomly sampled somewhere for this to work.
for z=1
    
    %%% initialize some values %%%
    end_t = 100; % how long the simulation will run for
    %check_t =150; % some time before that to make sure the steady state has been achieved
    %times = [check_t end_t];
    times=0.05:0.3:50;
    %times = [10 20 30 40 50]
    n_reps = 400; % number of initial cells, or "repetitions".  A true steady
    % state with no noise should only need 1, but this can help
    % reduce noise.
    
    %%
    %For Eukaryotic cells
    %G1 = 6-12 hours
    %S = 6-8 hours
    %G2 = 3-4 hours
    % M = 1 hour
    %%
    

    %cycle=[.1 .1 .1 .1]; % 14 h, 10 h, 5 h, 1.42 h
    %cycle=[.15 .1 .45 .5]; % 14 h, 10 h, 5 h, 1.42 h
    cycle=[0.07 0.1 0.2 0.7]; 
    disp('Cell cycle transition rates in hours (G1->S, S->G2, G2->M, M-> G1) :')
    disp(1./cycle)
 
    
    %%% Stoichiometry and Reactant Matrices %%%
    % Each row in this matrix corresponds to a reaction
    % Each column corresponds to a species in the system
    % The stoich matrix corresponds to the changes in each species for each
    % reaction
    % The rxn_reactants matrix corresponds to which proteins should be used to
    % calculate propensities for each reaction
    % They are flexibile; design your own system here!!
    
 
     %[A A*]
     %reaction 1: phosphorylation kp
     %reaction 2: dephosphorylation kd;
     %reaction 3: Synthesis of A;
     %reaction 4-5: degradation of A, A*
     
     
    %% Linear protein model where A->A*,A*->A,0->A,A->0,A*->0 the stoich and rxn_reactants will be different
    stoich=[-1 1;
            1 -1;
            1  0;
            -1 0;
            0 -1];

    disp(stoich)
    
    rxn_reactants=[1 0;% propensity = k1*pi(x_i^r_i)
                    0 1;
                    0 0;
                    1 0;
                    0 1];


 
    %%
    % Reaction stoich for proteins that don't obey heritibility go here,
    % Present code only assumed heritability
    uninherit_rxn = [];
    % Example: uninherit_rxn = [+1 -1]; %corresponds to phosphorylated
    % protein becoming substrate upon division, 
    
    %%% Assign reaction rates for each reaction %%%
    % reaction rates are a rxp matrix, where r is number of reactions and p
    % is number of cell cycle stages.  If you would like to alter a
    % reaction rate for a particular cell cycle stage, alter the
    % appropriate column.  For instance, to increase binding in the S
    % phase, use rxn_rates(1,2) = rxn_rates(1,2)*10; (or whatever works)
    
    
    %%Linear model
     rxn_rates = ones(size(stoich,1),4) ;  %signaling reactions
     rxn_rates=rxn_rates*0.5
     
     
     %rxn_rates(1,2)=5.0 %phosphorylation rates
     
     rxn_rates(1,:)=1.0 %phosphorylation rates
     rxn_rates(2,:)=1.0 %dephosph rate
     rxn_rates(3,:)=1.0 %synthesis rate
     rxn_rates(4:5,:)=0.01 %degradation rate
     
     %rxn_rates(1,2)=10.0 %phosphorylation rates
     


    % This code is for matching ODE simulations to gillespie.  It is commented
    % out because this code involves volume corrections, which lead to a
    % mismatch between the stochastic solutions and deterministic.  Feel free
    % to uncomment for any exploration, especially if there is no volume
    % correction.
    
    init = zeros((length(stoich(:,1))+1)*length(cycle),1)';
    init(1:2) = [300,0];
    init(end-length(cycle)+1) = 1
    
    stoich
    rxn_reactants
    rxn_rates
    cycle
    uninherit_rxn
    end_t

    sol = convert_gillespie_setup_to_ode(stoich,rxn_reactants,[300 0 1 0 0 0 0],rxn_rates,cycle,uninherit_rxn,end_t);
    odesol = deval(sol,times)
    
    %print(odesol)
    save('odesol.mat', 'odesol');
    
    cell_values_ode = odesol(9:12,:);
    tot_A_levels_ode = odesol(1:2:8,:);
    tot_A_phosph_levels_ode = odesol(2:2:8,:);
    
    
    A_levels_ode = tot_A_levels_ode./cell_values_ode;
    A_phosph_levels_ode = tot_A_phosph_levels_ode./cell_values_ode;
    
    save('A_levels_ode.mat', 'A_levels_ode');
    save('A_phosph_levels_ode.mat', 'A_phosph_levels_ode');
    save('times.mat', 'times');

%     %plotting protein A with time
    figure()
    plot(times, A_levels_ode(1,:), 'DisplayName', 'Protein G1')
    hold on
    plot(times, A_levels_ode(2,:), 'DisplayName', 'Protein S')
    plot(times, A_levels_ode(3,:), 'DisplayName', 'Protein G2')
    plot(times, A_levels_ode(4,:), 'DisplayName', 'Protein M')

    % Add a title
    title('ODE: a Across Cell Cycle Stages')
    % Add a legend
    legend('show')
    % Optionally, label the axes
    xlabel('Time')
    ylabel('Avg. Protein Level')

    %plotting protein A* with time
    figure()
    plot(times, A_phosph_levels_ode(1,:), 'DisplayName', 'Protein G1')
    hold on
    plot(times, A_phosph_levels_ode(2,:), 'DisplayName', 'Protein S')
    plot(times, A_phosph_levels_ode(3,:), 'DisplayName', 'Protein G2')
    plot(times, A_phosph_levels_ode(4,:), 'DisplayName', 'Protein M')

    % Add a title
    title('ODE: a* Across Cell Cycle Stages')
    % Add a legend
    legend('show')
    % Optionally, label the axes
    xlabel('Time')
    ylabel('Avg. Protein Level')
% %     
%     % initial protein and cell distribution counts.  This shouldn't matter for
    % steady state calculations; but you may need to edit it if there are
    % different numbers of species or cell cycle stages. %NOTE: I'm seeing this
    % doesn't match the 'init' for the ODE solutions.  Might need to be edited
    %Initially there is 3 cells and protein per cell is 80, total E and S
    %=240 initially
    % init = [#A,#A*,#G1,#S,#G2,#M,start_time]
    

    
    init = [300 0 1 0 0 0 0];
    %init = [100 0 1 0];
    
    
    for t=1:length(times) % loop for t trajectory
        cells_at_t{t} = []; %empty array initialization
        disp(t) %
        tic;
        for k=1:n_reps % loop for repetitions at t
            %appends the resulting output to the cell array cells_at_t{t}
            cells_at_t{t} = [cells_at_t{t}; KMC_for_ODEs(init,times(t),uninherit_rxn,cycle,stoich,rxn_rates,rxn_reactants)]; % evaluate KMC
            % I want to check whether steady state is reached and I can
            % terminate the loop.  This is mostly when the number of
            % reps to reduce noise is large. 
            cell_mean(k,:) = mean(cells_at_t{t},1);
            %This assigns the row vector of mean values to the k-th row of the matrix cell_mean. Each row of cell_mean corresponds to a repetition (k) within the current time step (t).
            % if more reps, check whether steady state is already reached
            % and end the loop if so. (really ugly convergence check)
            if k>3
                if max(abs((cell_mean(k,:)-cell_mean(k-2,:))./cell_mean(k,:)))<0.01
                    break
                end
            end
        end
        
        toc;
    end
    
    % save the output according to the outer param scan loop. This line constructs a filename string by concatenating the prefix 'init_cond_all_rate_scan_' with the string representation of the variable z.
    savestring = strcat('init_cond_all_rate_scan_',string(z));
    %This line saves the variables cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, and init to a MAT-file with the filename specified by 
    
    %This line saves the variables cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, and init to a MAT-file with the filename specified by 
    save(savestring,'cells_at_t','stoich','rxn_reactants','times','n_reps','cycle','rxn_rates','init')
    
end

%end