function NL_models_paper_plot()
%Neutral peak, Dec17, 11:45 PM, reducing rp in S by 100 times
all clear;
% Author: Darren Wethington, Indrani Nayak
% Developer's note: Check for %%% ___ %%% for spots where edits can be made
% for individual changes

% This is designed to run a gillespie simulation for a very long time in
% order to get an "initial condition".  This is "before" stimulation occurs
% which would perturb the cells from a steady state condition.

    %%% initialize some values %%%
    check_t =150; % some time before that to make sure the steady state has been achieved
    end_t =170; % how long the simulation will run for
    
    times = [check_t end_t];
    n_reps = 20; % number of initial cells, or "repetitions".  A true steady
    % state with no noise should only need 1, but this can help

    cycle=[.07 .1 .2 .7]; % 14 h, 10 h, 5 h, 1.42 h
    %disp('Cell cycle transition rates in hours (G1->S, S->G2, G2->M, M-> G1) :')
    disp(1./cycle)

    %%% Stoichiometry and Reactant Matrices %%%
    % Each row in this matrix corresponds to a reaction
    % Each column corresponds to a species in the system
    % The stoich matrix corresponds to the changes in each species for each
    % reaction
    % The rxn_reactants matrix corresponds to which proteins should be used to
    % calculate propensities for each reaction
    % They are flexibile; design your own system here!!
    
    % Here I've designed a system where enzymes bind/unbind to substrates to
    % form a complex C, and phorphorylation/unbinding happens to the substrate
    % [E S C P]
    
        %E+S=C %signaling reaction 1: binding
        %C=E+S %signaling reaction 2: unbinding
        %C=E+P % signaling reaction 3: protein phopshorylation
        %P=S %signaling reaction 4 : protein dephosphorylation
        
        %phi=E % Protein synthesis (5-6)
        %phi=S % S synthesis
        %E=phi %Protein degrdation (7-8)
        %S=phi %S degradation
        %C=phi % Complex degradation
        %P=phi %pVAV degradation
        

    %Here, E S indicates free Enzyme and Complex
    % Total Enzyme = Etot=E+C
    %Total substrate= Stot=S+C+P    
    % [E S C P]
    % [E vav E-vav pvav]
    
        
    stoich = [-1 -1 1 0; % rate k1
        1 1 -1 0;
        1 0 -1 1;
        0 1 0 -1;
        1 0 0 0;
        0 1 0 0
        -1 0 0 0
        0 -1 0 0
        0 0 -1 0
        0 0 0 -1];
    disp(stoich)
    rxn_reactants = [1 1 0 0; % propensity = k1*pi(x_i^r_i)
        0 0 1 0;
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1];

    
    
    % Reaction stoich for proteins that don't obey heritibility go here
    uninherit_rxn = [];
    % Example: uninherit_rxn = [0 1 0 -1]; %corresponds to phosphorylated
    % protein becoming substrate upon division
    
    %%% Assign reaction rates for each reaction %%%
    % reaction rates are a rxp matrix, where r is number of reactions and p
    % is number of cell cycle stages.  If you would like to alter a
    % reaction rate for a particular cell cycle stage, alter the
    % appropriate column.  For instance, to increase binding in the S
    % phase, use rxn_rates(1,2) = rxn_rates(1,2)*10; (or whatever works)
    
    rxn_rates = ones(size(stoich,1),4) ;  %signaling reactions

    %disp('1. binding vav+E = C, 2. Unbinding C=vav+E, 3. Phosphorylation C=pVav+E 4. Dephosphorylation pvav=vav, 5-6. Synthesis , 7-8. Degradation E, vav, 9. pVav degradation') 
    %disp('Rates in minutes:')



   %%This is non-linear model paper plot for non-monotonic trend by changing phosphorylation rate in synthsis and mitosis stage+ with volume effect
    %%Monotonic increase
    % Reaction rates
        rxn_rates = ones(size(stoich, 1), 4); % signaling reactions
        % Define reaction rates (as per your system)
        
        r= 20.0; %(r) [18-5350]
        delta=0.5; %(r) [0.03-0.9]
        
        kcat=10000.0; %(kcat) [36-7200000]
        rd= 36.0; % De-phosphorylation of pVAV=VAV (rd)  [36-7200000]
        
        kon=10.0; % Complex production (kon) [5-4500]
        koff = 3000.0; % Complex dissociation (koff) [450-7200]
       
        
        rxn_rates(1,:) = kon; % Complex production (kon) [5-4500]
        rxn_rates(2,:) = koff; % Complex dissociation (koff) [450-7200]
        
        rxn_rates(3,:) = kcat; % Phosphorylation of complex C=pVAV+E (kcat) [36-7200000]
        
        

        rxn_rates(4,:) = rd; % De-phosphorylation of pVAV=VAV (rd)  [36-7200000]
         
        rxn_rates(5,:) = r; % Protein production (Enzyme) (r) [18-5350]
        rxn_rates(6,:) = r; % Protein production (S or Vav1) (r) [18-5350]
        rxn_rates(7:8,:) = delta; % Protein degradation (delta) [0.03-0.9]

        rxn_rates(9,:) = delta; % Complex degradation (delta) [0.03-0.9]
        rxn_rates(10,:) =delta; % Degradation of pVAV (delta) [0.03-0.9]
        
        % Account for volume correction
        V = [1 1.3 1.6 2.0];
        rxn_rates(1,:) = rxn_rates(1,:) ./ V;

        disp('alpha ')
        alpha=kcat/(delta+rd)
        
        
        disp('E0 and S0')
        E0=rxn_rates(5,1)/rxn_rates(7,1)
        S0=rxn_rates(6,1)/rxn_rates(8,1)
        bE0_plus_S0=((1+alpha)*E0)+S0   
        E0_times_S0=E0*S0
   
        disp('KDprime')
        kdp=(koff+delta+kcat)/kon

% Delete any existing parallel pool        
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool(20);
        
        
%running z number of independent parallel runs
parfor z = 1:20
   
    % initial protein and cell distribution counts.  This shouldn't matter for
    % steady state calculations; but you may need to edit it if there are
    % different numbers of species or cell cycle stages. %NOTE: I'm seeing this
    % doesn't match the 'init' for the ODE solutions.  Might need to be edited
    %Initially there is 3 cells and protein per cell is 80, total E and S
    %=240 initially
    % init = [#E,#S,#C,#P,#G1,#S,#G2,#M,start_time]
    init = [80 80 0 0 1 0 0 0 0];
    
    % Create worker-specific variables
    local_cells_at_t = cell(1, length(times));
    local_cell_mean = zeros(n_reps, size(init, 2));
    
    for t=1:length(times)
        local_cells_at_t{t} = [];
        disp(['Worker ' num2str(z) ' processing time ' num2str(t)]);
        tic;
        for k=1:n_reps
            % Use local variables instead of shared ones
            local_cells_at_t{t} = [local_cells_at_t{t}; KMC_for_ODEs(init,times(t),uninherit_rxn,cycle,stoich,rxn_rates,rxn_reactants)];
            
            local_cell_mean(k,:) = mean(local_cells_at_t{t},1);
            
            if k>3
                if max(abs((local_cell_mean(k,:)-local_cell_mean(k-2,:))./local_cell_mean(k,:)))<0.01
                    break
                end
            end
        end
        toc;
    end
    
    % Save results for this specific worker/iteration
    savestring = strcat('init_cond_all_rate_scan_',string(z));
    % Use local variables when saving
    parsave(savestring, local_cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, init, V);
end


end

% Helper function to save results (needs to be in separate file)
function parsave(savestring, cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, init, V)
    save(savestring, 'cells_at_t', 'stoich', 'rxn_reactants', 'times', 'n_reps', 'cycle', 'rxn_rates', 'init', 'V');
end

