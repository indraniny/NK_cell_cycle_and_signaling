function NL_model_configuration_generator(t_end,initial_cell_no,number_of_workers,kcat, all_init)
%cells_at_t file [E, S, C, P, G1, S, G2, M, time]

all clear;
% Author: Darren Wethington, Indrani Nayak
% Developer's note: Check for %%% ___ %%% for spots where edits can be made
% for individual changes

% This is designed to run a gillespie simulation for a very long time in
% order to get an "initial condition".  This is "before" stimulation occurs
% which would perturb the cells from a steady state condition.


    times=[t_end];
    n_reps=initial_cell_no;
    no_of_worker=number_of_workers;


    % state with no noise should only need 1, but this can help

    cycle=[.07 .1 .2 .7]; % 14 h, 10 h, 5 h, 1.42 h
    %disp('Cell cycle transition rates in hours (G1->S, S->G2, G2->M, M-> G1) :')
    disp(1./cycle)

    %%% Stoichiometry and Reactant Matrices %%%

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
    % [E S C=ES P=S*]
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
        
        r=5350.0; %(r) [18-5350]
        delta=0.03; %(r) [0.03-0.9]
        
        %kcat=10.0; %(kcat) [36-7200000]
        rd= 2500.0; % De-phosphorylation of pVAV=VAV (rd)  [36-7200000]
        
        kon=4500.0; % Complex production (kon) [5-4500]
        koff = 450.0; % Complex dissociation (koff) [450-7200]
       
        
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


% Delete any existing parallel pool        
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool();
        
        
%running z number of independent parallel runs
parfor z = 1:no_of_worker
   
    % initial protein and cell distribution counts.  This shouldn't matter for
    % steady state calculations; but you may need to edit it if there are
    % different numbers of species or cell cycle stages. 
   
    % init = [#E,#S,#C,#P,#G1,#S,#G2,#M,start_time]
    %init = [80 80 0 0 1 0 0 0 0]; # it was an old configuration
    

    init=all_init{z};

    size(init);
    init;
    % Create worker-specific variables
    local_cells_at_t = cell(1, length(times));
    local_cell_mean = zeros(n_reps, size(init, 2));
    
    if size(init, 1) ~= n_reps || size(init, 2) ~= 9
        error('init must be a matrix of size [n_reps x 9]');
    end

    for t=1:length(times)
        local_cells_at_t{t} = [];
        disp(['Worker ' num2str(z) ' processing time ' num2str(t)]);
        tic;
        for k=1:n_reps
            
            local_init=init(k,:);

            % disp('---')
            % disp(k)
            % disp(local_init)
            % disp(rxn_rates)
            % disp(n_reps)
            % Use local variables instead of shared ones
            local_cells_at_t{t} = [local_cells_at_t{t}; KMC_for_ODEs(local_init,times(t),uninherit_rxn,cycle,stoich,rxn_rates,rxn_reactants)];
            
            local_cell_mean(k,:) = mean(local_cells_at_t{t},1);
            
            % if k>3
            %     if max(abs((local_cell_mean(k,:)-local_cell_mean(k-2,:))./local_cell_mean(k,:)))<0.01
            %         break
            %     end
            % end
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

