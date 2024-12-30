function NL_models_paper_plot_gpu()

    all clear;
    % Author: Darren Wethington, Indrani Nayak
    % Developer's note: Check for %%% ___ %%% for spots where edits can be made
    % for individual changes

    % This is designed to run a gillespie simulation for a very long time in
    % order to get an "initial condition".  This is "before" stimulation occurs
    % which would perturb the cells from a steady state condition.

    %%% initialize some values %%%
    check_t = 150; % some time before that to make sure the steady state has been achieved
    end_t = 170; % how long the simulation will run for

    times = [check_t end_t];
    n_reps = 2; % number of initial cells, or "repetitions".  A true steady
    % state with no noise should only need 1, but this can help

    cycle = [.07 .1 .2 .7]; % 14 h, 10 h, 5 h, 1.42 h
    %disp('Cell cycle transition rates in hours (G1->S, S->G2, G2->M, M-> G1) :')
    disp(1./cycle)

    %%% Stoichiometry and Reactant Matrices %%%
    stoich = [-1 -1 1 0; % rate k1
        1 1 -1 0;
        1 0 -1 1;
        0 1 0 -1;
        1 0 0 0;
        0 1 0 0;
        -1 0 0 0;
        0 -1 0 0;
        0 0 -1 0;
        0 0 0 -1];

    rxn_reactants = [1 1 0 0; % propensity = k1*pi(x_i^r_i)
        0 0 1 0;
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0;
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];

    % Reaction rates and parameters
    rxn_rates = ones(size(stoich, 1), 4); % signaling reactions
    r = 20.0; %(r) [18-5350]
    delta = 0.5; %(r) [0.03-0.9]
    kcat = 10000.0; %(kcat) [36-7200000]
    rd = 36.0; % De-phosphorylation of pVAV=VAV (rd)  [36-7200000]
    kon = 10.0; % Complex production (kon) [5-4500]
    koff = 3000.0; % Complex dissociation (koff) [450-7200]

    rxn_rates(1,:) = kon; % Complex production (kon) [5-4500]
    rxn_rates(2,:) = koff; % Complex dissociation (koff) [450-7200]
    rxn_rates(3,:) = kcat; % Phosphorylation of complex C=pVAV+E (kcat) [36-7200000]
    rxn_rates(3,2) = kcat/100;
    rxn_rates(4,:) = rd; % De-phosphorylation of pVAV=VAV (rd) [36-7200000]
    rxn_rates(5,:) = r; % Protein production (Enzyme) (r) [18-5350]
    rxn_rates(6,:) = r; % Protein production (S or Vav1) (r) [18-5350]
    rxn_rates(7:8,:) = delta; % Protein degradation (delta) [0.03-0.9]
    rxn_rates(9,:) = delta; % Complex degradation (delta) [0.03-0.9]
    rxn_rates(10,:) = delta; % Degradation of pVAV (delta) [0.03-0.9]

    % Account for volume correction
    V = [1 1.3 1.6 2.0];
    rxn_rates(1,:) = rxn_rates(1,:) ./ V;

    disp('alpha ')
    alpha = kcat/(delta + rd);

    disp('E0 and S0')
    E0 = rxn_rates(5, 1) / rxn_rates(7, 1);
    S0 = rxn_rates(6, 1) / rxn_rates(8, 1);
    bE0_plus_S0 = ((1 + alpha) * E0) + S0;
    E0_times_S0 = E0 * S0;

    disp('KDprime')
    kdp = (koff + delta + kcat) / kon;

    % Delete any existing parallel pool        
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    parpool(20);

    % Move stoichiometry, reaction rates, and other variables to GPU
    stoich_gpu = gpuArray(stoich);
    rxn_reactants_gpu = gpuArray(rxn_reactants);
    rxn_rates_gpu = gpuArray(rxn_rates);
    V_gpu = gpuArray(V);

    % Running z number of independent parallel runs
    parfor z = 1:20
        % initial protein and cell distribution counts
        init = gpuArray([80 80 0 0 1 0 0 0 0]);  % Move initial condition to GPU

        local_cells_at_t = cell(1, length(times));
        local_cell_mean = zeros(n_reps, size(init, 2), 'gpuArray');  % Move mean matrix to GPU

        for t = 1:length(times)
            local_cells_at_t{t} = [];
            disp(['Worker ' num2str(z) ' processing time ' num2str(t)]);
            tic;
            for k = 1:n_reps
                % Use local variables instead of shared ones
                local_cells_at_t{t} = [local_cells_at_t{t}; KMC_for_ODEs(init, times(t),[], cycle, stoich_gpu, rxn_rates_gpu, rxn_reactants_gpu)];

                local_cell_mean(k, :) = mean(local_cells_at_t{t}, 1);

                if k > 3
                    if max(abs((local_cell_mean(k,:) - local_cell_mean(k-2,:)) ./ local_cell_mean(k,:))) < 0.01
                        break
                    end
                end
            end
            toc;
        end

        % Save results for this specific worker/iteration
        savestring = strcat('init_cond_all_rate_scan_', string(z));
        % Use local variables when saving
        parsave(savestring, local_cells_at_t, stoich_gpu, rxn_reactants_gpu, times, n_reps, cycle, rxn_rates_gpu, init, V_gpu);
    end
end

% Helper function to save results (needs to be in separate file)
function parsave(savestring, cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, init, V)
    save(savestring, 'cells_at_t', 'stoich', 'rxn_reactants', 'times', 'n_reps', 'cycle', 'rxn_rates', 'init', 'V');
end
