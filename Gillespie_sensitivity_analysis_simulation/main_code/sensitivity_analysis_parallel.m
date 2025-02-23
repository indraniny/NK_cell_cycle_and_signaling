function sensitivity_analysis_parallel()
%This code is for sensitivity analysis, Dec 22, 6:20 PM , running for 40h,
%nrep=20 
all clear;
% Author: Darren Wethington, Indrani Nayak
% Developer's note: Check for %%% ___ %%% for spots where edits can be made
% for individual changes

% This is designed to run a gillespie simulation for a very long time in
% order to get an "initial condition".  This is "before" stimulation occurs
% which would perturb the cells from a steady state condition.

% Add number of runs for the sensitivity analysis
num_of_runs=200
output_file = 'sensitivity_simulation_results.txt';

% Number of paramter sets want to run for the sensitivity analysis
for idx = 1:num_of_runs
    
    
    n_reps = 20; % number of initial cells, or "repetitions".  A true steady
    % state with no noise should only need 1, but this can help
  
    %generate uniform random number for sensitivity analysis  
    r = 20.0
    delta = rand*(0.9-0.03)+0.03;
        
    kcat=10000.0
    rd=36.0
    
    kon=rand*(100.0-5.0)+5.0
    koff=3000.0

    
    



    cycle=[.07 .1 .2 .7]; % 14 h, 10 h, 5 h, 1.42 h
    disp('Cell cycle transition rates in hours (G1->S, S->G2, G2->M, M-> G1) :')
    disp(1./cycle)

        %%% Stoichiometry and Reactant Matrices %%%
        
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

%         rxn_rates = ones(size(stoich,1),4) ;  %signaling reactions
% 
%         disp('1. binding vav+E = C, 2. Unbinding C=vav+E, 3. Phosphorylation C=pVav+E 4. Dephosphorylation pvav=vav, 5-6. Synthesis , 7-8. Degradation E, vav, 9. pVav degradation') 
%         disp('Rates in minutes:')

 
    
    % With:
    % Delete any existing parallel pool
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    % Create new parallel pool
    parpool(20);
    
   
    
    parfor z = 1:20
        %%% initialize some values %%%
        check_t =150; % some time before that to make sure the steady state has been achieved
        end_t =170; % how long the simulation will run for
        
        %times = [check_t end_t];
        times = [end_t];

       %%This is non-linear model paper plot for non-monotonic trend by changing phosphorylation rate in synthsis and mitosis stage+ with volume effect
        %%Monotonic increase
        % Reaction rates
            rxn_rates = ones(size(stoich, 1), 4); % signaling reactions
            % Define reaction rates (as per your system)
            

           
            
            rxn_rates(1,:) = kon; % Complex production (kon) [5-4500]
            rxn_rates(2,:) = koff; % Complex dissociation (koff) [450-7200]
            
            rxn_rates(3,:) = kcat; % Phosphorylation of complex C=pVAV+E (kcat) [36-7200000]
            %rxn_rates(3,2) = kcat/120;
            %rxn_rates(3,3:4) = kcat/120;

            rxn_rates(4,:) = rd; % De-phosphorylation of pVAV=VAV (rd)  [36-7200000]
             
            rxn_rates(5,:) = r; % Protein production (Enzyme) (r) [18-5350]
            rxn_rates(6,:) = r; % Protein production (S or Vav1) (r) [18-5350]
            rxn_rates(7:8,:) = delta; % Protein degradation (delta) [0.03-0.9]

            rxn_rates(9,:) = delta; % Complex degradation (delta) [0.03-0.9]
            rxn_rates(10,:) =delta; % Degradation of pVAV (delta) [0.03-0.9]
            
            % Account for volume correction
            V = [1 1.3 1.6 2.0];
            rxn_rates(1,:) = rxn_rates(1,:) ./ V;



      
       
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
    
    % After simulation completes, save results

    save_simulation_results(output_file,r,delta,kon,koff,kcat,rd);
    
    % Clean up parallel pool before next iteration
    delete(gcp('nocreate'));
end


end

% Helper function to save results (needs to be in separate file)
function parsave(savestring, cells_at_t, stoich, rxn_reactants, times, n_reps, cycle, rxn_rates, init, V)
    save(savestring, 'cells_at_t', 'stoich', 'rxn_reactants', 'times', 'n_reps', 'cycle', 'rxn_rates', 'init', 'V');
end

function save_simulation_results(filename,r,delta,kon,koff,kcat,rd)
    % Create format string for all 10 values (6delta + 11 data points)
    format_str = '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4e\t%.4e\t%.4e\n';
    
    % Prepare data row
    [data] = process_protein_averages_for_sensitivity();
    data2 = [r,delta,kon,koff,kcat,rd,data];
    
    % Open file in append mode
    fileID = fopen(filename, 'a');
    
    % Write header if file is empty
    if fileID == -1
        fileID = fopen(filename, 'w');
        fprintf(fileID, 'r\tDelta\tkon\tkoff\tkcat\trd\tG1\tS\tG2\tM\tG1_sem\tS_sem\tG2_sem\tM_sem\t p1\tp2\tp3\n');
    end
    
    % Write all data in one go using the format string
    fprintf(fileID, format_str, data2);
    
    % Close file
    fclose(fileID);
end

