function finished_cells = KMC_for_ODEs(op_m,end_t,uninherit_rxn,cycle_rates,rxn_stoich,rxn_rates,rxn_reactants)

% Author: Darren Wethington
% This code performs the Gillespie (or kinetic monte carlo KMC) simulations
% for the cell cycle signaling project.

% Inputs: op_m: the cells ("operating matrix")
%         end_t: the time to collect the output
%         uninherit_rxn: reaction stoichiometry for what happens to
%         proteins when cells divide
%         cycle_rates: cell cycle rates
%         rxn_stoich: stoichiometry
%         rxn_rates: reaction rates
%         rxn_reactants: reactants to be used in calculating propensities

% initialize some values
k = cycle_rates;
r = rxn_rates';
stoich = rxn_stoich;
% add cell cycle transitions to the possible reactions
stoich = [stoich zeros(length(stoich(:,1)),length(k))];
for i=2:length(k)
    stoich(end+1,length(rxn_stoich(1,:))+i) = 1;
    stoich(end,length(rxn_stoich(1,:))+i-1) = -1;
    rxn_reactants(end+1,end+1) = 1;
end
rxn_reactants(end+1,end+1) = 1;
stoich(end+1,end) = -1;
stoich(end,length(rxn_stoich(1,:))+1) = 1;
% because I account for time, I need to add an empty stoich column here
stoich(:,end+1) = zeros(length(stoich(:,1)),1);
k = repmat(k,length(k),1);

% The big loop. The logic is that each cell progresses through time
% independently. I calculate propensities and delta_t for each cell
% according to the reaction rules in the input. When a cell would enter my
% time of interest end_t, I remove it from the simulation and place it in
% finished_cells.

% initialize
finished_cells = [];
% while any cell still exists in the simulation
while any(op_m,'all')
    
    %determine each cell cycle stage
    state = zeros(length(op_m(:,1)),1);
    for i=1:length(state)
        state(i) = find(op_m(i,end-length(k):end-1));
    end
    %calculate all the possible reactions and their propensities
    options = propensities_from_reactants(rxn_reactants,op_m(:,1:end-1),[r k']',state);
    %cumulative sum of propensities so we can effectively choose one
    options = cumsum(options,2);
    %calculate delta t
    dt = log(1./rand(length(op_m(:,1)),1)).*(1./options(:,end));
    %if any cells would progress past end_t, remove them from simulation
    finished_cells = [finished_cells; op_m(op_m(:,end)+dt>end_t,:)];
    op_m(:,end) = op_m(:,end) + dt;
    options(op_m(:,end)>end_t,:) = [];
    op_m(op_m(:,end)>end_t,:) = [];
    %generate random numbers; these will dictate which reaction each cell
    %chooses from its options
    w = rand(length(op_m(:,1)),1);
    %translate random number into a value corresponding to sum of
    %propensities
    choices = w.*(options(:,end));
    %find which reaction corresponds to each choice for each cell
    rxn = zeros(length(op_m(:,1)),1);
    for i=1:length(op_m(:,1))
        rxn(i) = find(choices(i)<options(i,:),1);
    end
    %update each cell according to the stoichiometry of the chosen reaction
    op_m = op_m + stoich(rxn,:);
    
    %for cells that divide, we must randomly assign protein values to each
    %daughter cell
    % initialize
        first_cell_n_protein = op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:)));
        % for each cell
        for i=1:length(first_cell_n_protein(:,1))
            % for each protein
            for j=1:length(first_cell_n_protein(1,:))
                % assign a number of proteins to the first cell
                % binomial probability; 50% chance of each cell getting
                % each protein
                first_cell_n_protein(i,j) = sum(rand(first_cell_n_protein(i,j),1)>0.5);
            end
        end
        % assign remaining to second cell
        second_cell_n_protein = op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:)))-first_cell_n_protein;
    % for each protein that doesn't obey inheritence, revert all the
    % phospho-protein to unphosphorylated
    for i=1:length(uninherit_rxn)
        first_cell_n_protein(uninherit_rxn(i,:)==1,:) = first_cell_n_protein(uninherit_rxn(i,:)==1,:)+first_cell_n_protein(uninherit_rxn(i,:)==-1,:);
        second_cell_n_protein(uninherit_rxn(i,:)==1,:) = second_cell_n_protein(uninherit_rxn(i,:)==1,:)+second_cell_n_protein(uninherit_rxn(i,:)==-1,:);
        first_cell_n_protein(uninherit_rxn(i,:)==-1,:) = 0;
        second_cell_n_protein(uninherit_rxn(i,:)==-1,:) = 0;
    end
    % add the daughter cells to op_m
    op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:))) = first_cell_n_protein;
    op_m = [op_m; second_cell_n_protein op_m(rxn==length(stoich(:,1)),length(rxn_stoich(1,:))+1:end)];

end

end