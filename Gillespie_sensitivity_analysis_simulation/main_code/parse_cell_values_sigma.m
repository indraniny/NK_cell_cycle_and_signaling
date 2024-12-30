function [proteins_dist_state,mean_in_state,std_in_state,protein_levels,cell_levels] = parse_cell_values_sigma(cells,n_states,n_proteins)

% Author: Darren Wethington
% This function parses numbers of protein and cells in each cell cycle
% stage.
% Inputs: cells: the output from the KMC script
            % n_states: # cell cycle stages
            % n_proteins: # species
% Output protein_levels is a matrix where rows correspond to cell cycle
% stages and columns correspond to species

% % DW Edit 5/30/2024
% if nargin < 4

% reformat the input data if there are many cells (data format) of data
n_datasets = length(cells)
rearranged_cells = cells{1}
for i=2:n_datasets
    rearranged_cells = [rearranged_cells; cells{i}];
end
% count the number of proteins corresponding to each stage of the cell
% cycle
for i=1:n_states % for each stage
    protein_levels(i,:) = sum(rearranged_cells(rearranged_cells(:,end-n_states-1+i)==1,1:n_proteins),1)/n_datasets;
    cell_levels(i) = sum(rearranged_cells(:,end-n_states-1+i))/n_datasets;
    % I divide by n_datasets assuming we are looking to average across the
    % noise, not sum across them
end

% %example DW code 5/30/2024
% 
protein_levels
cell_levels

for i=1:n_states
    protein_levels_in_state = rearranged_cells(rearranged_cells(:,end-n_states-1+i)==1,1:n_proteins);
    proteins_dist_state{i}= protein_levels_in_state;
    mean_in_state(i,:) = mean(protein_levels_in_state,1);
    std_in_state(i,:) = std(protein_levels_in_state,1);
%     mean_in_state = mean(rearranged_cells,1)
%     std_in_state = std(rearranged_cells,1)
end

 proteins_dist_state{4}
 proteins_dist_state{4}(:,2)
 mean_in_state(4,:)
% class(proteins_dist_state)

% Each row represents E S C P respectively and column represents cell cycle
% stages

mean_in_state
std_in_state

%Date Nov 11, 2024, I found that many cases there is no cell with proteins
%in some stages, so, cell number in denominator is giving NAN and thus we
%replaces average protein is 0 in that stage

%mean_in_state(isnan(mean_in_state)) = 0

end