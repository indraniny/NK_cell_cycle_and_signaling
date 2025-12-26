    
function [X1_all,X2_all,X3_all,X4_all]=process_protein_averages(number_of_workers,qq)

    num_files = number_of_workers;  % Adjust this to the number of files you have
    n_stages=4;
    n_proteins=4;

    % Initialize an empty array to store the results
    P = [];  % Matrix to store p1, p2, p3 as columns
    all_mean=[];
    all_cell_levels=[];
    %{E S ES S*]
    % Loop over the files


    % Initialize empty matrices before the loop
    X1_all = [];
    X2_all = [];
    X3_all = [];
    X4_all = [];

    for i = 1:num_files
        % Construct the filename for the current iteration
        filename = sprintf('init_cond_all_rate_scan_%d.mat', i);
        load(filename,'cells_at_t');
        cell_matrix=cells_at_t(1);
     
        [X1,X2,X3,X4] = get_protein_dist(cell_matrix,n_stages,n_proteins,qq);

        %protein distribution for protein index qq (e.g., qq=4 represents
        %phospho-protein"

        X1_all = [X1_all; X1];  % G1
        X2_all = [X2_all; X2];  % S
        X3_all = [X3_all; X3];  % G2
        X4_all = [X4_all; X4];  % M

        % disp('-------For each worker protein dist-------')
        % length(X1)
        % length(X2)
        % length(X3)
        % length(X4)
        
    end

    
    % disp('show the protein distribution combined all workers for ')
    % qq
    % mean(X1_all)
    % mean(X2_all)
    % mean(X3_all)
    % mean(X4_all)   
end