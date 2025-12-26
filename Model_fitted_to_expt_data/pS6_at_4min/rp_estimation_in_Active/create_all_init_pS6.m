function [all_init]=create_all_init_pAkt(number_of_workers,factor,n_select)

all_init={}; %creating an empty cell array

% Load initial configuration for the cells from a data file at t=0 min


% %Load experimental data at 2 min to input into the model
data_1 = readmatrix('pS6_expression_2min_G1.csv');
data_2 = readmatrix('pS6_expression_2min_S.csv');
data_3 = readmatrix('pS6_expression_2min_G2.csv');


%data preprocessing on input data at 2 min %Dec 14, 2025
data_1(data_1 < 0) = 1e-11; %G1
data_2(data_2 < 0) = 1e-11; %S
data_3(data_3 < 0)= 1e-11; %M


data_G1 = [data_1(2:end,2), ones(length(data_1)-1,1), zeros(length(data_1)-1,2)];
data_S = [data_2(2:end,2), zeros(length(data_2)-1,1), ones(length(data_2)-1,1), zeros(length(data_2)-1,1)];
data_G2=[data_3(2:end,2), zeros(length(data_3)-1,2),ones(length(data_3)-1,1)];

tot_cell=length(data_G1)+length(data_S)+length(data_G2);

cell_config_0=zeros(tot_cell,9); 
cell_config_0(:,4:7)=[data_G1;data_S;data_G2];



cell_config_0(:,4) = floor(cell_config_0(:,4) * factor); % Scale the pAkt protein fluroscence to protein number


total_cells=length(cell_config_0);

G1_tot=sum(cell_config_0(:,5))
S_tot=sum(cell_config_0(:,6))
G2_tot=sum(cell_config_0(:,7))
M_tot=sum(cell_config_0(:,8))

(G1_tot/total_cells)*100
(S_tot/total_cells)*100
(G2_tot/total_cells)*100
(M_tot/total_cells)*100


sigma = 0.05;                  % log-space standard deviation

muG1=300
muS=muG1*1.2065
muG2=muG1*1.3045
muM=muG1*1.3351

% For each worker assign the cell distribution [E S C P], where E, S, C
% follows lognormal distribution with mu and sigma, P (here, pS6) is
% collected from the experimental data at 2 min
all_init = cell(number_of_workers, 1);

%---------------------------------------------------------------
for i = 1:number_of_workers
    % Randomly choose n_select rows from experimental data
    temp= randperm(total_cells, n_select);  % <== keep this line as you had
    init = zeros(n_select, size(cell_config_0,2));

    % Assign experimental data for fluorescence and stage flags
    init(:,4)     = cell_config_0(temp, 4);      % pS6 fluorescence
    init(:,5:end) = cell_config_0(temp, 5:end);  % stage flags (G1, S, G2, M)


    % Initialize synthetic variables for E, S, C (columns 1,3,4)
    for k = 1:n_select
        if init(k,5) == 1          % G1
            init(k,1) = floor(lognrnd(log(muG1), sigma));
            init(k,2) = floor(lognrnd(log(muG1), sigma));
            init(k,3) = floor(lognrnd(log(muG1), sigma));
        elseif init(k,6) == 1      % S
            init(k,1) = floor(lognrnd(log(muS), sigma));
            init(k,2) = floor(lognrnd(log(muS), sigma));
            init(k,3) = floor(lognrnd(log(muS), sigma));
        elseif init(k,7) == 1      % G2
            init(k,1) = floor(lognrnd(log(muG2), sigma));
            init(k,2) = floor(lognrnd(log(muG2), sigma));
            init(k,3) = floor(lognrnd(log(muG2), sigma));
        elseif init(k,8) == 1      % M
            init(k,1) = floor(lognrnd(log(muM), sigma));
            init(k,2) = floor(lognrnd(log(muM), sigma));
            init(k,3) = floor(lognrnd(log(muM), sigma));
        end
    end

    % Store initialization for this worker
    all_init{i} = init;
end

end