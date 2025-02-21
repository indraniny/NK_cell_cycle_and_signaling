
n_states = 4; % number of cell cycle stages
n_proteins = 2; % number of proteins
all_mean_in_state = []; % Empty array to store all mean_in_state results (40 x 2)

load('A_levels_ode.mat', 'A_levels_ode');
load('A_phosph_levels_ode.mat', 'A_phosph_levels_ode');
load('times.mat','times');
%times=0.05:0.3:50;
%times=0.05:0.3:0.6;
length(times)
% Loop over cells_at_t from 1 to 10
for k = 1:length(times)
    % Call the function with current cells_at_t(t)
    [proteins_dist_state, mean_in_state, protein_levels, cell_levels] = parse_cell_values_sigma(cells_at_t(k), n_states, n_proteins);

    % Append mean_in_state (4x2) to all_mean_in_state
    all_mean_in_state = [all_mean_in_state; mean_in_state]; % concatenate 4x2 matrices
end

% Save the result to a .mat file
save('mean_in_state_results.mat', 'all_mean_in_state');

max_length=(length(times)*4)
%First column is A and second column is A*
protein_G1_levels=all_mean_in_state(1:4:max_length,:)
protein_S_levels=all_mean_in_state(2:4:max_length,:)
protein_G2_levels=all_mean_in_state(3:4:max_length,:)
protein_M_levels=all_mean_in_state(4:4:max_length,:)

length(times)
save('Gillespie_protein_G1_levels.mat','protein_G1_levels')
save('Gillespie_protein_S_levels.mat','protein_S_levels')
save('Gillespie_protein_G2_levels.mat','protein_G2_levels')
save('Gillespie_protein_M_levels.mat','protein_M_levels')



%plotting protein A with time
figure()
plot(times, protein_G1_levels(:,1), 'DisplayName', 'Protein G1')
hold on
plot(times, protein_S_levels(:,1), 'DisplayName', 'Protein S')
plot(times, protein_G2_levels(:,1), 'DisplayName', 'Protein G2')
plot(times, protein_M_levels(:,1), 'DisplayName', 'Protein M')

% Plot ODE protein levels with increased line thickness and dotted lines
plot(times, A_levels_ode(1,:),'b-', 'DisplayName', 'Protein G1-ODE', 'LineStyle', ':', 'LineWidth', 2)
hold on
plot(times, A_levels_ode(2,:), 'DisplayName', 'Protein S-ODE', 'LineStyle', ':', 'LineWidth', 2)
plot(times, A_levels_ode(3,:), 'DisplayName', 'Protein G2-ODE', 'LineStyle', ':', 'LineWidth', 2)
plot(times, A_levels_ode(4,:), 'DisplayName', 'Protein M-ODE', 'LineStyle', ':', 'LineWidth', 2)

% Add a title
title('Gillespie: a Across Cell Cycle Stages')
% Add a legend
legend('show')
% Optionally, label the axes
xlabel('Time')
ylabel('Avg. Protein Level')








%plotting protein A* with time
figure()
plot(times, protein_G1_levels(:,2),'DisplayName','Protein G1', 'LineWidth', 2)
hold on
plot(times, protein_S_levels(:,2), 'DisplayName', 'Protein S', 'LineWidth', 2)
plot(times, protein_G2_levels(:,2), 'DisplayName', 'Protein G2', 'LineWidth', 2)
plot(times, protein_M_levels(:,2), 'DisplayName', 'Protein M', 'LineWidth', 2)

hold on;
lw=6
% Plot ODE protein levels with increased line thickness and dotted lines
plot(times, A_phosph_levels_ode(1,:),'b-', 'DisplayName', 'Protein G1-ODE', 'LineStyle', ':', 'LineWidth', lw)
hold on
plot(times, A_phosph_levels_ode(2,:), 'DisplayName', 'Protein S-ODE', 'LineStyle', ':', 'LineWidth', lw)
plot(times, A_phosph_levels_ode(3,:), 'DisplayName', 'Protein G2-ODE', 'LineStyle', ':', 'LineWidth', lw)
plot(times, A_phosph_levels_ode(4,:), 'DisplayName', 'Protein M-ODE', 'LineStyle', ':', 'LineWidth', lw)

% Add a title
title('Gillespie vs. ODE: a* Across Cell Cycle Stages')
% Add a legend
legend('show')
% Optionally, label the axes
xlabel('Time')
ylabel('Avg. Protein Level')


