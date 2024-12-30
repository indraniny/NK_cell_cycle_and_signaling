
function [data]=process_protein_averages_for_sensitivity()
    % Define the number of files you want to process
    num_files = 20;  % Adjust this to the number of files you have
    n_stages=4
    n_proteins=4
    qq=4 %which protein you want to check mean or significance
    
    % Initialize an empty array to store the results
    P = [];  % Matrix to store p1, p2, p3 as columns
    all_mean=[]
    
    % Loop over the files
    for i = 1:num_files
        % Construct the filename for the current iteration
        filename = sprintf('init_cond_all_rate_scan_%d.mat', i);
        load(filename,'cells_at_t')
        cell_matrix=cells_at_t(1)
        [mean_in_state, std_in_state,protein_levels, cell_levels, avgs,total_avgs,p1,p2,p3] = get_protein_averages(cell_matrix,n_stages,n_proteins,qq);
        protein=mean_in_state(:,qq)
        all_mean=[all_mean,protein]
                
        % Append p1, p2, and p3 as a new row in matrix P
        P = [P; p1, p2, p3];  % Each row is p1, p2, p3 for one iteration
    end
    
    % Save the collected data into a single .mat file
    save('processed_protein_averages.mat', 'P');
    format shortE
    P
    all_mean
    
    avg_mean=mean(all_mean,2)
    std_of_mean=std(all_mean,0,2)
    
    % Adjust font size and axis width
    ax = gca;
    ax.FontSize = 22;
    ax.LineWidth = 2;
    xlim([0.8 4.2]);  % Set x-axis range from 0.5 to 4
    x_values = [1, 2, 3, 4];

    % Create the plot with error bars
    figure;  % Create a new figure
    errorbar(x_values, avg_mean, std_of_mean, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');

    xticks([1 2 3 4]);  % Set tick positions corresponding to data points
    xticklabels({'G1', 'S', 'G2', 'M'});  % Manually set labels
    legend( '<mean>','FontSize', 16, 'LineWidth', 2);
    
    combined= [avg_mean,std_of_mean]
%     writematrix(combined, 'avg_and_std.txt');
%     
%     writematrix(all_mean, 'all_mean.txt'); %Column represents an independent run
%     
%     writematrix(P, 'p_values.txt'); %row represents an independent run
    
    % Assign each stage to a variable
    G1 = all_mean(1, :);
    S = all_mean(2, :);
    G2 = all_mean(3, :);
    M = all_mean(4, :);
    
     %Perform Welch's t-test (unequal variances) between the specified stages
    % G1 vs S
    [h1, p1, ci1, stats1] = ttest2(G1, S, 'Vartype', 'unequal');

    % S vs G2
    [h2, p2, ci2, stats2] = ttest2(S, G2, 'Vartype', 'unequal');

    % G2 vs M
    [h3, p3, ci3, stats3] = ttest2(G2, M, 'Vartype', 'unequal');

    fprintf('T-test results between G1 and S:\n');
    fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats1.tstat, p1);

    fprintf('T-test results between S and G2:\n');
    fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats2.tstat, p2);

    fprintf('T-test results between G2 and M:\n');
    fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats3.tstat, p3);

    avg_meant=transpose(avg_mean)
    std_of_meant=transpose(std_of_mean)
    data=[avg_meant,std_of_meant,p1,p2,p3]
    
end% Save to a .mat file

