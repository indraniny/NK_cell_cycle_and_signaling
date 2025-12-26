%This code is to estimate parameters for the neutral scenario, where
%phosphorylation rate of phospho-protein(pAkt) (rp) is cell cycle-independent (kcat=rp)
%Date: Dec 14, 2025 (extra functions deleted, return variables in the
%functions are shortened)


% Clear all variables, globals, functions, and MEX links
delete all;
all_init={}; %creating an empty cell array
number_of_workers = 20;% Number of parallel workers to use for simulations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Experimental data%%%%%%%%%%%%%%%%%%%%
factor = 100.0; % Factor to scale fluroscence to protein number in the experiments
N=3 %Number of data points to compare (G1, S, G2-M)
K=1 %number of estimated parameters, here 1 (constant rp for neutral case)
n_select=20; % randomly selecting only n rows to reduce the computational runtime
initial_cell_no=n_select;

%load Experimental data at 2 min, distributed between parallel workers 
[all_init]=create_all_init_pAkt(number_of_workers,factor,n_select);


%celldisp(all_init);


% %Load experimental data at 8 min to be fitted
data_1 = readmatrix('pAkt_expression_8min_G1.csv');
data_2 = readmatrix('pAkt_expression_8min_S.csv');
data_3 = readmatrix('pAkt_expression_8min_G2.csv');

%negative data preprocessing on data at 8 min 
data_1(data_1 < 0) = 1e-11;
data_2(data_2 < 0) = 1e-11;
data_3(data_3 < 0)= 1e-11;

%multiplying the experimental data with a factor to convert it from fluroscence to
%protein number and then taking logarithm (here we took base 10) of it
log_data_G1=log((data_1(2:end,2)*factor)); 
log_data_S=log((data_2(2:end,2)*factor));
log_data_G2=log((data_3(2:end,2)*factor));



% These log values are calculated to minimize the cost function 
disp('data to be fitted at the later time point')
X = {log_data_G1, log_data_S, log_data_G2};

mean_expt = cellfun(@mean, X)
se_expt = cellfun(@(x) std(x)/sqrt(numel(x)), X)
var_expt  = cellfun(@var, X)  % variance directly



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimating rp by in silico model optimization %%%%%%%%%%%%%%%%%%%%
%------unblock the bottom part while doing RSS optimization%-------------
% t_end = 0.1; % 6 min/ 60  h
% % Define initial guess
% lossfun = @(params) loss_function(params, t_end, initial_cell_no, number_of_workers, all_init, mean_expt, var_expt,factor)
% 
% 
% 
% % Initial guess for parameters (kcat values for G1, S, G2M)
% x0 =[600.0];
% 
% % Lower and upper bounds
% lb = [350.0];
% ub = [2000.00];
% 
% %% ----------------------------------------------------
% % Stage 1: Global search with patternsearch
% % ----------------------------------------------------
% options_ps = optimoptions('patternsearch', ...
%     'Display', 'iter', ...        % show progress
%     'MeshTolerance', 1e-4, ...    % mesh refinement stopping criterion
%     'StepTolerance', 1e-4, ...    % convergence stopping criterion
%     'UseParallel', false);         % use parallel workers if available
% 
% fprintf('--- Running patternsearch (global search) ---\n');
% [param_ps, fval_ps] = patternsearch(lossfun, x0, [], [], [], [], lb, ub, options_ps);
% 
% %% ----------------------------------------------------
% % Stage 2: Local refinement with fmincon
% % ----------------------------------------------------
% options_fmin = optimoptions('fmincon', ...
%     'Display', 'iter', ...
%     'Algorithm', 'sqp', ...
%     'FiniteDifferenceStepSize', 0.1, ...  % make steps bigger
%     'StepTolerance', 1e-6, ...
%     'OptimalityTolerance', 1e-6, ...
%     'MaxFunctionEvaluations', 5000);
% 
% fprintf('--- Refining with fmincon (local search) ---\n');
% [param_est, fval] = fmincon(lossfun, param_ps, [], [], [], [], lb, ub, [], options_fmin);
% 
% %% ----------------------------------------------------
% % Results
% % ----------------------------------------------------
% fprintf('\nFinal parameter estimates:\n');
% disp(param_est);
% 
% best_kcat   = param_est(1);
% 
% 
% 
% fprintf('Best parameters:\n kcat = %.4g\n', best_kcat);
% 
% 
% 
% %plot for optimized parameters 
% kcat =best_kcat; % Set the catalytic constant (reaction rate)[36-7200000]
% 
% 
% %% Save to a text file
% fileID = fopen('optimized_parameters_neutral.txt','w');
% fprintf(fileID, 'Optimized parameters:\n');
% fprintf(fileID, 'kcat   = %.6f\n', kcat);
% fclose(fileID);
% ----------------------------------------------------------------------


% % Set the simulation end time
t_end =0.1 ;% 6 min/ 60  h
qq=4 %protein index [E S C P] here qq=4 represents phospho-protein
%-----block the bottom part when optimizing parameter, open for plotting -----
%estimated parameters
kcat = 1464.00 %optimized parameter
% %------------------------------------------------------------------------


% Run the main simulation model with parameters, calculated for all
% parallel workers
NL_model_configuration_generator(t_end, initial_cell_no, number_of_workers, kcat,all_init);

%Process the simulation output to compute protein averages, extract [E S C
%P] protein distribution for each workers and concatenate them
[X1_all,X2_all,X3_all,X4_all] = process_protein_averages(number_of_workers,qq);

%If model gives the protein abundances to be 0
X1_all(X1_all==0) = 1e-11*factor %small rate or factor can cause 0 phospho molecule
X2_all(X2_all==0) = 1e-11*factor 
X3_all(X3_all==0) = 1e-11*factor 
X4_all(X4_all==0) = 1e-11*factor 


%take log of simulation values to match experiments
X_G1=log(X1_all); %proteins at G1
X_S=log(X2_all); %proteins at S

if all(isnan(X4_all))
    X_G2M = log(X3_all); % use only G2 if M is all NaN
    disp('if')
else
    X_G2M = [log(X3_all); log(X4_all)];
    disp('else')
end


%statistics of log(copy numbers) in model
X2 = {X_G1, X_S, X_G2M}
mean_model=cellfun(@mean,X2)
std_model = cellfun(@std, X2)
se_model  = cellfun(@(x) std(x)/sqrt(numel(x)), X2)

%%%%%%%%%%%%%%%%%%%Finding loss%%%%%%%%%%%%%%%%%%%%%%%%

% Compute element-wise squared differences or loss for optimizing
RSS = sum(((mean_expt - mean_model).^2)./ var_expt)


%AIC score calculation based on loss
AIC=N*log(RSS/N)+(2*K)+((2*K*(K+1))/(N-K-1))


%plotting
% Create bar plot comparing experimental vs simulated means with error bars
% Category labels for the first 3 proteins
categories = {'G1', 'S', 'G2-M'};



% Grouped bar plot data: 2x3 (experimental vs simulated)
data = [mean_expt; mean_model];
%data = [mean_model];

% Plot grouped bars
figure();
hold on;
b = bar(data', 'grouped');
b(1).FaceColor = [0.2, 0.6, 0.8];  % Light blue for experimental
b(2).FaceColor = [0.8, 0.4, 0.4];  % Light red for simulated


x(1,:) = [0.9, 1.9, 2.9];   % left positions
x(2,:) = [1.1, 2.1, 3.1];
numGroups = size(data, 2); % Should be 3 (G1-G2)


% Add error bars
errorbar(x(1,:), mean_expt, se_expt, ...
    'k.', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(x(2,:), mean_model, se_model, ...
    'k.', 'LineWidth', 1.5, 'CapSize', 10);

% Customize plot
xticks(1:numGroups);
xticklabels(categories);
%ylabel('log(Protein abundances)');
legend({'Experiment', 'Model'}, 'FontSize', 28, 'Location', 'NorthWest');
set(gca, 'FontName', 'Arial', 'FontSize', 36, 'LineWidth', 2);
hold off;



disp('-------Model------------')
%Perform Welch's t-test (unequal variances) between the specified stages
%two-tailed t-test whether the means of two independent samples are significantly different.
% G1 vs S
[h1, p1, ci1, stats1] = ttest2(X_G1, X_S, 'Vartype', 'unequal');

% S vs G2
[h2, p2, ci2, stats2] = ttest2(X_S, X_G2M, 'Vartype', 'unequal');


fprintf('T-test results between G1 and S:\n')
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats1.tstat, p1)

fprintf('T-test results between S and G2-M:\n');
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats2.tstat, p2)


disp('-------Expt. data (two-tailed t test)------------')
%Perform Welch's t-test (unequal variances) between the specified stages
%two-tailed t-test whether the means of two independent samples are significantly different.
% G1 vs S
[h1, p1, ci1, stats1] = ttest2(log_data_G1, log_data_S, 'Vartype', 'unequal');

% S vs G2
[h2, p2, ci2, stats2] = ttest2(log_data_S, log_data_G2, 'Vartype', 'unequal');


fprintf('T-test results between G1 and S:\n')
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats1.tstat, p1)

fprintf('T-test results between S and G2-M:\n');
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats2.tstat, p2)



disp('-------Expt. data (one tailed t test)------------')%p values will
%be half of the two-tailed t-test
%Perform Welch's t-test (unequal variances) between the specified stages
%two-tailed t-test whether the means of two independent samples are significantly different.
% G1 vs S
[h1, p1, ci1, stats1] = ttest2(log_data_G1, log_data_S, 'Vartype', 'unequal','Tail', 'left');

% S vs G2
[h2, p2, ci2, stats2] = ttest2(log_data_S, log_data_G2, 'Vartype', 'unequal','Tail', 'left');


fprintf('T-test results between G1 and S:\n')
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats1.tstat, min(p1, 1-p1))

fprintf('T-test results between S and G2-M:\n');
fprintf('T-statistic: %.4e, P-value: %.4e\n\n', stats2.tstat, min(p2, 1-p2))




function RSS = loss_function(params, t_end, initial_cell_no, number_of_workers, all_init, mean_expt, var_expt,factor)

    % -------------------- Unpack parameters --------------------
    kcat    = params(1);


    qq = 4; % protein index (phospho-protein)

    % -------------------- Run simulation --------------------
    NL_model_configuration_generator(t_end, initial_cell_no, number_of_workers, kcat,all_init);
    % -------------------- Process outputs --------------------
    
    [X1_all,X2_all,X3_all,X4_all] = process_protein_averages(number_of_workers,qq);
    
    %If model gives the protein abundances to be 0
    X1_all(X1_all==0) = 1e-11*factor %small rate or factor can cause 0 phospho molecule
    X2_all(X2_all==0) = 1e-11*factor 
    X3_all(X3_all==0) = 1e-11*factor 
    X4_all(X4_all==0) = 1e-11*factor 

    % -------------------- Log-transform safely --------------------
    X_G1  = log(X1_all);
    X_S   = log(X2_all);


    if all(isnan(X4_all))
        X_G2M = log(X3_all); % use only G2 if M is all NaN
    else
        X_G2M = [log(X3_all); log(X4_all)];
    end



    mean_model = [mean(X_G1), mean(X_S), mean(X_G2M)];

    % -------------------- Compute weighted squared error --------------------
    % Compute element-wise squared differences or loss for optimizing
    RSS = sum(((mean_expt - mean_model).^2)./ var_expt)

    % -------------------- Optional: log for debugging --------------------
    fid = fopen('loss_track_neutral.txt', 'a');  % append mode
    fprintf(fid, 'kcat=%.4f,  RSS=%.6f, mean_model=[%.3f %.3f %.3f]\n', ...
        kcat, RSS, mean_model(1), mean_model(2), mean_model(3));
    fclose(fid);
    

end
