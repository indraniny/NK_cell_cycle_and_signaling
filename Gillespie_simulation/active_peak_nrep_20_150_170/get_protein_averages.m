function [mean_in_state, std_in_state,protein_levels, cell_levels, avgs,total_avgs,p1,p2,p3] = get_protein_averages(cell_matrix,n_stages,n_proteins,qq)
%cell_matrix=cells_at_t(1) higher time or, cells_at_t(2) lower time 
%n_stages=4 [G1,S,G2, M]
%n_proteins=4 [E S C P]
%qq =Protein index to see statistics t-test, histogram



disp( '[E S C P]')

close all;

%This function calculates the number of cells in each stages and their
%number of E,S,C,P. To plot protein distribution one has to play with
%proteins_dist_state, mean, std. Along Y axis is cell cycle stages and
%along X axis is the protein [E S C P]
[proteins_dist_state, mean_in_state, std_in_state, protein_levels,cell_levels] = parse_cell_values_sigma(cell_matrix,n_stages,n_proteins)


disp(protein_levels)
total_protein_levels=protein_levels; 
% [E S C P]
% [E vav E-vav pvav]
%Here, E S indicates free Enzyme and Complex
% Total Enzyme = Etot=E+C
%Total substrate= Stot=S+C+P

%Unblock it if you want to see total E= Free E + C, total S= Free S+C + 
% Total Enzyme = Etot=E+C
%Total substrate= Stot=S+C+P
%to plot the total E, S plot from total protein level

total_protein_levels(:,1)=total_protein_levels(:,1)+total_protein_levels(:,3);
total_protein_levels(:,2)=total_protein_levels(:,2)+total_protein_levels(:,3)+total_protein_levels(:,4); %p-vav1 is added on Aug 31, 2024

disp('total_protein_levels=')

disp( '[E S C P]')
avgs = protein_levels./cell_levels' %This is for [free E, free S, C, P]
diff1=avgs(2,:)-avgs(1,:)
disp(diff1)
total_avgs = total_protein_levels./cell_levels' %This is [total E, total S, C, P]
% columns will correspond to a protein's progression through the cell
% cycle.  Each column represents a different protein, while each row
% represents a different cell cycle stage.
% x=['G1','S','G2','M']
% figure()
% plot(x, avgs(:,1), 'r*--', 'LineWidth', 2);
% hold on
% plot(x, avgs(:,2), 'b*--', 'LineWidth', 2);
% hold on
% plot(x, total_avgs(:,1), 'ro-', 'LineWidth', 2);
% hold on
% plot(x, total_avgs(:,2), 'bo-', 'LineWidth', 2);
% hold on
% plot(x, avgs(:,3), 'm*--', 'LineWidth', 2);
% hold on
% plot(x, avgs(:,4), 'g*--', 'LineWidth', 2);
% hold on
% ax = gca; % Get current axes
% ax.FontSize = 22; % Adjust font size of axes tick labels
% ax.LineWidth = 2; % Adjust axis width


x = [1, 2, 3, 4]; % Numeric values for data points



% ... rest of your plotting code (refer to previous enhanced version)


figure('DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);  % Set default font for all axes in the figure
dark_green = [0, 0.5, 0]; % RGB values for dark green
% Plot lines with bigger point size
plot(x, avgs(:,1), 'r*:', 'MarkerSize', 12, 'LineWidth', 2);%Free E
hold on
plot(x, avgs(:,2), 'b*:', 'MarkerSize', 12, 'LineWidth', 2); %Free S
hold on
plot(x, total_avgs(:,1), 'rs--', 'MarkerSize', 12, 'LineWidth', 2);%Total E
hold on
plot(x, total_avgs(:,2), 'bs--', 'MarkerSize', 12, 'LineWidth', 2);%Total S
hold on
plot(x, avgs(:,3), 'm^:', 'MarkerSize', 12, 'LineWidth', 2);%Complex
hold on
plot(x, avgs(:,4), 'gO:', 'MarkerSize', 12, 'LineWidth', 2);  % pVav1- Green dotted line



% Add labels and title
% xlabel('Cell Cycle Phase');
% ylabel('Average Values');
% title('Plot of Average Values');

% Adjust font size and axis width
ax = gca;
ax.FontSize = 22;
ax.LineWidth = 2;
xlim([0.8 4.2]);  % Set x-axis range from 0.5 to 4
% Additional formatting options (optional)
% grid on;  % Add grid lines
% Set custom x-axis labels
xticks([1 2 3 4]);  % Set tick positions corresponding to data points
xticklabels({'G1', 'S', 'G2', 'M'});  % Manually set labels
legend('free <E>', ' free <S>','total <E>', ' total <S>', '<C>', '<P>','FontSize', 16, 'LineWidth', 2);
%legend('<C>', '<P>','FontSize', 22, 'LineWidth', 2);

hold on


figure('DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);  % Set default font for all axes in the figure
% Plot lines with bigger point size
plot(x, avgs(:,3), 'm^:', 'MarkerSize', 12, 'LineWidth', 2);


figure('DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);
plot(x, avgs(:,4), 'go:', 'MarkerSize', 12, 'LineWidth', 2);

% Adjust font size and axis width
ax = gca;
ax.FontSize = 22;
ax.LineWidth = 2;
xlim([0.8 4.2]);  % Set x-axis range from 0.5 to 4
% Additional formatting options (optional)
% grid on;  % Add grid lines
% Set custom x-axis labels
xticks([1 2 3 4]);  % Set tick positions corresponding to data points
xticklabels({'G1', 'S', 'G2', 'M'});  % Manually set labels
legend( '<C>', '<P>','FontSize', 16, 'LineWidth', 2);
%legend('<C>', '<P>','FontSize', 22, 'LineWidth', 2);

diff1=avgs(2,:)-avgs(1,:)
diff2=avgs(3,:)-avgs(2,:)
diff3=avgs(4,:)-avgs(3,:)

%E, S, C , P
%represents the index for protein, For e.g., qq=4 is the phospho-protein or pVav1 
X1=proteins_dist_state{1}(:,qq); %G1 Use this data to plot the prob dist
X2=proteins_dist_state{2}(:,qq); %S Use this data to plot the prob dist
X3=proteins_dist_state{3}(:,qq); %G2 Use this data to plot the prob dist
X4=proteins_dist_state{4}(:,qq)%M Use this data to plot the prob dist

%plotting histogram for protein [E S C P] with qq index at G1, S, G2 and M stages
% figure()
% bin=20
% histogram(X1,bin,FaceAlpha=0.6, FaceColor='blue')%G1
% hold on;
% histogram(X2,bin, FaceAlpha=0.6,FaceColor='yellow')%S
% hold on;
% hist(X3,bin,FaceAlpha=0.6,FaceColor='red')%G2
% hold on;
% hist(X4,bin,FaceAlpha=0.6,FaceColor='maroon')%M
% hold on;

% Protein trend
% t-test 1: Compare mean between G1 and S
disp('-------*********-------')
disp('compare G1 and S')
[h1,p1] = ttest2(X1,X2,'Vartype','unequal','Alpha',0.05)

% t-test 2: Compare mean between S and G2
disp('-------*********-------')
disp('compare S and G2')
[h2,p2] = ttest2(X2,X3,'Vartype','unequal','Alpha',0.05)
% t-test 3: Compare mean between G2 and M
disp('-------*********-------')
disp('compare G2 and M')
[h3,p3] = ttest2(X3,X4,'Vartype','unequal','Alpha',0.05)
disp('-------*********-------')
fprintf('Mean proteins:\n %f \n %f \n %f \n %f \n', mean(X1),mean(X2),mean(X3),mean(X4) );
fprintf('[Hypothesis p]:\n %f \t %e\n %f \t %e \n %f \t %e \n', h1, p1, h2, p2, h3, p3);
end