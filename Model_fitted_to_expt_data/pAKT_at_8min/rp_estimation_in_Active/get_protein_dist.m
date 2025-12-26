function [X1,X2,X3,X4] = get_protein_dist(cell_matrix,n_stages,n_proteins,qq)

%n_stages=4 [G1,S,G2, M]
%n_proteins=4 [E S C P]
%qq =Protein index to see statistics t-test



%disp( '[E S C P]');

close all;

%This function calculates the number of cells in each stages and their
%number of E,S,C,P. To plot protein distribution one has to play with
%proteins_dist_state, mean, std. Along Y axis is cell cycle stages and
%along X axis is the protein [E S C P]
[proteins_dist_state, mean_in_state, std_in_state, protein_levels,cell_levels] = parse_cell_values_sigma(cell_matrix,n_stages,n_proteins);


%disp(protein_levels);
total_protein_levels=protein_levels; 
% [E S C P]
% [E vav E-vav pvav]
%Here, E S indicates free Enzyme and Complex
% Total Enzyme = Etot=E+C
%Total substrate= Stot=S+C+P



%E, S, C , P
%represents the index for protein, For e.g., qq=4 is the phospho-protein or pVav1 
X1=proteins_dist_state{1}(:,qq); %G1 Use this data to plot the prob dist
X2=proteins_dist_state{2}(:,qq); %S Use this data to plot the prob dist
X3=proteins_dist_state{3}(:,qq); %G2 Use this data to plot the prob dist
X4=proteins_dist_state{4}(:,qq);%M Use this data to plot the prob dist
% %------------------------------------------------------------
% %Date: Oct 6, 2025 to check the total protein
% disp('S')
% proteins_dist_state{1}(:,4)
% proteins_dist_state{1}(:,3)
% proteins_dist_state{1}(:,2)
% %To check total substrate Stot=S+Complex+S*
% disp('S_tot')
% X1=proteins_dist_state{1}(:,4)+proteins_dist_state{1}(:,3)+proteins_dist_state{1}(:,2) %G1 Use this data to plot the prob dist
% X2=proteins_dist_state{2}(:,4)+proteins_dist_state{2}(:,3)+proteins_dist_state{2}(:,2); %S Use this data to plot the prob dist
% X3=proteins_dist_state{3}(:,4)+proteins_dist_state{3}(:,3)+proteins_dist_state{3}(:,2); %G2 Use this data to plot the prob dist
% X4=proteins_dist_state{4}(:,4)+proteins_dist_state{4}(:,3)+proteins_dist_state{4}(:,2);%M Use this data to plot the prob dist

end