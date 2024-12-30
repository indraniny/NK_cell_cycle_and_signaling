function sol = convert_gillespie_setup_to_ode(stoich,rxn_reacts,C,rates,cycle,heritibility,end_t)

% Author: Darren Wethington
% Developer's note: I apologize for how complex this code is.  It's just a
% complex thing I'm trying to do.  Currently, the input C is not robust; it
% assumes there is only one cell in the input and that it's in the first
% cell cycle stage.  If this needs to change you will have to work on it.
% Also, heritibility is currently always assumed; if you want to program a
% case with no heritability it will impact the allocation of the protein
% transitions from the last cell cycle stage to the first.  I'll designate
% where this is in the code.

% The input stoich will be nxm, with n reactions and m species
% The output should be (nxk+k+m)x(mxk+k), to represent that each reaction happens
% in each cell cycle stage, that there is a measurement of protein in
% each stage, and to measure cell stage abundances.

% if anyone down the line wants to improve this, a quick and easy way for
% readability would be to replace all my size/length statements with a
% definition of n, k, and m at the top of the script.

% initialize the output stoichiometry matrix
out_stoich = zeros((size(stoich)+1)*length(cycle));
out_stoich = [out_stoich; zeros(length(stoich(1,:))*length(cycle),length(out_stoich(1,:)))];
% reactions happen in each cycle stage
for i=1:length(cycle)
    out_stoich((1:length(stoich(:,1)))+(i-1)*length(stoich(:,1)),(1:length(stoich(1,:)))+(i-1)*length(stoich(1,:))) = stoich;
end

% transitions between stages transfer protein between stages
temp = eye(length(stoich(1,:))*length(cycle));
out_stoich(length(stoich(:,1))*length(cycle)+(1:(length(stoich(1,:))*length(cycle))),1:length(temp(1,:))) = -temp;
temp(:,end-(length(stoich(1,:))-1:-1:0)) = [];
out_stoich(length(stoich(:,1))*length(cycle)+(1:(length(stoich(1,:))*length(cycle))),length(stoich(1,:))+(1:length(temp(1,:)))) = out_stoich(length(stoich(:,1))*length(cycle)+(1:(length(stoich(1,:))*length(cycle))),length(stoich(1,:))+(1:length(temp(1,:)))) + temp;
%this next line is where you would want to designate if heritibility were
%to not exist.
out_stoich(length(stoich(:,1))*length(cycle)+((length(stoich(1,:))*(length(cycle)-1)+1):(length(stoich(1,:))*length(cycle))),1:length(stoich(1,:))) = out_stoich(length(stoich(:,1))*length(cycle)+((length(stoich(1,:))*(length(cycle)-1)+1):(length(stoich(1,:))*length(cycle))),1:length(stoich(1,:))) + eye(size(out_stoich(length(stoich(:,1))*length(cycle)+((length(stoich(1,:))*(length(cycle)-1)+1):(length(stoich(1,:))*length(cycle))),1:length(stoich(1,:)))));

% cell cycle transitions
for i=1:length(cycle)-1
    out_stoich(end-(length(cycle)-i),length(stoich(1,:))*length(cycle)+i) = -1;
    out_stoich(end-(length(cycle)-i),length(stoich(1,:))*length(cycle)+i+1) = 1;
end
out_stoich(end,length(stoich(1,:))*length(cycle)+i+1) = -1;
out_stoich(end,length(stoich(1,:))*length(cycle)+1) = 2;

% repeat all of the above, but for the reactant designations

out_reac = zeros(size(out_stoich));
% reactions happen in each cycle stage
for i=1:length(cycle)
    out_reac((1:length(stoich(:,1)))+(i-1)*length(stoich(:,1)),(1:length(stoich(1,:)))+(i-1)*length(stoich(1,:))) = rxn_reacts;
end

% cell cycle transitions
for i=1:length(cycle)
    out_reac(end-(length(cycle)-i),length(stoich(1,:))*length(cycle)+i) = 1;
end
% protein transitions due to cell cycle
temp = eye(length(stoich(1,:))*length(cycle));
out_reac(length(stoich(:,1))*length(cycle)+(1:length(stoich(1,:))*length(cycle)),1:length(temp(1,:))) = temp;

% some reactions are synthesis reactions and have no reactants.  Find them
% and make the reactant the corresponding cell type.
synth = find(sum(rxn_reacts,2) == 0);
for i=1:length(cycle)
    for j=1:length(synth)
        out_reac(length(stoich(:,1))*(i-1)+synth(j),end-(length(cycle)-i)) = 1;
    end
end

% make the rates matrix match the appropriate dimensions of out_stoich
temp = repmat(cycle,length(stoich(1,:)),1);
out_rates = [rates(:); temp(:); cycle'];

% this is the shakiest part, and sets the initial abundances
out_C = zeros(length(out_reac(1,:)),1);
out_C(1:length(stoich(1,:))) = C(1:length(stoich(1,:)));
out_C(length(stoich(1,:))*length(cycle)+1) = C(length(stoich(1,:))+1);
% if this condition trips, it's because the input is not a single cell or
% that cell is not in the first cell cycle stage.  Check header for
% description of this assumption.
if length(C(:,1))>1 || C(length(stoich(1,:))+1)==0
    disp('You did not expect this initial cell input, check the code')
end

sol = ode45(@(t,x)odes_from_reactants(t,x,out_reac,out_stoich,out_rates),[0 end_t],out_C);

end