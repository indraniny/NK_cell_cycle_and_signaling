function prop = propensities_from_reactants(react_bool,react_C,rates,states)

% This is a code that calculates propensities for n reactions governing m species in p cells.
% Each of the n reactions can have a rate constant that varies depending on
% the cell state given in states.  Each cell is in one of q states. Thus
% the inputs should be as follows:
% react_bool is a nxm matrix
% react_C is a pxm matrix
% rates is a nxq matrix
% states is a px1 vector

prop = ones(length(rates(:,1)),length(states));
for i=1:length(prop(:,1)) %reaction
    for j=1:length(prop(1,:)) %cell
        prop(i,j) = rates(i,states(j));
        for k=1:length(react_bool(i,:)) %species
            if react_bool(i,k) ~= 0
                prop(i,j) = prop(i,j)*react_C(j,k)^react_bool(i,k);
            end
        end
    end
end
prop = prop';

end