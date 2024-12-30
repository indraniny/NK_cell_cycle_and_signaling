function dx = odes_from_reactants(t,react_C,react_bool,stoich,rates)

% Author: Darren Wethington
% this code is passed to an ODE solver to update protein concentrations

%initialize
dx = zeros(length(react_C),1);
% for each protein
for i=1:length(dx)
    % for each reaction
    for j=1:length(rates)
        % We will want to calculate the product of the rate and all
        % reactant concentrations.  Do that recursively here.
        % initialize with reaction rate
        rxn_contribution = rates(j);
        % for each protein (again)
        for k=1:length(react_bool(j,:))
            % if the protein is a reactant
            if react_bool(j,k) ~= 0
                % multiply the rate by the concentration of that protein
                rxn_contribution = rxn_contribution*react_C(k)^react_bool(j,k);
                % note: I take ^react_bool in case down the line somebody
                % wants to do a reaction like A + A -> B, where the
                % propensity would be rate*A^2.  In that case,
                % rxn_reactants can have values that aren't 0 or 1!
            end
        end
        % Update protein expression change according to stoichiometry and
        % propensity of that reaction
        dx(i) = dx(i) + stoich(j,i)*rxn_contribution;
    end
end

end