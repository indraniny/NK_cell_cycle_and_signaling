
% Read the data from the file
filename = 'sensitivity_simulation_results.txt';  % Replace with your actual filename
data = dlmread(filename);

% Extract columns 7 to 10 for the difference calculation
avg_protein = data(:, 7:10);

% Calculate the protein abundance differences betweeb (S-G1), (G2-S), (M-S)
a = avg_protein(:, 2) - avg_protein(:, 1);  % column 8 - column 7
b = avg_protein(:, 3) - avg_protein(:, 2);  % column 9 - column 8
c = avg_protein(:, 4) - avg_protein(:, 3);  % column 10 - column 9

% Extract the last 3 columns for the binary values
p_val = data(:, end-2:end);

% Apply the condition to the p-values: if < 0.05, assign 1; otherwise, 0
p_val_binary = p_val < 0.05;  % Logical matrix where values < 0.05 are 1, else 0

% Convert logical matrix to double (1 and 0)
p_val_binary = double(p_val_binary);

% Initialize the new binary vector q
q = zeros(size(p_val_binary));

% Loop through each row and apply the transformation
for i = 1:size(data, 1)
    % Get the binary vector for the current row
    binVec = p_val_binary(i, :);
    
    % Multiply (a, b, c) by the binary values
    x = a(i) * binVec(1);
    y = b(i) * binVec(2);
    z = c(i) * binVec(3);
    
    % Convert the results into a new binary vector q
    q(i, 1) = sign(x);  % +1 if x > 0, -1 if x < 0, 0 if x == 0
    q(i, 2) = sign(y);  % +1 if y > 0, -1 if y < 0, 0 if y == 0
    q(i, 3) = sign(z);  % +1 if z > 0, -1 if z < 0, 0 if z == 0
end

combined=[data(:,1:6),q];
% Save the result in a new file
dlmwrite('param_protein_var.txt', combined, 'delimiter', '\t', 'precision', 6);

