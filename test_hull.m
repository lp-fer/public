function test_hull_rational()
% TEST_HULL_RATIONAL
% This test verifies that cddmex('hull', ...) correctly converts a V-representation
% (provided in rational (numeric) mode) to an H-representation for the unit square in 2D.
%
% For the unit square with vertices:
%   (0,0), (1,0), (0,1), (1,1)
% the homogeneous V-representation (each row: [1, x, y]) is:
%
%   [1, 0, 0]
%   [1, 1, 0]
%   [1, 0, 1]
%   [1, 1, 1]
%
% In numeric (rational) mode the wrapper expects split fields for the non-homogenized
% part of the vertices (the homogenizing coordinate is fixed to 1).
% Thus, we define:
%
%   Vnum = [0 0; 1 0; 0 1; 1 1]  (4x2 matrix)
%   Vden = ones(4,2)              (4x2 matrix)
%
% One acceptable H–representation (inequalities in the form b - A*x >= 0) is:
%
%   [0, -1,  0]    %  x >= 0
%   [0,  0, -1]    %  y >= 0
%   [1,  1,  0]    %  x <= 1
%   [1,  0,  1]    %  y <= 1
%
% Since H–representations are not unique, we compare the sorted rows.

fprintf('\n===== Testing cddmex(''hull'') on the unit square in 2D (rational input) =====\n');

% Define the non-homogeneous part of the vertices as rational numbers.
% (Here, we assume each entry is given as an integer with denominator 1.)
Vnum = int64([0 0; 
              1 0; 
              0 1; 
              1 1]);  % 4x2: each row is [x, y]
Vden = int64(ones(4,2));  % Denominators (all ones)

% Build the input V-structure (numeric mode) for the hull command.
vStruct = struct('Vnum', Vnum, 'Vden', Vden);

% Call the hull command via the wrapper.
try
    result = cddmex('hull', vStruct);
catch ME
    warning('Error calling cddmex(''hull''): %s', ME.message);
    return;
end

disp('Result from cddmex(''hull'') (numeric output):');
disp(result);

% Check that the output has the required split fields.
if ~isfield(result, 'Anum') || ~isfield(result, 'Aden') || ...
   ~isfield(result, 'Bnum') || ~isfield(result, 'Bden')
    error('Test failed: The output H-structure is missing required fields.');
end

% Reconstruct the H matrix:
% The H-structure represents each inequality as [B, A] where B is the right-hand side
% and A contains the coefficients. (The field Bnum/Bden is a column vector and
% Anum/Aden is an m x (n-1) matrix.)
B = double(result.Bnum) ./ double(result.Bden);  % m x 1
A = double(result.Anum) ./ double(result.Aden);   % m x (n-1)
H_computed = [B, A];  % m x n

% Define one acceptable expected H-representation.
% Here, we use the following set of inequalities:
%   x >= 0  ->  0 - (-1)*x >= 0  -> row: [0, -1,  0]
%   y >= 0  ->  0 - (-1)*y >= 0  -> row: [0,  0, -1]
%   x <= 1  ->  1 - (1*x)   >= 0  -> row: [1,  1,  0]
%   y <= 1  ->  1 - (1*y)   >= 0  -> row: [1,  0,  1]
H_expected = [0, -1,  0;
              0,  0, -1;
              1,  1,  0;
              1,  0,  1];

% Sort the rows for comparison (since ordering may vary).
H_computed_sorted = sortrows(H_computed);
H_expected_sorted = sortrows(H_expected);

tol = 1e-12;
if all(abs(H_computed_sorted - H_expected_sorted) < tol, 'all')
    fprintf('PASSED: The computed H-representation matches the expected inequalities for the unit square.\n');
else
    fprintf('FAILED: The computed H-representation does not match the expected values.\n');
    disp('Computed H (sorted):');
    disp(H_computed_sorted);
    disp('Expected H (sorted):');
    disp(H_expected_sorted);
end

end
