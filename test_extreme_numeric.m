function test_extreme()
% test_extreme
% ------------
% This test script defines the H-representation (in numeric rational split form)
% of a triangle described by:
%
%   1) -x  <= -1/2   (i.e., x >= 1/2)
%   2) -y  <= -1/3   (i.e., y >= 1/3)
%   3)  x + y <= 5/4
%
% In split-fraction form:
%   Anum = [ -1,   0;
%              0,  -1;
%              1,   1 ]
%   Aden = ones(3,2)
%   Bnum = [ -1; -1;  5 ]
%   Bden = [  2;  3;  4 ]
%
% The expected extreme vertices are:
%      (1/2, 1/3), (1/2, 3/4), (11/12, 1/3)
%
% The script calls cddmex('extreme', Hstruct), prints each vertex as a fraction
% (num/den), prints the sorted vertices as fractions, and then checks (via a
% double conversion) that the computed vertices match the expected values.
    
    fprintf('\n===== Testing cddmex(''extreme'') on a triangle (numeric rational input) =====\n');
    
    % --- Define the H-representation in numeric rational (split) form ---
    % The inequalities:  B - A*x >= 0
    % Row 1:  -x <= -1/2   <=>  x >= 1/2
    % Row 2:  -y <= -1/3   <=>  y >= 1/3
    % Row 3:   x + y <= 5/4
    Anum = int64([ -1,  0;
                    0, -1;
                    1,  1 ]);
    Aden = int64(ones(3,2));  % All denominators for A are 1.
    Bnum = int64([ -1; -1;  5 ]);
    Bden = int64([  2;  3;  4 ]);
    
    Hstruct = struct('Anum', Anum, 'Aden', Aden, 'Bnum', Bnum, 'Bden', Bden);
    
    fprintf('Input H-structure (numeric rational):\n');
    disp(Hstruct);
    
    % --- Call the MEX function to compute extreme points ---
    try
        result = cddmex('extreme', Hstruct);
    catch ME
        error('Error calling cddmex(''extreme''): %s', ME.message);
    end
    
    fprintf('\nReturned V-structure:\n');
    disp(result);
    
    if ~isfield(result, 'Vnum') || ~isfield(result, 'Vden')
        error('Result does not contain Vnum/Vden fields.');
    end
    
    Vnum = result.Vnum;  % (#vertices x dimension)
    Vden = result.Vden;
    nV   = size(Vnum,1);
    nDim = size(Vnum,2);
    
    % --- Print raw vertices as fractions ---
    fprintf('\n--- Vertices as fractions (raw) ---\n');
    for i = 1:nV
        fprintf('  vertex %d: ', i);
        for j = 1:nDim
            fprintf('%d/%d  ', Vnum(i,j), Vden(i,j));
        end
        fprintf('\n');
    end
    
    % --- Sort vertices for testing ---
    % Convert vertices to double (only for sorting/comparison).
    vertices_double = double(Vnum) ./ double(Vden);
    [vertices_sorted, sortIdx] = sortrows(vertices_double);
    
    % Also compute sorted fraction data (without floating conversion)
    Vnum_sorted = Vnum(sortIdx,:);
    Vden_sorted = Vden(sortIdx,:);
    
    fprintf('\n--- Sorted vertices (as fractions) ---\n');
    for i = 1:nV
        fprintf('  sorted vertex %d: ', i);
        for j = 1:nDim
            fprintf('%d/%d  ', Vnum_sorted(i,j), Vden_sorted(i,j));
        end
        fprintf('\n');
    end
    
    fprintf('\n--- Sorted vertices (as double for testing) ---\n');
    disp(vertices_sorted);
    
    % --- Define expected vertices for the triangle ---
    % Expected vertices (as double):
    %   (1/2, 1/3)       -> [0.5, 0.333333333333333]
    %   (1/2, 3/4)       -> [0.5, 0.75]
    %   (11/12, 1/3)     -> [0.916666666666667, 0.333333333333333]
    expected = [ 0.5,                0.333333333333333;
                 0.5,                0.75;
                 0.916666666666667,   0.333333333333333 ];
    expected_sorted = sortrows(expected);
    
    tol = 1e-12;
    if all(abs(vertices_sorted - expected_sorted) < tol, 'all')
        fprintf('TEST PASSED: Computed vertices match the expected triangle extreme points.\n');
    else
        fprintf('TEST FAILED: Computed vertices do not match the expected values.\n');
        fprintf('Computed sorted vertices:\n');
        disp(vertices_sorted);
        fprintf('Expected sorted vertices:\n');
        disp(expected_sorted);
    end
end
