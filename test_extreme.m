function test_extreme_unit_square()
    % TEST_EXTREME_UNIT_SQUARE  Test script for cddmex('extreme') on a simple square.
    %
    % The square is defined in 2D by
    %   0 <= x <= 1
    %   0 <= y <= 1
    %
    % In cdd's convention (b - A x >= 0 => A x <= b), the H-representation is:
    %   1 - x >= 0   (x <= 1)
    %   1 - y >= 0   (y <= 1)
    %   0 - (-x) >= 0  =>  x >= 0
    %   0 - (-y) >= 0  =>  y >= 0
    %
    % We define these constraints in symbolic form and pass them to cddmex('extreme').
    % The code should return 4 vertices (0,0), (1,0), (0,1), (1,1) with no rays.

    % Define the constraints in symbolic form:
    %   row1:  b=1, A=[1, 0] => 1 - x >= 0 => x <= 1
    %   row2:  b=1, A=[0, 1] => 1 - y >= 0 => y <= 1
    %   row3:  b=0, A=[-1, 0] => -x <= 0 => x >= 0
    %   row4:  b=0, A=[0, -1] => -y <= 0 => y >= 0
    %
    % We'll store them in "A" and "B" as per the usual "symbolic H-mode" in the cddmex wrapper.
    % Then we call cddmex('extreme', Hstruct).

    % Symbolically define the 4x2 matrix A:
    A = sym([1  0;  % row1
             0  1;  % row2
            -1  0;  % row3
             0 -1]);% row4

    % Symbolic B (4x1):
    B = sym([1; 1; 0; 0]);

    % Create an input struct. For an H-representation,
    % cddmex typically looks for fields named 'A' and 'B' (symbolic),
    % or 'Anum','Aden','Bnum','Bden' (numeric). We'll do symbolic here.
    Hstruct = struct('A', A, 'B', B);

    % We want the extreme points, so we call:  cddmex('extreme', Hstruct)
    fprintf('\n===== Testing cddmex(''extreme'') on the unit square in 2D =====\n');

    try
        result = cddmex('extreme', Hstruct);
    catch ME
        warning('Error calling cddmex(''extreme''): %s', ME.message);
        return;
    end

    % The MEX call should return a V-structure with fields:
    %  - Vnum, Vden, Rnum, Rden, vpos, rpos, etc. (if numeric)
    % or, if your wrapper merges them into symbolic, something like
    %  - result.V, result.R, ...
    % Exactly how it’s returned depends on your wrapper's design.
    %
    % In the typical "symbolic" path, your wrapper might reassemble
    % the rational data into fields like result.V, result.R (both symbolic).
    % Let’s check how many vertices and rays we got:

    disp('Result from cddmex(''extreme''):');
    disp(result);

    % Usually, for a 2D bounded polytope, we expect R to be empty and
    % V to be the set of corners.  Let's assume your wrapper merges them
    % into a field result.V which is Nx2 (symbolic).
    %
    % If your wrapper merges them differently, adjust accordingly.

    if ~isfield(result, 'V') 
        warning('No field "V" found in result. Check your wrapper output.');
        return;
    end

    Vmat = result.V;   % Suppose result.V is a Nx2 symbolic matrix
    disp('Computed vertices (symbolic):');
    disp(Vmat);

    % The unit square's corners are (0,0), (1,0), (0,1), (1,1).
    % Let’s see if the set of rows in Vmat matches that set.

    % Because cdd might return them in any order, we do a set comparison:
    expected = sym([0 0; 
                    1 0; 
                    0 1; 
                    1 1]);

    % Sort them (just to compare in a canonical order).
    VmatSorted      = sortrows(double(Vmat));      % convert symbolic -> double for sorting
    expectedSorted  = sortrows(double(expected));

    % Check if they're the same within a tolerance:
    tol = 1e-12;
    diffVal = abs(VmatSorted - expectedSorted);
    maxDiff = max(diffVal(:));
    if maxDiff < tol
        fprintf('PASSED: The extreme points match the expected unit square corners.\n');
    else
        fprintf('FAILED: The computed vertices differ from the expected corners.\n');
        disp('Difference:');
        disp(diffVal);
    end
end
