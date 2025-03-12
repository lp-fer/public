function test_reduce_v_symbolic()
% TEST_REDUCE_V_SYMBOLIC
% ---------------------------------
% This script tests the 'reduce_v' command of the cddmex MEX interface
% using symbolic (rational) input.  We define a 2D set:
%
%   * A triangle formed by the points (0,0), (1,0), (0,1)
%   * An additional interior point (1/2, 1/2)
%
% All points are stored in homogeneous V-representation rows of the form [1, x, y].
% Because the interior point does not affect the convex hull, we expect 'reduce_v'
% to remove it, leaving only the three triangle corners.

    fprintf('======================================\n');
    fprintf(' Test: cddmex(''reduce_v'') with symbolic V-struct\n');
    fprintf('======================================\n\n');

    %---------------------------------------------------------
    % 1) Define a V-representation in symbolic (rational) form
    %---------------------------------------------------------
    % We have four points in 2D:
    %   (0,0), (1,0), (0,1), (1/2,1/2).
    % Each row: [1, x, y] stored as symbolic rationals.

    V = sym([
        1,  0,   0;      % corner #1
        1,  1,   0;      % corner #2
        1,  0,   1;      % corner #3
        1,  1/2, 1/2     % an interior point
    ]);

    % Create the input structure for cddmex (symbolic mode)
    vStruct.V = V;

    % Display the input
    fprintf('Initial V-struct (symbolic):\n');
    disp(vStruct);

    %---------------------------------------------------------
    % 2) Call the MEX function: cddmex('reduce_v', ...)
    %---------------------------------------------------------
    fprintf('\nCalling cddmex(''reduce_v'', vStruct)...\n');
    try
        out = cddmex('reduce_v', vStruct);
    catch ME
        warning('Error calling cddmex(''reduce_v''): %s', ME.message);
        return;
    end

    fprintf('\nResult from cddmex(''reduce_v''):\n');
    disp(out);

    %---------------------------------------------------------
    % 3) Verify the redundant interior point was removed
    %---------------------------------------------------------
    if ~isfield(out, 'V')
        error('Output structure missing field "V". The MEX call did not return a valid V-struct.');
    end

    % The returned 'V' is symbolic if our input was symbolic
    V_reduced = out.V;
    [numRows, numCols] = size(V_reduced);

    fprintf('\nReduced V has size %dx%d.\n', numRows, numCols);
    disp(V_reduced);

    % Expect exactly 3 rows (the triangle corners).
    % Let's check that row count is 3:
    if numRows ~= 3
        fprintf('TEST FAILED: Expected 3 corner points after reduce_v, but got %d.\n', numRows);
        return;
    end

    % Each row is of the form [1, x, y]. Let's parse out x,y.
    corners = double(V_reduced(:,2:3));  % convert to double for easier comparison
    disp('Corners from the reduced V-struct (converted to numeric):');
    disp(corners);

    % Our expected corner set (in any order):
    expected = [0 0; 1 0; 0 1];

    % Check that each row in 'corners' matches something in 'expected'
    isMatch = true;
    for i = 1:size(corners,1)
        row = corners(i,:);
        found = any(all(abs(expected - row) < 1e-14, 2));
        if ~found
            isMatch = false;
            break;
        end
    end

    % Print final PASS/FAIL
    if isMatch && size(corners,1) == 3
        fprintf('TEST PASSED: Redundant interior point was removed, leaving the 3 triangle vertices.\n');
    else
        fprintf('TEST FAILED: The points in the reduced V-struct do not match the expected corners.\n');
    end

end
