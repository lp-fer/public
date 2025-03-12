function test_find_interior_numeric()
% TEST_FIND_INTERIOR_NUMERIC
% -----------------------------------------
% This test calls cddmex('find_interior') on the region:
%
%     0 <= x <= 1,   0 <= y <= 1,
%
% by encoding each inequality in the cdd "b - A*x >= 0" format,
% storing all coefficients as integers in Anum/Aden, Bnum/Bden.

    fprintf('================================================\n');
    fprintf(' Test: cddmex(''find_interior'') with numeric rationals\n');
    fprintf('================================================\n\n');

    %----------------------------------------------------------
    % 1) Define the constraints in split form
    %----------------------------------------------------------
    %
    % We want:
    %    x >= 0   =>   0 - [(-1)*x + 0*y] >= 0
    %    x <= 1   =>   1 - [(+1)*x + 0*y] >= 0
    %    y >= 0   =>   0 - [0*x + (-1)*y] >= 0
    %    y <= 1   =>   1 - [0*x + (+1)*y] >= 0
    %
    % So we have 4 rows, each with 2 columns. 
    %   Row1 => b=0, A=[-1,  0] 
    %   Row2 => b=1, A=[+1,  0] 
    %   Row3 => b=0, A=[ 0, -1]
    %   Row4 => b=1, A=[ 0, +1]
    %
    % We'll store A and B as int64 arrays for the "numeric rational" mode.

    Anum = int64([
       -1,  0;   % x >= 0
        1,  0;   % x <= 1
        0, -1;   % y >= 0
        0,  1    % y <= 1
    ]);
    Aden = int64([
        1, 1;
        1, 1;
        1, 1;
        1, 1
    ]);  % All denominators = 1

    Bnum = int64([0; 1; 0; 1]);
    Bden = int64([1; 1; 1; 1]);

    % Assemble them into an H-struct
    hStruct.Anum = Anum;
    hStruct.Aden = Aden;
    hStruct.Bnum = Bnum;
    hStruct.Bden = Bden;

    % Show the input
    fprintf('H-struct in numeric rational form:\n');
    disp(hStruct);

    %----------------------------------------------------------
    % 2) Call cddmex('find_interior', hStruct)
    %----------------------------------------------------------
    fprintf('\nCalling cddmex(''find_interior'', hStruct)...\n');
    try
        out = cddmex('find_interior', hStruct);
    catch ME
        warning('Error calling cddmex(''find_interior''): %s', ME.message);
        return;
    end

    fprintf('\nResult from cddmex(''find_interior''):\n');
    disp(out);

    %----------------------------------------------------------
    % 3) Extract numeric solution fields
    %----------------------------------------------------------
    if ~isfield(out, 'solnum') || ~isfield(out, 'solden')
        error('No "solnum"/"solden" fields found. Possibly the region is empty or something went wrong.');
    end

    solnum = double(out.solnum);  % convert int64 -> double
    solden = double(out.solden);

    % The dimension d = #vars + 1. For 2D, we expect length(sol)=3.
    if length(solnum) < 3
        fprintf('\nTEST FAILED: The solution vector is shorter than expected.\n');
        return;
    end

    % Homogeneous coordinate is sol(1), then x=sol(2), y=sol(3).
    xVal = solnum(2)/solden(2);
    yVal = solnum(3)/solden(3);

    fprintf('\nInterior point found:\n');
    fprintf('  x = %.4f,  y = %.4f\n', xVal, yVal);

    %----------------------------------------------------------
    % 4) Verify that (xVal,yVal) is indeed in the interior
    %----------------------------------------------------------
    epsTol = 1e-14;
    if (xVal > 0+epsTol) && (xVal < 1-epsTol) && ...
       (yVal > 0+epsTol) && (yVal < 1-epsTol)
        fprintf('TEST PASSED: Found a strictly interior point.\n');
    else
        fprintf('TEST FAILED: The point is not strictly in the interior.\n');
    end

end
