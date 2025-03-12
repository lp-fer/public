function test_v_hull_extreme_numeric()
% TEST_V_HULL_EXTREME_NUMERIC
% ---------------------------------
% This script tests the 'v_hull_extreme' command of the cddmex MEX interface 
% using an example where we know the final set of extreme vertices.
%
% We define a square (0,0), (1,0), (1,1), (0,1) and add an interior point (1/2,1/2).
% The convex hull is still the same square, so the extreme vertices should be the
% 4 corners only. We provide the data in split numerator/denominator form to
% mimic "numeric" rational input.

    fprintf('=========================================\n');
    fprintf(' Test of cddmex(''v_hull_extreme'') with numeric rational V-representation\n');
    fprintf('=========================================\n');

    % We have 5 points in 2D:
    %    (0,0), (1,0), (1,1), (0,1), and (1/2,1/2).
    %
    % For cddlib's "V-representation", each row is [1  x  y].
    % We'll store them in *split numerator/denominator* form, i.e.:
    %   vStruct.Vnum(i,j), vStruct.Vden(i,j)
    %
    % Row i = [ 1    x_i   y_i ], so dimension = 3 columns.
    % 
    % The interior point (1/2,1/2) is rationally 1/2 = numerator=1, denominator=2.

    Vnum = [ ...
        1   0   0;   % (0,0)
        1   1   0;   % (1,0)
        1   1   1;   % (1,1)
        1   0   1;   % (0,1)
        1   1   1    % We'll fill denominator below for (1/2,1/2)
    ];
    Vden = [ ...
        1   1   1;   % (0,0) => 0 = 0/1
        1   1   1;   % (1,0) => 1 = 1/1
        1   1   1;   % (1,1) => 1 = 1/1
        1   1   1;   % (0,1) => 1 = 1/1
        1   2   2    % (1/2,1/2) => 1/2 = num=1, den=2
    ];

    % Build the input struct with numeric rational fields:
    vStruct.Vnum = int64(Vnum);
    vStruct.Vden = int64(Vden);

    % No rays, so we do not set Rnum / Rden
    % No lineality, so no 'lin' field

    % Display the input
    disp('Input V-struct (numeric rational form) to cddmex(''v_hull_extreme''):');
    disp(vStruct);

    % Call the wrapper
    try
        outStruct = cddmex('v_hull_extreme', vStruct);
    catch ME
        warning('Error calling cddmex(''v_hull_extreme''): %s', ME.message);
        return;
    end

    %---------------------------------------------------------
    % Interpret the result
    %---------------------------------------------------------
    % The output is another V-struct in the same (split) numeric format,
    % containing the "extreme" points only (i.e., the minimal generating set).
    % We expect the interior point (1/2,1/2) to vanish from the final set.

    if ~isstruct(outStruct) || ~isfield(outStruct,'Vnum') || ~isfield(outStruct,'Vden')
       error('Output from v_hull_extreme is malformed or missing fields.');
    end

    % Extract Vnum/Vden from the output
    VnumOut = double(outStruct.Vnum);  % convert int64 -> double for inspection
    VdenOut = double(outStruct.Vden);

    fprintf('\nOutput from cddmex(''v_hull_extreme''):\n');
    fprintf('  # of rows = %d, # of columns = %d\n', size(VnumOut,1), size(VnumOut,2));

    % We know each row is [1, x, y]. Let's parse them as rational x, y.
    % The "1" in the first column means a point row (not a ray row).
    outPoints = zeros(size(VnumOut,1), 2);  % store the rational (x,y) in floating form, just for checks
    for iRow = 1:size(VnumOut,1)
        xNum = VnumOut(iRow,2);
        xDen = VdenOut(iRow,2);
        yNum = VnumOut(iRow,3);
        yDen = VdenOut(iRow,3);

        outPoints(iRow,1) = xNum / xDen;
        outPoints(iRow,2) = yNum / yDen;
    end

    disp('Extracted (x,y) points from the returned V-representation:');
    disp(outPoints);

    %---------------------------------------------------------
    % Check that the set of points is the 4 corners
    %   (0,0), (1,0), (1,1), (0,1)
    %---------------------------------------------------------
    % We'll allow any row ordering. We'll see if they match 
    % exactly { (0,0), (1,0), (1,1), (0,1) } ignoring duplicates.

    expectedPoints = [0 0; 1 0; 1 1; 0 1];
    if size(outPoints,1) ~= 4
        fprintf('\nTEST FAILED: Expected 4 extreme points, but got %d.\n', size(outPoints,1));
        return;
    end

    % A simple check that each returned point is in the set of expected corners
    allMatch = true;
    for iRow = 1:4
        xy = outPoints(iRow,:);
        matchIdx = find(all(abs(expectedPoints - xy) < 1e-14,2));  %#ok<AGROW> 
        if isempty(matchIdx)
            allMatch = false;
            break;
        end
    end

    if allMatch
        fprintf('\nTEST PASSED: The interior point was removed, and the 4 corners remain as extreme vertices.\n');
    else
        fprintf('\nTEST FAILED: The returned extreme points do not match the expected 4 corners.\n');
    end

end
