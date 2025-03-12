function test_reduce_h_corrected()
    fprintf('====================================\n');
    fprintf(' Test: cddmex(''reduce_h'') with symbolic H (corrected)\n');
    fprintf('====================================\n\n');

    % Define the H-representation in symbolic form.
    % We represent the region:
    %   x >= 0, y >= 0, x+y <= 1, x+y <= 2 (redundant).
    %
    % In cddmex the H-structure is given in the form: 
    %       B - A*x >= 0.
    % For:
    %   x >= 0   we use:   0 - [-1  0]*x >= 0   (i.e. A = [-1  0], B = 0)
    %   y >= 0   we use:   0 - [ 0 -1]*x >= 0   (i.e. A = [ 0 -1], B = 0)
    %   x+y <= 1 we use:   1 - [ 1  1]*x >= 0   (i.e. A = [ 1  1], B = 1)
    %   x+y <= 2 we use:   2 - [ 1  1]*x >= 0   (i.e. A = [ 1  1], B = 2)
    %
    % Note: Later, FT_set_H_MatrixPtr negates A to recover the original constraint coefficients.
    
    A_sym = sym([
         -1,  0;  % corresponds to x >= 0
          0, -1;  % corresponds to y >= 0
          1,  1;  % corresponds to x+y <= 1
          1,  1;  % corresponds to x+y <= 2 (redundant)
    ]);
    B_sym = sym([0; 0; 1; 2]);
    
    hStruct.A = A_sym;
    hStruct.B = B_sym;
    
    fprintf('Initial H-struct (symbolic):\n');
    disp(hStruct);

    % Call cddmex('reduce_h', hStruct)
    fprintf('\nCalling cddmex(''reduce_h'', hStruct)...\n');
    try
        out = cddmex('reduce_h', hStruct);
    catch ME
        error('Error calling cddmex(''reduce_h''): %s', ME.message);
    end

    fprintf('\nResult from cddmex(''reduce_h''):\n');
    disp(out);

    % Check that the output contains fields A and B.
    if ~isfield(out, 'A') || ~isfield(out, 'B')
        error('Output structure is missing fields "A" or "B".');
    end

    A_red = out.A;
    B_red = out.B;
    
    [mRows, nCols] = size(A_red);
    fprintf('\nReduced H has size %dx%d in A.\n', mRows, nCols);
    disp(A_red);
    disp(B_red);

    % In the minimal (irredundant) representation we expect 3 constraints:
    %    x >= 0, y >= 0, and x+y <= 1.
    if mRows ~= 3
        fprintf('TEST FAILED: Expected 3 rows after reduce_h, but got %d.\n', mRows);
        return;
    end

    % Check that the nonnegativity constraints are preserved:
    % They should appear as rows with B=0 and A = [-1,0] and [0,-1] respectively.
    % Also, the x+y constraint should have B = 1 and A = [1,1].
    % Because FT_set_H_MatrixPtr applies an extra negation to the stored A, the output A is:
    %    [-1,  0]
    %    [ 0, -1]
    %    [ 1,  1]
    % and B is [0; 0; 1].
    expectedA = sym([ -1,  0;
                      0,  -1;
                      1,  1]);
    expectedB = sym([0; 0; 1]);
    
    if ~isequal(A_red, expectedA) || ~isequal(B_red, expectedB)
        fprintf('TEST FAILED: The reduced constraints do not match expected values.\n');
        fprintf('Expected A:\n');
        disp(expectedA);
        fprintf('Expected B:\n');
        disp(expectedB);
    else
        fprintf('TEST PASSED: The redundant row was correctly removed.\n');
    end
end
