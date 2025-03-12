function test_copy_h()
    % test_copy_h
    % ------------
    % This script defines a simple H‐structure (inequalities)
    % in symbolic (rational) form, calls the cddmex('copy_h', ...) command,
    % and then verifies that the returned structure is identical to the input.
    %
    % For this example, the H‐structure represents the two inequalities:
    %    1 - ( x + y ) >= 0   (i.e.  x + y <= 1 )
    %    1 - (-x + y ) >= 0   (i.e. -x + y <= 1 )
    %
    % The H‐structure is provided in symbolic form with fields A and B.
    
    fprintf('----- Testing cddmex(''copy_h'') with symbolic H-structure -----\n');
    
    % Define the inequality coefficients (A) and right-hand side (B)
    % We write the inequalities in the form: B - A*x >= 0.
    % For the inequality: x + y <= 1, we write: 1 - [1 1]*x >= 0.
    % For the inequality: -x + y <= 1, we write: 1 - [-1 1]*x >= 0.
    A = sym([ 1/1,  1/1;
             -1/1,  1/1]);
    B = sym([1/1; 1/1]);
    
    % Create the input H-structure (symbolic mode)
    hStruct.A = A;
    hStruct.B = B;
    
    % Display the input H-structure
    fprintf('\nInput H-structure (symbolic):\n');
    disp(hStruct);
    
    % Call the copy_h command using the wrapper
    try
        hCopy = cddmex('copy_h', hStruct);
    catch ME
        error('Error calling cddmex(''copy_h''): %s', ME.message);
    end
    
    % Display the copied H-structure
    fprintf('\nOutput H-structure from cddmex(''copy_h''):\n');
    disp(hCopy);
    
    % Verify that the copied structure matches the original
    % (We compare the symbolic matrices A and B.)
    isEqualA = isequal(hStruct.A, hCopy.A);
    isEqualB = isequal(hStruct.B, hCopy.B);
    
    if isEqualA && isEqualB
        fprintf('\nTest PASSED: The copied H-structure matches the input.\n');
    else
        fprintf('\nTest FAILED: The copied H-structure does not match the input.\n');
        if ~isEqualA
            fprintf('Field A differs.\n');
        end
        if ~isEqualB
            fprintf('Field B differs.\n');
        end
    end
end
