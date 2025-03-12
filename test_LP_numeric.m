function test_my_mex()
    % test_my_mex
    % ------------
    % This script defines a small LP problem in symbolic form (using rational
    % numbers), converts the problem data into numeric split-fraction form,
    % calls the MEX function with these numeric values, and checks whether
    % the returned solution matches the known optimum.
    %
    % The LP problem here is:
    %    Maximize    0 + 3 x1 + 4 x2 
    %    Subject to       - 2 x1 -   x2  <= -4/3
    %                            -   x2  <= -2/3
    %                         x1         >= 0
    %                                x2  >= 0
    %
    %       b - Ax >= 0 is standard for cddlib
    %
    %                 4/3 - 2 x1 -   x2  >= 0
    %                 2/3        -   x2  >= 0
    %                         x1         >= 0
    %                                x2  >= 0  
    %
    % The expected solution is:
    %      x1 = 1/3, x2 = 2/3, objective value = 11/3.

    Anum = int64([-2, -1; 
                0, -1; 
                1,  0; 
                0,  1]);
    % All entries have denominator 1.
    Aden = int64(ones(4,2));

    % For vector B (4×1): the values are 4/3, 2/3, 0, 0.
    Bnum = int64([4; 2; 0; 0]);
    Bden = int64([3; 3; 1; 1]);

    % For the objective vector (3×1): the coefficients are 0, 3, and 4.
    objNum = int64([0; 3; 4]);
    % All objective coefficients have denominator 1.
    objDen = int64(ones(3,1));

    % Bundle into a struct with numeric fields for cddmex.
    % Note: The objective fields must be named 'objNum' and 'objDen' (capital N and D).
    lpStruct.Anum   = Anum;
    lpStruct.Aden   = Aden;
    lpStruct.Bnum   = Bnum;
    lpStruct.Bden   = Bden;
    lpStruct.objNum = objNum;
    lpStruct.objDen = objDen;

    % Display some basic info
    fprintf('----- Testing MEX function with a simple LP (Numeric Mode) -----\n');
    fprintf('Problem:\n');
    fprintf('  Maximize 0 + 3 x1 + 4 x2 \n');
    fprintf('  Subject to:\n');
    fprintf('    4/3 - 2 x1 -   x2  >= 0\n');
    fprintf('    2/3        -   x2  >= 0\n');
    fprintf('            x1         >= 0\n');
    fprintf('                   x2  >= 0\n');

    % Call your MEX function with action 'solve_lp_ds'
    try
        solution = cddmex('solve_lp_ds', lpStruct);
    catch ME
        warning('Error encountered while calling cddmex: %s', ME.message);
        return;
    end

    % Extract numeric results directly:
    % The C code returns the solution as structures with fields "num" and "den".
    % For the objective value:
    optValue = double(solution.optvalue.num) / double(solution.optvalue.den);
    % For the primal solution:
    sol = double(solution.sol.num) ./ double(solution.sol.den);
    % For the dual solution:
    dsol = double(solution.dsol.num) ./ double(solution.dsol.den);
    
    % Display the results
    fprintf('\nExtracting primal solution (sol):\n');
    fprintf('%d/%d\n', solution.sol.num, solution.sol.den);
    % disp(sol);
    
    fprintf('\nExtracting dual solution (dsol):\n');
    fprintf('%d/%d\n', solution.dsol.num, solution.dsol.den);
    % disp(dsol);
    
    fprintf('\nExtracting objective value:\n');
    fprintf('%d/%d\n', solution.optvalue.num, solution.optvalue.den);
    % disp(optValue);
    
    fprintf('\nExtracting LP status value:\n');
    lps = solution.how;
    % disp(lps);

    % For checking, extract the relevant values:
    % Note: based on previous test script indexing, x1 is sol(2) and x2 is sol(3).
    x1Val = sol(2);
    x2Val = sol(3);
    
    % Display computed solution
    fprintf('\nComputed solution:\n');
    fprintf('  x1 = %d/%d, x2 = %d/%d\n', solution.sol.num(2), solution.sol.den(2), solution.sol.num(3), solution.sol.den(3));
    fprintf('  objective = %d/%d\n\n', solution.optvalue.num, solution.optvalue.den);

    % Define the known/expected solution as doubles
    x1Expected   = 1/3; 
    x2Expected   = 2/3; 
    objExpected  = 11/3;

    % Compare results within a tolerance
    tol = 1e-12;
    isCorrect = (abs(x1Val - x1Expected) < tol) && ...
                (abs(x2Val - x2Expected) < tol) && ...
                (abs(optValue - objExpected) < tol);

    if isCorrect
        fprintf('Test PASSED. The solution matches the expected optimum.\n');
    else
        fprintf('Test FAILED!\n');
        fprintf('  Expected x1=%g, x2=%g, obj=%g.\n', x1Expected, x2Expected, objExpected);
        fprintf('  Got x1=%g, x2=%g, obj=%g.\n', x1Val, x2Val, optValue);
    end
end
