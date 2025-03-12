function test_my_mex()
    % test_my_mex
    % ------------
    % This script defines a small LP problem in symbolic form (using rational
    % numbers), calls your MEX function, and checks whether the returned solution 
    % matches the known optimum.
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
    %
    % The known optimal solution is x=4, y=0, with objective value 12.


    % Define a symbolic matrix with rational values
    A = sym([
        -2/1, -1/1;    
        0/1, -1/1;    
        1/1, 0/1;
        0/1, 1/1    
    ]);
    B = sym([4/3; 2/3; 0/1; 0/1]);  

    % Define the objective function: 0 + 3 x1 + 4 x2 
    obj = sym([0/1; 3/1; 4/1]);   % Stored as a column vector for cddmex

    % Bundle into a struct (the wrapper cddmex.m expects 'A', 'B', 'obj' fields)
    lpStruct.A = A;
    lpStruct.B = B;
    lpStruct.obj = obj;

    % Display some basic info
    fprintf('----- Testing MEX function with a simple LP -----\n');
    fprintf('Problem:\n');
    fprintf('  Maximize 0 + 3 x1 + 4 x2 \n');
    fprintf('  Subject to:\n');
    fprintf('    4/3 - 2 x1 -   x2  >= 0\n');
    fprintf('    2/3        -   x2  >= 0\n');
    fprintf('            x1         >= 0\n');
    fprintf('                   x2  >= 0\n');

    % Call your MEX function with action 'solve_lp'
try
    solution = cddmex('solve_lp_ds', lpStruct);
catch ME
    warning('Error encountered while calling cddmex: %s', ME.message);
    return;
end


% Extract and display primal solution (sol)
fprintf('\nExtracting primal solution (sol):\n');
sol = solution.sol; % sol is symbolic
disp(sol);

% Extract and display dual solution (dsol)
fprintf('\nExtracting dual solution (dsol):\n');
dsol = solution.dsol; % dsol is symbolic
disp(dsol);

% Extract and display objective value
fprintf('\nExtracting objective value:\n');
optvalue = solution.optvalue; % optvalue is symbolic
disp(optvalue);

% Extract and display LP status value
fprintf('\nExtracting LP status value:\n');
lps = solution.how;
disp(lps);



x1Val = sol(2);
x2Val = sol(3);
objVal = optvalue;

% Display computed solution
fprintf('\nComputed solution:\n');
fprintf('  x1 = %s, x2 = %s\n', x1Val, x2Val);
fprintf('  objective = %s\n\n', objVal);

% Define the known/expected solution
x1Expected   = sym(1/3); 
x2Expected   = sym(2/3); 
objExpected = sym(11/3);

% Compare results with a small tolerance
isCorrect = ( x1Val == x1Expected ) && ...
            ( x2Val == x2Expected ) && ...
            ( objVal == objExpected);

if isCorrect
    fprintf('Test PASSED. The solution matches the expected optimum.\n');
else
    fprintf('Test FAILED!\n');
    fprintf('  Expected x1=%s, x2=%s, obj=%s.\n', x1Expected, x2Expected, objExpected);
    fprintf('  Got x1=%s, x2=%s, obj=%s.\n', x1Val, x2Val, objVal);
end
