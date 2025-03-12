function test_copy_v_symbolic()
% TEST_COPY_V_SYMBOLIC
% -----------------------
% This script tests the 'copy_v' command of the cddmex MEX interface using
% a V-representation provided in symbolic (rational) form.
%
% We define a V-structure with a 3×3 matrix in homogeneous coordinates,
% where the entries are given as rational numbers:
%
%       V = sym([1/1, 0/1, 0/1;
%                1/1, 1/1, 0/1;
%                1/1, 0/1, 1/1]);
%
% This corresponds to the points (0,0), (1,0) and (0,1).
% The copy_v command should simply return a structure whose field V 
% exactly matches the input.
%
% The test compares the output field 'V' with the input field 'V'.

    fprintf('============================================\n');
    fprintf(' Test: cddmex(''copy_v'') with symbolic rational input\n');
    fprintf('============================================\n\n');
    
    % Define the V-representation using rational numbers.
    % Each row is in homogeneous coordinates [1, x, y]
    V = sym([1/1, 0/1, 0/1;
             1/1, 1/1, 0/1;
             1/1, 0/1, 1/1]);
         
    % Create the input structure for the V-representation
    vStruct.V = V;
    
    % Display the input structure
    fprintf('Input V-structure (symbolic rational):\n');
    disp(vStruct);
    
    % Call the MEX function with action 'copy_v'
    fprintf('\nCalling cddmex(''copy_v'', vStruct)...\n');
    try
        out = cddmex('copy_v', vStruct);
    catch ME
        warning('Error calling cddmex(''copy_v''): %s', ME.message);
        return;
    end
    
    % Display the output structure
    fprintf('\nOutput from cddmex(''copy_v''):\n');
    disp(out);
    
    % Verify that the output field "V" exists and matches the input.
    if isfield(out, 'V')
        if isequal(out.V, vStruct.V)
            fprintf('\nTEST PASSED: The copied V matches the input V.\n');
        else
            fprintf('\nTEST FAILED: The copied V does not match the input V.\n');
            fprintf('Input V:\n');
            disp(vStruct.V);
            fprintf('Copied V:\n');
            disp(out.V);
        end
    else
        fprintf('\nTEST FAILED: The output structure does not contain a field "V".\n');
    end
end
