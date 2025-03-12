function test_adjacency_symbolic_corrected()
% test_adjacency_symbolic_corrected
% -------------------------
% This script creates a V-representation of the unit square in symbolic mode,
% with vertices given in affine coordinates, so that the homogenizing coordinate
% is added automatically by the MEX function. It then calls cddmex('adjacency', Vstruct)
% and verifies that the computed adjacency matches the expected values.
%
% For the unit square, the affine vertices are:
%    (0,0), (1,0), (0,1), (1,1)
%
% The expected adjacency (1-indexed) is:
%    Vertex 1: [2, 3]
%    Vertex 2: [1, 4]
%    Vertex 3: [1, 4]
%    Vertex 4: [2, 3]

    fprintf('----- Testing cddmex(''adjacency'') with a unit square V-structure (Symbolic Mode, corrected) -----\n');

    % Define the affine coordinates of the unit square vertices (4x2)
    V_affine = [0/1 0/1;
                1/1 0/1;
                0/1 1/1;
                1/1 1/1];
    V_affine = sym(V_affine);  % Convert to symbolic exact rationals

    % Create the V-structure (field 'V') using affine coordinates.
    vStruct.V = V_affine;

    fprintf('Input V-structure (affine, symbolic):\n');
    disp(vStruct);

    % Call the MEX function using the wrapper.
    try
        adj_out = cddmex('adjacency', vStruct);
    catch ME
        error('Error calling cddmex(''adjacency''): %s', ME.message);
    end

    fprintf('\nOutput adjacency cell array (symbolic):\n');
    disp(adj_out);

    % Expected adjacency (1-indexed):
    %   Vertex 1: adjacent to [2, 3]
    %   Vertex 2: adjacent to [1, 4]
    %   Vertex 3: adjacent to [1, 4]
    %   Vertex 4: adjacent to [2, 3]
    expected = { sym([2 3]), sym([1 4]), sym([1 4]), sym([2 3]) };

    % Compare each vertex's computed adjacency with expected values.
    numVertices = numel(expected);
    isCorrect = true;
    for i = 1:numVertices
        if i > numel(adj_out)
            isCorrect = false;
            fprintf('Output cell array has fewer cells than expected.\n');
            break;
        end
        % Convert each cell to a double array for comparison.
        out_i = double(adj_out{i});
        out_i_sorted = sort(out_i);
        exp_i_sorted = sort(double(expected{i}));
        if ~isequal(out_i_sorted, exp_i_sorted)
            isCorrect = false;
            fprintf('Mismatch for vertex %d: expected [%s], got [%s]\n', ...
                i, num2str(exp_i_sorted), num2str(out_i_sorted));
        end
    end

    if isCorrect
        fprintf('TEST PASSED: Adjacency output matches expected values.\n');
    else
        fprintf('TEST FAILED: Adjacency output does not match expected values.\n');
    end
end
