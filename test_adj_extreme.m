function test_adj_extreme_rational_corrected()
% test_adj_extreme_rational_corrected
% ------------------------------------
% This test verifies that cddmex('adj_extreme', ...) correctly computes
% both the V-representation (extreme generators) and the adjacency information
% from an H-representation provided in rational (numeric) mode.
%
% For the unit square (H-representation given as b - A*x >= 0):
%   Row 1:  1 - [ 1  0] * x >= 0    (x <= 1)
%   Row 2:  1 - [ 0  1] * x >= 0    (y <= 1)
%   Row 3:  0 - [-1  0] * x >= 0    (x >= 0)
%   Row 4:  0 - [ 0 -1] * x >= 0    (y >= 0)
%
% The expected extreme generators are (up to ordering):
%   (0,0), (1,0), (0,1), (1,1)
%
% For the square, the vertex adjacency (in cyclic order) is:
%   (0,0) adjacent to (0,1) and (1,0)
%   (0,1) adjacent to (0,0) and (1,1)
%   (1,0) adjacent to (0,0) and (1,1)
%   (1,1) adjacent to (0,1) and (1,0)
%
% In the test below we:
%  1. Call cddmex('adj_extreme') to get vStruct and adj (unsorted).
%  2. Reconstruct the vertices from vStruct and sort them (using sortrows).
%  3. Compute the inverse permutation to remap the adjacency cell array so that
%     its indices refer to the sorted order.
%  4. Check that the sorted vertices match the expected coordinates, and that
%     the re-labeled adjacency cell array matches the expected values.
%
% Expected sorted order (via sortrows) is:
%   [0 0]; [0 1]; [1 0]; [1 1]
% with expected sorted adjacency:
%   Vertex 1: [2 3]
%   Vertex 2: [1 4]
%   Vertex 3: [1 4]
%   Vertex 4: [2 3]

    fprintf('\n===== Testing cddmex(''adj_extreme'') on the unit square (rational input) =====\n');

    % Construct the H-representation (numeric, rational mode)
    Bnum = int64([1; 1; 0; 0]);         % RHS vector (4x1)
    Bden = int64(ones(4,1));             % Denominators (all 1)
    Anum = int64([ 1,  0;
                   0,  1;
                  -1,  0;
                   0, -1]);           % Coefficients (4x2)
    Aden = int64(ones(4,2));             % Denominators for A
    Hstruct = struct('Bnum', Bnum, 'Bden', Bden, 'Anum', Anum, 'Aden', Aden);

    fprintf('Input H-structure (rational):\n');
    disp(Hstruct);

    % Call the adj_extreme command with two outputs.
    try
        [vStruct, adj] = cddmex('adj_extreme', Hstruct);
    catch ME
        error('Error calling cddmex(''adj_extreme''): %s', ME.message);
    end

    fprintf('\nResulting V-structure (extreme generators):\n');
    disp(vStruct);
    fprintf('\nResulting adjacency cell array (original order):\n');
    disp(adj);

    % Reconstruct vertices from the V-structure.
    % In numeric mode, the fields Vnum and Vden are returned.
    vertices = double(vStruct.Vnum) ./ double(vStruct.Vden);  % each row is a vertex

    % Sort the vertices (using sortrows) and obtain the permutation.
    [vertices_sorted, sortedIdx] = sortrows(vertices);
    
    % Compute the inverse permutation: for each original index, find its new sorted position.
    invIdx = zeros(size(sortedIdx));
    for k = 1:length(sortedIdx)
        invIdx(sortedIdx(k)) = k;
    end

    % Remap the adjacency cell array:
    % For each vertex in sorted order (i=1:n), the original vertex index is sortedIdx(i).
    % Then map each neighbor j in adj{sortedIdx(i)} to its new index using invIdx(j).
    n = numel(adj);
    adj_sorted = cell(n, 1);
    for i = 1:n
        origIndex = sortedIdx(i);
        origNeighbors = double(adj{origIndex});
        % Remap each neighbor index.
        newNeighbors = arrayfun(@(j) invIdx(j), origNeighbors);
        newNeighbors = sort(newNeighbors);  % sort for consistency
        adj_sorted{i} = newNeighbors;
    end

    fprintf('\nRe-labeled adjacency cell array (sorted order):\n');
    disp(adj_sorted);

    % Expected sorted vertices (via sortrows) for the unit square.
    expected_vertices = [0 0; 0 1; 1 0; 1 1];
    tol = 1e-12;
    if all(abs(vertices_sorted - expected_vertices) < tol, 'all')
        fprintf('PASSED: The computed extreme generators match the expected vertices.\n');
    else
        fprintf('FAILED: The computed vertices do not match the expected values.\n');
        disp('Computed vertices (sorted):');
        disp(vertices_sorted);
        disp('Expected vertices (sorted):');
        disp(expected_vertices);
    end

    % Expected adjacency for the sorted order:
    % With sorted vertices as: 1: (0,0), 2: (0,1), 3: (1,0), 4: (1,1)
    expected_adj = { [2 3], [1 4], [1 4], [2 3] };

    adjacency_passed = true;
    for i = 1:n
        computed_neighbors = adj_sorted{i};
        expected_neighbors = expected_adj{i};
        if ~isequal(computed_neighbors, expected_neighbors)
            fprintf('Mismatch for sorted vertex %d: expected [%s], got [%s]\n', ...
                i, num2str(expected_neighbors), num2str(computed_neighbors));
            adjacency_passed = false;
        end
    end

    if adjacency_passed
        fprintf('PASSED: Re-labeled adjacency output matches expected values.\n');
    else
        fprintf('FAILED: Re-labeled adjacency output does not match expected values.\n');
    end
end
