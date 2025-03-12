function test_adj_extreme_fraction_display()
% test_adj_extreme_fraction_display
% ----------------------------------
% Demonstrates how to supply rational H-data as Anum, Aden, Bnum, Bden,
% call cddmex('adj_extreme', ...), and then print each returned vertex as
% "num/den" (without converting to floating point). At the end the test
% checks if the computed vertices and adjacency match the expected results.
%
% The region is described by b - A*x >= 0 (i.e., A*x <= b) with 3 constraints:
%   1) -x  <= -1/2   (i.e., x >= 1/2)
%   2) -y  <= -1/3   (i.e., y >= 1/3)
%   3)  x + y <= 5/4
%
% Expected extreme vertices (in any order) are:
%   (1/2, 1/3), (1/2, 3/4), (11/12, 1/3)
% When sorted (by rows) they become:
%      [0.5,                0.333333333333333];
%      [0.5,                0.75];
%      [0.916666666666667,   0.333333333333333]
%
% In a triangle every vertex is adjacent to the other two. When the vertices
% are sorted as above, the expected adjacency is:
%      Vertex 1: [2 3]
%      Vertex 2: [1 3]
%      Vertex 3: [1 2]

    fprintf('\n===== Testing cddmex(''adj_extreme'') with fraction display =====\n');

    % --- Define the H-representation in numeric rational (split) form ---
    % For row i, the inequality is: B(i) - A(i,:)*x >= 0.
    %
    % Row 1:  -x <= -1/2   <=>  x >= 1/2
    % Row 2:  -y <= -1/3   <=>  y >= 1/3
    % Row 3:   x + y <= 5/4
    %
    % Thus:
    %   Anum = [ -1,   0;
    %              0,  -1;
    %              1,   1 ]
    %   Aden = ones(3,2)
    %   Bnum = [ -1; -1;  5 ]
    %   Bden = [  2;  3;  4 ]
    Anum = int64([ -1,  0;
                    0, -1;
                    1,  1 ]);
    Aden = int64(ones(3,2));
    Bnum = int64([ -1; -1;  5 ]);
    Bden = int64([  2;  3;  4 ]);
    
    % Bundle into an H-struct
    Hstruct = struct('Anum', Anum, 'Aden', Aden, ...
                     'Bnum', Bnum, 'Bden', Bden);
    
    fprintf('Input H-structure:\n');
    disp(Hstruct);
    
    % --- Call the MEX function ---
    try
        [vStruct, adjacency] = cddmex('adj_extreme', Hstruct);
    catch ME
        error('Error calling cddmex(''adj_extreme''): %s', ME.message);
    end
    
    fprintf('\nReturned adjacency (un‐permuted):\n');
    disp(adjacency);
    
    fprintf('Returned V-structure:\n');
    disp(vStruct);
    
    % --- Print each vertex as a fraction (num/den) ---
    Vnum = vStruct.Vnum;  % (#vertices x dimension)
    Vden = vStruct.Vden;
    nV   = size(Vnum,1);
    
    fprintf('\n--- Vertices as fractions (num/den) ---\n');
    for i = 1:nV
        fprintf('  vertex %d: ', i);
        for j = 1:size(Vnum,2)
            num = Vnum(i,j);
            den = Vden(i,j);
            fprintf('%d/%d  ', num, den);
        end
        fprintf('\n');
    end
    
    % --- For testing, sort the vertices and remap adjacency ---
    % Convert vertices to double for sorting and comparison purposes only.
    vertices_double = double(Vnum) ./ double(Vden);
    [vertices_sorted, sortIdx] = sortrows(vertices_double);
    
    % Build inverse permutation mapping:
    invIdx = zeros(nV,1);
    for k = 1:nV
        invIdx(sortIdx(k)) = k;
    end
    
    % Remap adjacency according to sorted order.
    % We ensure each cell is a row vector.
    adj_sorted = cell(nV,1);
    for i = 1:nV
        newNbrs = arrayfun(@(x) invIdx(x), adjacency{i});
        newNbrs = sort(newNbrs); % sort for consistency
        adj_sorted{invIdx(i)} = newNbrs(:)'; % force row vector
    end
    
    fprintf('\n--- Sorted vertices (as floating point for testing) ---\n');
    disp(vertices_sorted);
    fprintf('--- Sorted adjacency ---\n');
    disp(adj_sorted);
    
    % --- Define expected sorted vertices (as double) and expected adjacency ---
    expected_vertices = [ 0.5,                0.333333333333333;
                          0.5,                0.75;
                          0.916666666666667,   0.333333333333333 ];
    % For a triangle, every vertex is adjacent to the other two.
    expected_adj = { [2 3], [1 3], [1 2] };
    
    tol = 1e-12;
    vertices_ok = all(abs(vertices_sorted - expected_vertices) < tol, 'all');
    if vertices_ok
        fprintf('TEST PASSED: Vertices match the expected fractional values.\n');
    else
        fprintf('TEST FAILED: Vertices do NOT match the expected values.\n');
        fprintf('Computed sorted vertices:\n');
        disp(vertices_sorted);
        fprintf('Expected sorted vertices:\n');
        disp(expected_vertices);
    end
    
    % Check each adjacency cell element by element.
    adjacency_ok = true;
    for i = 1:nV
        if ~isequal(adj_sorted{i}(:)', expected_adj{i})
            fprintf('Mismatch for vertex %d: computed [%s], expected [%s]\n', ...
                i, num2str(adj_sorted{i}(:)'), num2str(expected_adj{i}));
            adjacency_ok = false;
        end
    end
    if adjacency_ok
        fprintf('TEST PASSED: Adjacency matches the expected connectivity.\n');
    else
        fprintf('TEST FAILED: Adjacency does NOT match the expected connectivity.\n');
        fprintf('Computed sorted adjacency:\n');
        disp(adj_sorted);
        fprintf('Expected sorted adjacency:\n');
        disp(expected_adj);
    end
    
    % Overall test result: if both vertices and adjacency are correct.
    if vertices_ok && adjacency_ok
        fprintf('\nOVERALL TEST PASSED.\n');
    else
        fprintf('\nOVERALL TEST FAILED.\n');
    end
end
