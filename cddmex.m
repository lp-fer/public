function varargout = cddmex(action, varargin)
% CDDMEX  MATLAB interface for the new cddmex MEX.
% 
% This function can receive rational numbers either as:
%   1) Symbolic values in fields: A, B, obj, V, ...
%   2) Split numerators/denominators in fields: Anum, Aden, Bnum, Bden, objNum, objDen, Vnum, Vden, ...
%
% Commands supported (matching "combined.c"):
%   solve_lp, solve_lp_ds,
%   adj_extreme, extreme, hull,
%   adjacency, copy_h, copy_v,
%   find_interior, find_interior_DS,
%   reduce_h, reduce_v, v_hull_extreme,
%   version
%
% Usage:
%    out = cddmex('solve_lp',  inStruct);
%    out = cddmex('adj_extreme', inStruct);
%    ...
%
% See also test_numeric.m, test_symbolic.m

try
    narginchk(1, 2); 
    if ~ischar(action) && ~isstring(action)
        error('cddmex:InvalidAction', 'ACTION must be a string.');
    end
    action = lower(char(action));
    
    switch action
        
        %--------------------------------------------
        % No real input needed, just prints version
        %--------------------------------------------
        case 'version'
            [varargout{1:nargout}] = mex_interface('version', struct(), struct());
            return;
        
        %--------------------------------------------
        % Commands that expect an LP structure (H plus objective)
        %--------------------------------------------
        case {'solve_lp', 'solve_lp_ds'}
            inStruct = getInputStruct(varargin);
            lpStruct = parseLP(inStruct);  % parse A,B,obj or numeric versions
            [varargout{1:nargout}] = mex_interface(action, lpStruct, inStruct);
            return;
        
        %--------------------------------------------
        % Commands that expect an H structure (without objective)
        %   (extreme, adj_extreme, reduce_h, copy_h, find_interior, find_interior_ds)
        %--------------------------------------------
        case {'extreme', 'adj_extreme', 'reduce_h', 'copy_h', ...
              'find_interior', 'find_interior_ds'}
            inStruct = getInputStruct(varargin);
            hStruct = parseH(inStruct);  % parse A,B (symbolic or numeric)
            [varargout{1:nargout}] = mex_interface(action, hStruct, inStruct);
            return;
        
        %--------------------------------------------
        % Commands that expect a V structure (for hull, reduce_v, copy_v, v_hull_extreme, adjacency)
        %--------------------------------------------
        case {'hull', 'reduce_v', 'copy_v', 'v_hull_extreme', 'adjacency'}
            inStruct = getInputStruct(varargin);
            vStruct = parseV(inStruct);  % parse V (and optional R) (symbolic or numeric)
            [varargout{1:nargout}] = mex_interface(action, vStruct, inStruct);
            return;
        
        otherwise
            error('cddmex:UnknownAction','Unknown action: %s', action);
    end
    
catch ME
    % Show usage info on error
    fprintf('Error: %s\n\n', ME.message);
    fprintf('Usage: cddmex(action, inputStruct)\n');
    fprintf('Possible actions: solve_lp, solve_lp_ds, extreme, adj_extreme, ...\n');
    rethrow(ME);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = getInputStruct(args)
    if isempty(args)
        s = struct();
    else
        s = args{1};
        if ~isstruct(s)
            error('cddmex:InvalidInput','Second argument must be a structure.');
        end
    end
end

function lpStruct = parseLP(inStruct)
    % Distinguish symbolic vs numeric by checking for 'A' vs 'Anum'
    if isfield(inStruct, 'A') && isfield(inStruct, 'B') && isfield(inStruct, 'obj')
        % Symbolic input: convert symbolic fields to split numerator/denominator arrays.
        [Anum,Aden] = arrayfun(@(x) numden(x), inStruct.A);
        [Bnum,Bden] = arrayfun(@(x) numden(x), inStruct.B);
        [objNum,objDen] = arrayfun(@(x) numden(x), inStruct.obj);
        lpStruct = struct('Anum', int64(Anum), 'Aden', int64(Aden), ...
                          'Bnum', int64(Bnum), 'Bden', int64(Bden), ...
                          'objNum', int64(objNum), 'objDen', int64(objDen));
    elseif isfield(inStruct, 'Anum') && isfield(inStruct, 'Aden') && ...
           isfield(inStruct, 'Bnum') && isfield(inStruct, 'Bden') && ...
           isfield(inStruct, 'objNum') && isfield(inStruct, 'objDen')
        % Numeric mode: already split
        lpStruct = inStruct;
    else
        error('cddmex:InvalidLP','LP struct must have either A,B,obj (symbolic) or Anum,Aden,Bnum,Bden,objNum,objDen (numeric).');
    end
end

function hStruct = parseH(inStruct)
    % H-struct: either symbolic (fields A and B) or numeric (Anum, Aden, Bnum, Bden)
    if isfield(inStruct, 'A') && isfield(inStruct, 'B')
        [Anum,Aden] = arrayfun(@(x) numden(x), inStruct.A);
        [Bnum,Bden] = arrayfun(@(x) numden(x), inStruct.B);
        hStruct = struct('Anum', int64(Anum), 'Aden', int64(Aden), ...
                         'Bnum', int64(Bnum), 'Bden', int64(Bden));
        if isfield(inStruct, 'lin')
            hStruct.lin = inStruct.lin;
        end
    elseif isfield(inStruct, 'Anum') && isfield(inStruct, 'Aden') && ...
           isfield(inStruct, 'Bnum') && isfield(inStruct, 'Bden')
        hStruct = inStruct;
    else
        error('cddmex:InvalidH','H struct must have either A,B (symbolic) or Anum,Aden,Bnum,Bden (numeric).');
    end
end

function vStruct = parseV(inStruct)
    % V-struct: either symbolic (field V) or numeric (fields Vnum, Vden, etc.)
    if isfield(inStruct, 'V') && ~isfield(inStruct, 'Vnum')
        [Vnum,Vden] = arrayfun(@(x) numden(x), inStruct.V);
        vStruct = struct('Vnum', int64(Vnum), 'Vden', int64(Vden));
        if isfield(inStruct, 'R')
            [Rnum,Rden] = arrayfun(@(x) numden(x), inStruct.R);
            vStruct.Rnum = int64(Rnum);
            vStruct.Rden = int64(Rden);
        end
        if isfield(inStruct, 'lin')
            vStruct.lin = inStruct.lin;
        end
    elseif isfield(inStruct, 'Vnum') && isfield(inStruct, 'Vden')
        vStruct = inStruct;
    else
        error('cddmex:InvalidV','V struct must have either "V" (symbolic) or "Vnum","Vden" (numeric).');
    end
end

function varargout = mex_interface(action, dataStruct, origStruct)
    % Handle adj_extreme separately because it explicitly requires two outputs
    if strcmp(action, 'adj_extreme')
        % Call the compiled MEX explicitly with two outputs
        [mexResult, adjacency] = cddmex_mex(action, dataStruct);

        % Ensure two outputs exist
        if isempty(mexResult) || isempty(adjacency)
            error('adj_extreme did not return expected results.');
        end

        % Convert numeric results back to symbolic if input was symbolic
        if isfield(origStruct, 'A')
            mexResult.V = mergeRationalArrayToSym(struct('num', mexResult.Vnum, 'den', mexResult.Vden));
            mexResult = rmfield(mexResult, {'Vnum', 'Vden'}); % Remove numeric fields
        end

        varargout{1} = mexResult;
        varargout{2} = adjacency;
        return;
    end

    % Call the compiled MEX for all other actions
    mexResult = cddmex_mex(action, dataStruct);

    switch action
        case {'solve_lp', 'solve_lp_ds', 'find_interior', 'find_interior_ds'}
            if ~isfield(origStruct, 'Anum')  % symbolic mode
                if isfield(mexResult, 'optvalue')
                    mexResult.optvalue = mergeRationalToSym(mexResult.optvalue);
                end
                if isfield(mexResult, 'sol')
                    mexResult.sol = mergeRationalArrayToSym(mexResult.sol);
                end
                if isfield(mexResult, 'dsol')
                    mexResult.dsol = mergeRationalArrayToSym(mexResult.dsol);
                end
            else
                % numeric mode: keep split fields
                if isfield(mexResult, 'optvalue')
                    mexResult.optvaluenum = mexResult.optvalue.num;
                    mexResult.optvalueden = mexResult.optvalue.den;
                end
                if isfield(mexResult, 'sol')
                    mexResult.solnum = mexResult.sol.num;
                    mexResult.solden = mexResult.sol.den;
                end
                if isfield(mexResult, 'dsol')
                    mexResult.dsolnum = mexResult.dsol.num;
                    mexResult.dsolden = mexResult.dsol.den;
                end
            end

        case {'extreme'}
            if isfield(origStruct, 'A')
                if isfield(mexResult, 'Vnum') && isfield(mexResult, 'Vden')
                    mexResult.V = mergeRationalArrayToSym(struct('num', mexResult.Vnum, 'den', mexResult.Vden));
                end
                if isfield(mexResult, 'Rnum') && isfield(mexResult, 'Rden')
                    mexResult.R = mergeRationalArrayToSym(struct('num', mexResult.Rnum, 'den', mexResult.Rden));
                end
            end

        case {'hull','reduce_v','copy_v','v_hull_extreme','adjacency'}
            if isfield(origStruct, 'V')
                if isfield(mexResult, 'Vnum') && isfield(mexResult, 'Vden')
                    mexResult.V = mergeRationalArrayToSym(struct('num', mexResult.Vnum, 'den', mexResult.Vden));
                end
                if isfield(mexResult, 'Rnum') && isfield(mexResult, 'Rden')
                    mexResult.R = mergeRationalArrayToSym(struct('num', mexResult.Rnum, 'den', mexResult.Rden));
                end
            end

        case {'reduce_h','copy_h'}
            if isfield(origStruct, 'A')
                if isfield(mexResult, 'Anum') && isfield(mexResult, 'Aden')
                    mexResult.A = mergeRationalArrayToSym(struct('num', mexResult.Anum, 'den', mexResult.Aden));
                end
                if isfield(mexResult, 'Bnum') && isfield(mexResult, 'Bden')
                    mexResult.B = mergeRationalArrayToSym(struct('num', mexResult.Bnum, 'den', mexResult.Bden));
                end
            end

        % Add any other command-specific post-processing here
    end

    varargout{1} = mexResult;
end

function sVal = mergeRationalToSym(ratStruct)
    % Merge a structure with fields 'num' and 'den' into a single symbolic value.
    if isstruct(ratStruct) && isfield(ratStruct,'num') && isfield(ratStruct,'den')
        sVal = sym(ratStruct.num) / sym(ratStruct.den);
    else
        sVal = ratStruct;
    end
end

function sArray = mergeRationalArrayToSym(ratStruct)
    % Merge arrays from fields 'num' and 'den' into a single symbolic array.
    if ~isstruct(ratStruct) || ~all(isfield(ratStruct, {'num','den'}))
        sArray = ratStruct;
        return;
    end
    numA = ratStruct.num;
    denA = ratStruct.den;
    if ~isequal(size(numA), size(denA))
        error('Size mismatch between numerator and denominator arrays.');
    end
    sArray = arrayfun(@(n,d) sym(n)/sym(d), numA, denA);
end
