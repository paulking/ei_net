% NetModel_EvalExpression - Evaluate an expression within a model context
%
% Evaluate a MatLab expression within a protected scope containing
% local variables for accessing model properties.
%
% If sourceName is provided, then 'cg' and possibly 'in' will be initialized
% to point to the cell group and (optionally) input block indicated
% by sourceName. sourceName can be of the format 'cg_xxx' or 'cg_xxx.in_yyy'.
% In this mode, the following local variables are initialized:
%    model  = the model
%    cg     = the cell group referenced by sourceName
%    in     = the input block referenced by sourceName (if any)
%
% If sourceName is not provided, then all cell groups and input blocks
% are initialized to be referencable as follows:
%    model         = the model
%    cg_xxx        = cell group named xxx
%    cg_xxx.in_yyy = input block named yyy of cell group xxx
%
% Inputs:
%    model      = model
%    expr       = the MatLab expression to evaluate (or a cell array of expressions)
%    sourceName = the cell group or input block serving as the context (optional)
%
% Outputs:
%    val     = the resulting value
%
% Created:   12/28/2012, Paul King
%--------------------------------------------------------------------------
function val = NetModel_EvalExpression(model, expr, sourceName_opt)

    % initialize referencable variables
    if nargin < 3
        % assign cg_xxx and cg_xxx.in_yyy local variables from model
        for i = 1:numel(model.cellGroup)
            cg = model.cellGroup{i};
            for j = 1:numel(model.cellGroup{i}.inputBlock)
                cg.(['in_' cg.inputBlock(j).name]) = cg.inputBlock(j);
            end
            eval(['cg_' cg.name ' = cg;']);
            clear i j cg;              % delete temporary variables from local scope
        end
    elseif ~isempty(strfind(sourceName_opt, '.in_')) || ~isempty(strfind(sourceName_opt, '.ib_'))
        % initialize 'cg' and 'in' from format 'cg_xxx.in_yyy'
        [i_cg, i_ib] = NetModel_FindElement(model, sourceName_opt);
        cg    = model.cellGroup{i_cg};
        in    = model.cellGroup{i_cg}.inputBlock(i_ib);
    elseif strncmp(sourceName_opt, 'cg_', 3)
        % initialize 'cg' from format 'cg_xxx'
        i_cg  = NetModel_FindElement(model, sourceName_opt);
        cg    = model.cellGroup{i_cg};
    else
        error('sourceName "%s" must have format cg_xxx or cg_xxx.in_yyy', sourceName_opt);
    end

    % evaluate the user-supplied expression(s) within the current scope
    if ischar(expr)
        val = eval(expr);
    elseif iscell(expr)
        for i = 1:numel(expr)
            val{i} = eval(expr{i});
        end
    else
        error('unknown type for expr');
    end
end

