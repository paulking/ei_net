% NetModel_FindElement - Find a cell group or cell group / input block combination
%
% Format for cellgroup is 'cg_xxx'.
%
% If i_ib return parameter is provided, then input block is expected
% and format is 'cg_xxx.in_yyy'.
%
% 'cg_output' is a special case that will return the model's defined
% output cell group model.outputCellGroupId.
%
% Inputs:
%    model         = the network model
%    name          = name to find
%
% Outputs:
%    i_cg          = the index of the cell group
%    i_ib          = the index of the input block within the cell group (optional)
%
% Created:   2/5/2012, Paul King
%--------------------------------------------------------------------------
function [i_cg, i_ib] = NetModel_FindElement( model, name )

    if nargout < 2
        assert( strncmp(name,'cg_',3), 'cell group name "%s" must start with "cg_"', name );

        % special case: if name is a cell array, then return an array of cell group IDs
        if iscell(name)
            i_cg = [];
            for i = 1:numel(name)
                i_cg(i) = NetModel_FindElement(model, name{i});
            end
            return;
        end

        i_cg = FindByName(model.cellGroup, name(4:end));
        if ~isempty(i_cg)
            return;
        end

        % special case: cell group 'output' refers to outputCellGroupId
        if strcmp(name, 'cg_output') && model.outputCellGroupId > 0
            i_cg = model.outputCellGroupId;
            return;
        end

        error('cell group "%s" not found', name);

    else
        dots = strfind(name, '.');
        assert( strncmp(name,'cg_',3), ...
            'input block name "%s" must start with cellgroup "cg_"', name);
        assert( ~isempty(dots) && (strncmp(name(dots(1)+1:end),'in_',3) ...
                || strncmp(name(dots(1)+1:end),'ib_',3)), ...
            'input block name "%s" must look like cg_xxx.in_xxx', name);

        i_cg = NetModel_FindElement(model, name(1:dots(1)-1));
        i_ib = FindByName(model.cellGroup{i_cg}.inputBlock, name(dots(1)+4:end));
        if isempty(i_ib)
            error('input block "%s" not found', name);
        end
    end
end
