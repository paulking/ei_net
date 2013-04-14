% FindByName - Find an element in an array by name
%
% This function will search a struct array or a cell array of structs
% for an entry whose 'name' field has the value specified.  If
% no such entry can be found, [] is returned.  If there are multiple
% entries with the same name, only the first one is returned.
%
% Usage:
%    index = FindByName( array, name )
%
% Inputs:
%    array   = the array of cells or struct array to search
%    name    = the name to search for
%
% Output:
%    index   = the index where the item was found or []
%
% Created:   11/21/09, Paul King
%--------------------------------------------------------------------------
function index = FindByName( array, name )

    index = [];
    
    % search a cell array
    if iscell(array)
        for i = 1:numel(array)
            if strcmp(array{i}.name, name)
                index = i;
                return;
            end
        end

    % search a struct array
    elseif isstruct(array)
        for i = 1:length(array)
            if strcmp(array(i).name, name)
                index = i;
                return;
            end
        end
        % index = min(find(strcmp(name, {array.name}))); % slightly faster

    % no elements to search
    elseif isempty(array)
        return;                       % return null (the default value)

    else
        error('unsupported array type');
    end
end

