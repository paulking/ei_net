% ApplyDefaultValues - Apply default values to a parameter struct
%
% The parameter struct is overlayed on top of the default values struct,
% The field order of the defaultValues struct is preserved, but any field
% present in the param struct overwrites the default Value.
%
% If multiple default value structs are provided, they will be applied
% in order (as defaults to the defaults, etc.)
%
% If scalar parameters are provided, they will be handled reasonably.
%
% If 'underlay' is indicated as the 3rd parameter, default values
% will be added to the base params struct rather than the other way around,
% which has the net effect of preserving the field order of the param struct
% rather than the default value struct.
%
% Usage:
%   params = ApplyDefaultValues(params, defaultValues);
%   params = ApplyDefaultValues(params, defaultValues, defaultValues2, ...);
%   params = ApplyDefaultValues(params, defaultValues, 'underlay');
%
% Inputs:
%   params        = an input parameter struct
%   defaultValues = default values in a parallel struct to apply
%                    
% Output:
%   params        = the updated parameter struct

% Created:   2/6/11, Paul King
%--------------------------------------------------------------------------
function params_out = ApplyDefaultValues( params, defaultValues, varargin )

    % special case: params can be a cell array of alternating field names and values
    if iscell(params)
        newParams = struct();
        for i = 1:2:numel(params)
            newParams.(params{i}) = params{i+1};
        end
        params = newParams;
    end

    if numel(varargin) > 0
        % special case: default value 'underlay' requested
        if numel(varargin) == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'underlay')
            params_out = ApplyDefaultValues_Underlay(params, defaultValues);
            return;
        end
        
        % otherwise, if there are any extra arguments, apply them as backup defaults
        for k = 1:numel(varargin)
            defaultValues = ApplyDefaultValues(defaultValues, varargin{k});
        end
    end

    % special case: params are scalar, so do simple substitution
    if ~isstruct(defaultValues)
        if ~isempty(params)
            params_out = params;
        else
            params_out = defaultValues;
        end
        return;
    end

    % get the names of all fields in the params struct
    nameList = fieldnames(params);
    
    % apply each param individually to params_out
    params_out = defaultValues;
    for i = 1:numel(nameList)
        name = nameList{i};
        if isstruct(params.(name)) && isfield(params_out, name)
            % recursively apply subparameter values
            params_out.(name) = ApplyDefaultValues(params.(name), params_out.(name));
        else
            params_out.(name) = params.(name);
        end
    end
end


% ApplyDefaultValues_Underlay - Apply default values via 'underlay' method
%
% Same as above, but add missing defaults to the params struct (underlay)
% rather than adding extra params to the defaultValues struct (overlay).
%
% Usage:
%   params = ApplyDefaultValues_Underlay(params, defaultValues);
%   params = ApplyDefaultValues_Underlay(params, defaultValues, defaultValues2, ...);
%
% Inputs:
%   params        = an input parameter struct
%   defaultValues = default values in a parallel struct to apply
%                    
% Output:
%   params        = the updated parameter struct
function params = ApplyDefaultValues_Underlay( params, defaultValues, varargin )

    % if there are any extra arguments, apply them as backup defaults
    for k = 1:numel(varargin)
        defaultValues = ApplyDefaultValues_Underlay(defaultValues, varargin{k});
    end

    % get the names of all fields in the defaultValues struct
    nameList = fieldnames(defaultValues);
    
    % apply each default field individually to params
    for i = 1:numel(nameList)
        name = nameList{i};
        if ~ isfield(params, name)
            params.(name) = defaultValues.(name);
        elseif isstruct(defaultValues.(name))

            % recursively apply subparameter values, if any
            params.(name) = ApplyDefaultValues_Underlay(params.(name), defaultValues.(name));
        end
    end
end
