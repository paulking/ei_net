% NetModel_Figure - Draw a figure showing model state
%
% Works by calling NetModel_Plot().
%
% Inputs:
%    model         = the simulation model with current running state
%    fig           = figure specification parameters:
%        figureNum     = figure number to use (optional)
%        title         = title (optional)
%
% Usage:
%    NetModel_Figure(model, fig);
%    NetModel_Figure(model, plotType, param1, value1, param2, value2, ...)
%
% Created:   2/5/2012, Paul King
%--------------------------------------------------------------------------
function fig = NetModel_Figure( model, fig, varargin )

    % special case: 'fig' is a string, so treat as plotType plus varargin
    if ischar(fig)
        fig = struct('plotType',fig);
        for i = 1:2:numel(varargin)
            fig.(varargin{i}) = varargin{i+1};
        end
    end

    % initialize and apply default parameter values
    defaultValues.figureNum   = [];
    defaultValues.title       = [];
    fig = ApplyDefaultValues(fig, defaultValues);

    % display figure and possibly assign a figureNum
    if isempty(fig.figureNum)
        fig.figureNum = figure();   % show a new figure
    elseif fig.figureNum < 1
        return;                     % disabled figure
    else
        figure(fig.figureNum);      % use the designated figure number
    end

    % draw the plot inside the figure window
    spikeHistory = [];
    if isfield(model, 'snapshot') && isfield(model.snapshot, 'spikeHistory')
        spikeHistory = model.snapshot.spikeHistory;
    end
    fig = NetModel_Plot(fig, model, spikeHistory);

    % add a title to the figure window
    set(fig.figureNum, 'Name',fig.title, 'NumberTitle','off');

    drawnow();
end
