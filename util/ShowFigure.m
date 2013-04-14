% ShowFigure - Show a figure and set its name
%
% If the figure is already shown, the 'show' step is skipped, speeding
% up figure drawing by 10% - 40%.
%
% Usage:
%     ShowFigure(figureNum, figureName)
%     firstTime = ShowFigure(figureNum, figureName)
%
% Inputs:
%    figureNum      = the id of the figure to set
%    figureName     = the figure name to set
%
% Outputs:
%    firstTime      = true if the figure is being shown for the first
%                     time (optional).
%
% Created:   2/6/11, Paul King
%--------------------------------------------------------------------------
function firstTime_opt = ShowFigure( figureNum, figureName )

    % for speed, we only show the MatLab figure the first time through
    firstTime  = ~ ishandle(figureNum) || ~ strcmp(get(figureNum, 'Name'), figureName);

    if firstTime
        figure(figureNum);
        set(figureNum, 'Name',figureName, 'NumberTitle','off');

    elseif get(0,'CurrentFigure') ~= figureNum
        set(0, 'CurrentFigure', figureNum);
    end
    
    % return the firstTime flag if desired
    if nargout > 0
        firstTime_opt = firstTime;
    end
end
