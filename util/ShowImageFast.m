% ShowImageFast - Show an image using fast direct data access if possible
%
% This routine is especially useful for playing video or animation
% sequences with low overhead, when the same image frame is updated
% repeatedly.
%
% The first time this function is called, imshow() is used. So there is no
% performance gain on one-time calls.
%
% Notes:
%    - this must be called on a blank (or identically created) figure to work   
%    - drawnow() must be called following this function or the figure will
%      not be updated.
%
% Usage:
%    ShowImageFast(imageSrc)
%
% Inputs:
%    imageSrc   = image frame, type = double.  Can be either grayscale
%                 (a 2D matrix) or an RGB image (a 3D matrix with 3 channels).
%
% Created:   9/30/09, Paul King
%--------------------------------------------------------------------------
function ShowImageFast( imageSrc )

    % get the handle to the image, if one exists
    h = get(gca(), 'Children');

    % if this is the first call, show it using the normal means
    if isempty(h)
        imshow(imageSrc);
        return;
    end

    % update the frame directly via CData (much faster)
    set(h, 'CData', imageSrc);
end
