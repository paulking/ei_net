% ShowImagePatchGrid - Display a matrix of 2D image patches (weights, bases)
%
% A fast routine to show a matrix of 2D image patches. Fast rendering is
% achieved by drawing directly into a single offscreen image buffer.
% The image is drawn into the current window.
%
% Usage:
%    ShowImagePatchGrid(patchData)
%    imageBuf = ShowImagePatchGrid(patchData)
%    imageBuf = ShowImagePatchGrid(patchData. spacing)
%
% Inputs:
%    patchData         = array(height,width,numPatches) containing
%                        basis functions (image patches) to render.
%    spacing           = spacing between patches (default = 3)
%
% Outputs:
%    imageBuf_out  = output image buffer (optional)
%                    If output argument provided, no image will be shown.
%
% Created:   8/15/2011, Paul King
%--------------------------------------------------------------------------
function imageBuf_out = ShowImagePatchGrid( patchData, spacing )

    if nargin >= 2
        s = spacing;
    else
        s = 2;        % # pixel spacing between basis image patches
    end
    
    % initialize image buffer to middle gray
    [height, width, numPatches] = size(patchData);
    subplotColumns = ceil( sqrt(numPatches) );
    subplotRows    = ceil( numPatches / subplotColumns );
    nonNegBases    = all(patchData(:) >= 0);

    imageBuffer = 0.8*ones(s+subplotRows*(height+s), s+subplotColumns*(width+s));

    k = 1;
    for c = 1:subplotColumns
        for r = 1:subplotRows
            if k > numPatches
                continue;    % in case the number of patches does not evenly fill a grid
            end
            clim = max(max(abs(patchData(:,:,k))));
            clim(clim == 0) = 1;     % to prevent division by zero
            if nonNegBases
                vals = patchData(:,:,k) / clim;
            else
                vals = patchData(:,:,k) / clim / 2 + .5;
            end            
            imageBuffer(s+(r-1)*(height+s)+[1:height], s+(c-1)*(width+s)+[1:width]) = vals;
            k = k + 1;
        end
    end

    % either display the patch matrix or return it
    if nargout <= 0
        imagesc(imageBuffer, [0 1]); axis image off; colormap gray;
        %ShowImageFast(imageBuffer);
    else
        imageBuf_out = imageBuffer;
    end
end
