% PlayMovie - Play a movie at a specified frame rate
%
% Play a movie one time through at a specified frame rate.
%
% This does what the MatLab command movie(mov, 1, fps) should
% do, without the extra display of movie frames as the movie
% is being loaded.
%
% The maximum frame rate on 2.7 GHz core2 duo at 320x240 is about 220/sec.
%
% Usage:
%    PlayMovie( mov, fps, show_frame_num )
%
% Inputs:
%    mov            = the movie to play (MatLab "movie" format)
%                     struct array(frames) with element cdata
%    fps            = the desired frames per second (optional, default = 15)
%    show_frame_num = display frame numbers during movie?
%
% Outputs:
%    [none]
%
% Created:   9/30/09, Paul King
%--------------------------------------------------------------------------
function PlayMovie( mov, fps, show_frame_num )

    % TODO: add support for additional movie types
    %     - cell array of grayscale uint8
    %     - cell array of grayscale double
    %     - 3D array of grayscale
    %     - cell array of RGB double

    % set default parameter values
    if ~exist('fps', 'var');             fps            = 15;       end
    if ~exist('show_frame_num', 'var');  show_frame_num = false;    end

    % play a MatLab "movie" frame by frame
    t0 = clock();
    for i = 1:numel(mov);

        % draw the movie frame
        if i == 1
            % show the first frame using the normal means
            h = imshow(mov(i).cdata);
        else
            % update subsequent frames via CData directly (much faster)
            set(h, 'CData', mov(i).cdata);
        end
        % ShowImageFast( mov(i).cdata );   % not quite as fast as above

        % display the frame number (optional)
        if show_frame_num
            title( sprintf('frame #%3d/%3d', i, numel(mov)) );
        end
        drawnow();

        % pause if necessary in order to achieve the desired frame rate
        elapsedTime = etime(clock(), t0);
        additionalWaitTime = i/fps - elapsedTime;
        if additionalWaitTime > 0
            pause(additionalWaitTime);
        end
    end
end

