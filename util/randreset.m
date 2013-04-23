% randreset - A quick way to reset all random number generation
function randreset()
    % note: MatLab changed the method to use around release 2011 (v 7.13 = R2011b)
    if verLessThan('matlab', '7.13')
        reset(RandStream.getDefaultStream());
    else
        reset(RandStream.getGlobalStream());
    end
end
