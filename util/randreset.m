% randreset - A quick way to reset all random number generation
function randreset()
    reset(RandStream.getDefaultStream());
end
