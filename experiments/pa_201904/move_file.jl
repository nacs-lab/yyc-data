#!/usr/bin/julia

using MAT

names = matopen(ARGS[1]) do fd
    read(fd, "names")
end

from_dir = ARGS[2]
if length(ARGS) > 2
    to_dir = ARGS[3]
else
    to_dir = joinpath(@__DIR__, "data")
end

for name in names
    println("Moving $name")
    from = joinpath(from_dir, "data_$name.mat")
    to = joinpath(to_dir, "data_$name.mat")
    mv(from, to)
end
