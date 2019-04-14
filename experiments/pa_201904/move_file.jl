#!/usr/bin/julia

using MAT

names = matopen(ARGS[1]) do fd
    read(fd, "names")
end

for name in names
    println("Moving $name")
    from = joinpath(ARGS[2], "data_$name.mat")
    to = joinpath(@__DIR__, "data", "data_$name.mat")
    mv(from, to)
end
