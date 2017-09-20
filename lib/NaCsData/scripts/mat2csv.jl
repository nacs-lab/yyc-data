#!/usr/bin/julia --startup-file=no

push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using NaCsData

iname = ARGS[1:end - 1]
oname = ARGS[end]

if length(iname) == 0
    println(STDERR, "No input file given")
    exit(1)
end

if isdir(oname)
    oname = joinpath(oname, "$(splitext(basename(iname[1]))[1]).csv")
end

param_name, params, data = NaCsData.load_matscan(iname)

open(oname, "w") do fd
    print(fd, param_name)
    for j in 1:size(data, 1)
        write(fd, ",Image$j")
    end
    println(fd)
    for i in 1:length(params)
        print(fd, params[i])
        for j in 1:size(data, 1)
            write(fd, ",")
            print(fd, data[j, i])
        end
        println(fd)
    end
end
