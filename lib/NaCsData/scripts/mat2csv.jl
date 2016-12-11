#!/usr/bin/julia -f

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using NaCsData

iname = ARGS[1]
oname = ARGS[2]
param_name = length(ARGS) >= 3 ? ARGS[3] : ""

_param_name, params, data = NaCsData.load_matscan(iname)
isempty(param_name) && (param_name = _param_name)

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
