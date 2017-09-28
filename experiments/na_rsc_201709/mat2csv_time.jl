#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using MAT

iname = ARGS[1]
oname = ARGS[2]

if isdir(oname)
    oname = joinpath(oname, "$(splitext(basename(iname))[1]).csv")
end

const mf = matopen(iname)
const pl = read(mf, "ParamList")
const sa = read(mf, "Analysis")["SingleAtomLogical"]
const param_name = read(mf, "Scan")["ParamName"]

if size(sa, 2) != 1
    println(STDERR, "Multiple site not supported")
    exit(1)
end

open(oname, "w") do fd
    print(fd, param_name)
    for j in 1:size(sa, 1)
        write(fd, ",Image$j")
    end
    println(fd)
    for i in 1:length(pl)
        print(fd, pl[i], ",1")
        for j in 1:size(sa, 1)
            print(fd, ',', sa[j, 1, i])
        end
        println(fd)
    end
end
