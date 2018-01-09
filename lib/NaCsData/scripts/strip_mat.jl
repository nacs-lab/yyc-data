#!/usr/bin/julia --startup-file=no

using MAT

function process_args()
    local iname, oname
    images = false
    i = 0
    while i < length(ARGS)
        i += 1
        arg = ARGS[i]
        if arg[1] != '-'
            if !@isdefined(iname)
                iname = arg
            elseif !@isdefined(oname)
                oname = arg
            else
                throw(ArgumentError("Extra argument: \"$arg\""))
            end
            continue
        end
        if arg == "--images"
            images = true
        else
            throw(ArgumentError("Unknown option: \"$arg\""))
        end
    end
    if !@isdefined(iname)
        throw(ArgumentError("Missing input file name"))
    elseif !@isdefined(oname)
        throw(ArgumentError("Missing output file name"))
    end
    if stat(oname).size == 0
        mkpath(oname, 0o755)
    end
    if isdir(oname)
        oname = joinpath(oname, "$(splitext(basename(iname))[1]).mat")
    end
    return (iname=iname, oname=oname, images=images)
end

const opts = process_args()

matopen(opts.iname) do mf
    pl = read(mf, "ParamList")
    sa = read(mf, "Analysis")["SingleAtomLogical"]
    scan = read(mf, "Scan")
    param_name = scan["ParamName"]

    matopen(opts.oname, "w") do out
        write(out, "ParamList", pl)
        write(out, "SingleAtomLogical", sa)
        write(out, "ParamName", param_name)
        if opts.images
            write(out, "Images", scan["Images"])
        end
    end
end
