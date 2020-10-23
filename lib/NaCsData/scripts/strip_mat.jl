#!/usr/bin/julia --startup-file=no

using MAT

function process_args()
    local iname, oname
    images = false
    counts = false
    i = 0
    maxnum = 0
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
        elseif arg == "--maxnum"
            i += 1
            maxnum = parse(Int, ARGS[i])
            @assert maxnum >= 0
        elseif arg == "--counts"
            counts = true
        else
            throw(ArgumentError("Unknown option: \"$arg\""))
        end
    end
    if !@isdefined(iname)
        throw(ArgumentError("Missing input file name"))
    elseif !@isdefined(oname)
        throw(ArgumentError("Missing output file name"))
    end
    if stat(oname).size == 0 && !endswith(basename(oname), ".mat")
        mkpath(oname, mode=0o755)
    end
    if isdir(oname)
        oname = joinpath(oname, basename(iname))
    end
    return (iname=iname, oname=oname, images=images, counts=counts, maxnum=maxnum)
end

function compute_counts(scan)
    imgs = scan["Images"]::Array{Float64,3}

    total_imgs = size(imgs, 3)
    nimgs::Int = scan["NumImages"]
    nsites::Int = scan["NumSites"]
    nseqs = total_imgs ÷ nimgs

    center_y = cld(size(imgs, 1), 2)
    center_x = cld(size(imgs, 2), 2)
    box_r::Int = cld(scan["BoxSize"] - 1, 2)

    counts = zeros(nimgs, nsites, nseqs)

    all_sites = scan["SingleAtomSites"]

    for se = 1:nseqs
        for i = 1:nimgs
            imgidx = (se - 1) * nimgs + i
            sites = all_sites[i]::Matrix{Float64}
            for si = 1:nsites
                site_x::Int = center_x + sites[si, 2]
                site_y::Int = center_y + sites[si, 1]
                v = 0.0
                for xi = site_x - box_r:site_x + box_r
                    for yi = site_y - box_r:site_y + box_r
                        v += imgs[xi, yi, imgidx]
                    end
                end
                counts[i, si, se] = v
            end
        end
    end
    return counts
end

const opts = process_args()

# Work around MAT.jl bug
function clear_char!(dict::Dict)
    for k in keys(dict)
        v = dict[k]
        if isa(v, Char)
            dict[k] = string(v)
        else
            clear_char!(v)
        end
    end
end
function clear_char!(ary::Array)
    for i in 1:length(ary)
        v = ary[i]
        if isa(v, Char)
            ary[i] = string(v)
        else
            clear_char!(v)
        end
    end
end
function clear_char!(o)
end

matopen(opts.iname) do mf
    pl = read(mf, "ParamList")
    sa = read(mf, "Analysis")["SingleAtomLogical"]
    scan = read(mf, "Scan")
    param_name = scan["ParamName"]
    if isa(param_name, Char)
        param_name = string(param_name)
    end
    sg = haskey(scan, "ScanGroup") ? scan["ScanGroup"] : Dict{String,Any}()
    clear_char!(sg)
    if eltype(sa) == Bool
        sa = UInt8.(sa)
    end

    nseq = size(pl, 2)
    if opts.maxnum > 0 && opts.maxnum < nseq
        pl = pl[:, 1:opts.maxnum]
        sa = sa[:, :, 1:opts.maxnum]
        if opts.images || opts.counts
            nimgs::Int = scan["NumImages"]
            scan["Images"] = [:, :, 1:(opts.maxnum * nimgs)]
        end
    end

    matopen(opts.oname, "w") do out
        write(out, "ParamList", pl)
        write(out, "SingleAtomLogical", sa)
        write(out, "ParamName", param_name)
        write(out, "ScanGroup", sg)
        if opts.images
            write(out, "Images", scan["Images"])
        end
        if opts.counts
            write(out, "Counts", compute_counts(scan))
        end
    end
end
