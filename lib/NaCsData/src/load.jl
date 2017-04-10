#!/usr/bin/julia -f

using MAT
import NaCsCalc.Utils: binomial_estimate

"""
Load a MAT scan file
"""
function load_matscan(fnames)
    # FIXME: https://github.com/simonster/MAT.jl/issues/53
    param_name = "Parameter"
    data_sorter = Dict{Float64,Vector{Int}}()
    glob_numimages = Ref(0)
    for fname in fnames
        matopen(fname) do fd
            # * Values we want to save:
            #   * parameter name (error due to MAT.jl)
            #   * parameter values
            #   * load numbers
            # * Values we need to read:
            #   * single atom
            #      * single atom
            #      * count + (input) threshold
            #      * image + (input) ROI + (input?) threshold
            #   * parameter list
            #   * number of images per sequence
            local numimages = round(Int, read(fd, "memmap")["NumImages"])::Int
            if glob_numimages[] == 0
                glob_numimages[] = numimages
            elseif glob_numimages[] != numimages
                error("Image number mismatch. ($(glob_numimages[]) != $numimages)")
            end
            local single_atom = read(fd, "SingleAtom")::Matrix{Float64}
            local param_list = read(fd, "ParamList")::Matrix{Float64}
            for i in 1:length(param_list)
                @inbounds param = param_list[i]
                frame = if param in keys(data_sorter)
                    data_sorter[param]
                else
                    data_sorter[param] = zeros(Int, numimages + 1)
                end
                frame[1] += 1
                for j in 1:numimages
                    # Ignore "single atoms" in second image but not the first one.
                    # This should probably be an option/filter
                    if single_atom[numimages * (i - 1) + j] == 0
                        break
                    end
                    frame[j + 1] += 1
                end
            end
        end
    end
    numimages = glob_numimages[]
    params = sort(collect(keys(data_sorter)))
    data = Matrix{Int}(numimages + 1, length(params))
    for i in 1:length(params)
        frame = data_sorter[params[i]]
        for j in 1:(numimages + 1)
            data[j, i] = frame[j]
        end
    end
    param_name, params, data
end
load_matscan(fname::AbstractString) = load_matscan([fname])

function calc_survival(fnames)
    data_dict = Dict{Float64,Vector{Float64}}()
    local num_cnts::Int
    for fname in fnames
        data = readcsv(fname, Float64, skipstart=1)
        num_cnts = size(data, 2) - 1
        for i in 1:size(data, 1)
            param = data[i, 1]
            if param in keys(data_dict)
                frame = data_dict[param]
                for j in 1:num_cnts
                    frame[j] += data[i, j + 1]
                end
            else
                data_dict[param] = data[i, 2:end]
            end
        end
    end
    params = sort(collect(keys(data_dict)))
    len = length(params)
    ratios = Matrix{Float64}(len, num_cnts - 1)
    uncs = Matrix{Float64}(len, num_cnts - 1)
    for i in 1:len
        frame = data_dict[params[i]]
        base = frame[1]
        for j in 1:(num_cnts - 1)
            cur = frame[j + 1]
            ratios[i, j], uncs[i, j] = binomial_estimate(cur, base)
            base = cur
        end
    end
    params, ratios, uncs
end
calc_survival(fnames::AbstractString) = calc_survival([fnames])
