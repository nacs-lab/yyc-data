#!/usr/bin/julia -f

using MAT
import NaCsCalc.Utils: binomial_estimate
import NaCsCalc.Format: Unc
using DataStructures
using DelimitedFiles

abstract type AbstractValues{N} end

@inline to_arrayidx(idx::Integer) = idx:idx
@inline to_arrayidx(idx) = idx

function depth end
function combiner_type end
function combine end
function create_values end

struct SortedData{N,N2,K,Vs<:AbstractValues{N2}}
    params::Array{K,N}
    values::Vs
    function SortedData{N,N2,K,Vs}(params, values) where {N,N2,K,Vs}
        @assert N + 1 == N2
        return new(params, values)
    end
end
SortedData(params::AbstractArray{K,N}, values::Vs) where {N,N2,K,Vs<:AbstractValues{N2}} =
    SortedData{N,N2,K,Vs}(params, values)

SortedData1{K,Vs} = SortedData{1,2,K,Vs}

@generated function Base.getindex(data::SortedData{N,N2,K,Vs}, _args...) where {N,N2,K,Vs}
    nargs = length(_args)
    if N + 1 != N2
        return :(throw(ArgumentError($"Inconsistent type $(data)")))
    end
    if nargs > N2
        return :(throw(BoundsError(data, _args)))
    end
    if nargs == 0
        return :data
    end
    args = [Symbol("arg$i") for i in 1:N2]
    quote
        $(Expr(:meta, :inline))
        $((:($(args[i]) = to_arrayidx(_args[$i])) for i in 1:nargs)...)
        $((:($(args[i]) = :) for i in (nargs + 1):N2)...)
        return SortedData{N,N2,K,Vs}(data.params[$((args[i] for i in 1:N)...)],
                                     data.values[$((args[i] for i in 1:N2)...)])
    end
end

@inline Base.size(data::SortedData) = size(data.params)..., depth(data.values)
@inline Base.endof(data::SortedData) = endof(data.params)
@inline Base.size(data::SortedData{N}, dim) where {N} = if dim <= N
    return size(data.params, dim)
elseif dim <= N + 1
    return depth(data.values)
else
    return 1
end
@inline Base.ndims(::SortedData{N,N2}) where {N,N2} = N2
function Base.vcat(datas::SortedData1{K,Vs}...) where {K,Vs}
    # TODO
    CT = combiner_type(Vs)
    combiners = Dict{K,CT}()
    params = Vector{K}()
    for data in datas
        ps = data.params
        vals = data.values
        nps = length(ps)
        for i in 1:nps
            p = ps[i]
            if haskey(combiners, p)
                combine(combiners[p], vals, i)
            else
                push!(params, p)
                combiners[p] = CT(vals, i)
            end
        end
    end
    SortedData1{K,Vs}(params, create_values(params, combiners))
end
function Base.show(io::IO, data::SortedData)
    params, ratios, uncs = get_values(data)
    print(io, "SortedData", "(params=", params, ", data=", Unc.(ratios, uncs), ")")
end
function dump_raw(io::IO, data::SortedData{1})
    for i in 1:size(data, 1)
        print(io, data.params[i])
        for j in 1:size(data, 2)
            print(io, ',', data.values[i, j])
        end
        println(io)
    end
end
dump_raw(data::SortedData) = dump_raw(STDOUT, data)
dump_raw(fname::AbstractString, data::SortedData) = open(fname, "w") do io
    dump_raw(io, data)
end

struct CountValues{N2} <: AbstractValues{N2}
    counts::Array{Int,N2}
end
CountData{N,N2,K} = SortedData{N,N2,K,CountValues{N2}}
CountData(params::AbstractArray{T,N}, counts::AbstractArray{T2,N2}) where {N,N2,T,T2} =
    CountData{N,N2,T}(params, CountValues{N2}(counts))
@inline depth(vals::CountValues{N2}) where {N2} = size(vals.counts, N2)
@inline _maybe_countvalues(v::Number) = v
@inline _maybe_countvalues(v::AbstractArray{T,N}) where {T,N} = CountValues{N}(v)
@inline Base.getindex(vals::CountValues, args...) = _maybe_countvalues(vals.counts[args...])

function load_count_csv(fname)
    data = readdlm(fname, ',', Float64, skipstart=1)
    params = data[:, 1]
    counts = Int.(@view data[:, 2:end])
    return CountData(params, counts)
end

# TODO
const CountValues1 = CountValues{2}
CountData1{K} = CountData{1,2,K}
struct CountCombiner
    counts::Vector{Int}
end
combiner_type(::Type{CountValues1}) = CountCombiner
CountCombiner(vals::CountValues1, i) = CountCombiner(vals.counts[i, :])
function combine(comb::CountCombiner, vals::CountValues1, i)
    for j in 1:length(comb.counts)
        comb.counts[j] += vals.counts[i, j]
    end
end
function create_values(params, dict::Dict{<:Any,CountCombiner})
    nparams = length(params)
    @assert nparams > 0
    c0 = dict[params[1]].counts
    ncnts = length(c0)
    counts = Matrix{Int}(nparams, ncnts)
    for i in 1:nparams
        c = dict[params[i]].counts
        for j in 1:ncnts
            counts[i, j] = c[j]
        end
    end
    return CountValues1(counts)
end

struct SurvivalValues{N2} <: AbstractValues{N2}
    ratios::Array{Float64,N2}
    uncs::Array{Float64,N2}
end
function SurvivalValues(vals::CountValues{N2}) where {N2}
    counts = vals.counts
    csize = size(counts)
    psize = ntuple(i->csize[i], N2 - 1)
    num_cnts = csize[N2]
    ratios = Array{Float64,N2}(psize..., num_cnts - 1)
    uncs = Array{Float64,N2}(psize..., num_cnts - 1)
    for I in CartesianRange(psize)
        base = counts[I.I..., 1]
        for j in 1:(num_cnts - 1)
            cur = counts[I.I..., j + 1]
            ratios[I.I..., j], uncs[I.I..., j] = binomial_estimate(cur, base)
            base = cur
        end
    end
    return SurvivalValues{N2}(ratios, uncs)
end
@inline _maybe_survivalvalues(v1::Number, v2::Number) = Unc(v1, v2)
@inline function _maybe_survivalvalues(v1::AbstractArray{T1,N},
                                       v2::AbstractArray{T2,N}) where {T1,T2,N}
    return SurvivalValues{N}(v1, v2)
end
@inline function Base.getindex(vals::SurvivalValues, args...)
    return _maybe_survivalvalues(vals.ratios[args...], vals.uncs[args...])
end
@inline depth(vals::SurvivalValues{N2}) where {N2} = size(vals.ratios, N2)

# TODO
const SurvivalValues1 = SurvivalValues{2}
struct SurvivalCombiner
    sums::Vector{Float64}
    ws::Vector{Float64}
end
combiner_type(::Type{SurvivalValues1}) = SurvivalCombiner
function SurvivalCombiner(vals::SurvivalValues1, i)
    ratios = vals.ratios
    uncs = vals.uncs
    ncnts = size(ratios, 2)
    sums = Vector{Float64}(ncnts)
    ws = Vector{Float64}(ncnts)
    for j in 1:ncnts
        w = 1 / uncs[i, j]^2
        ws[j] = w
        sums[j] = ratios[i, j] * w
    end
    SurvivalCombiner(sums, ws)
end
function combine(comb::SurvivalCombiner, vals::SurvivalValues1, i)
    ratios = vals.ratios
    uncs = vals.uncs
    ncnts = size(ratios, 2)
    sums = comb.sums
    ws = comb.ws
    for j in 1:ncnts
        w = 1 / uncs[i, j]^2
        ws[j] += w
        sums[j] += ratios[i, j] * w
    end
end
function create_values(params, dict::Dict{<:Any,SurvivalCombiner})
    nparams = length(params)
    @assert nparams > 0
    s0 = dict[params[1]].sums
    ncnts = length(s0)
    ratios = Matrix{Float64}(nparams, ncnts)
    uncs = Matrix{Float64}(nparams, ncnts)
    for i in 1:nparams
        v = dict[params[i]]
        ss = v.sums
        ws = v.ws
        for j in 1:ncnts
            s = ss[j]
            w = ws[j]
            ratios[i, j] = s / w
            uncs[i, j] = 1 / √(w)
        end
    end
    return SurvivalValues1(ratios, uncs)
end

SurvivalData{N,N2,K} = SortedData{N,N2,K,SurvivalValues{N2}}
SurvivalData(data::CountData{N,N2,K}) where {N,N2,K} =
    SurvivalData{N,N2,K}(data.params, SurvivalValues(data.values))
get_values(data::SurvivalData) = data.params, data.values.ratios, data.values.uncs
get_values(data::CountData) = get_values(SurvivalData(data))

function map_params(f::F, data::SortedData) where {F}
    params = data.params
    rng = CartesianRange(size(params))
    SortedData([f(idx.I..., params[idx.I...]) for idx in rng], data.values)
end
map_params(f::F, data::Tuple) where {F} = map(d->map_params(f, d), data)
map_params(f::F, data::OrderedDict) where {F} = OrderedDict(k=>map_params(f, v) for (k, v) in data)

function _split_data(data, offset, dict, spec::OrderedDict)
    return OrderedDict(k=>_split_data(data, offset, dict, v) for (k, v) in spec)
end

function _split_data(data, offset, dict, spec::Tuple)
    return map(s->_split_data(data, offset, dict, s), spec)
end

function _split_data(data, _offset, dict, spec::AbstractArray{T}) where {T}
    nspec = length(spec)
    params = Vector{T}()
    idxs = Vector{Int}()
    offset = _offset[]
    _offset[] = offset + nspec
    for i in 1:nspec
        if !haskey(dict, i + offset)
            continue
        end
        push!(idxs, dict[i + offset])
        push!(params, spec[i])
    end
    return SortedData(params, data.values[idxs, :])
end

function split_data(_data::SortedData1, spec)
    # Remove duplicates (not sure if it's an better API to do this implicitly or explicitly)
    data = [_data;]
    params = data.params

    # Create a map from parameter value
    dict = Dict{Int,Int}()
    nparams = length(params)
    for i in 1:nparams
        dict[params[i]] = i
    end

    offset = Ref(0)
    return _split_data(data, offset, dict, spec)
end

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
            local single_atom
            local param_list
            param_list = read(fd, "ParamList")::Matrix{Float64}
            if exists(fd, "SingleAtom")
                # Old format
                single_atom = vec(read(fd, "SingleAtom")::Matrix{Float64})
            else
                # New multi atom format, only support one site for now
                analysis = read(fd, "Analysis")::Dict{String,Any}
                sal = analysis["SingleAtomLogical"]::Array{Float64,3}
                if size(sal, 2) != 1
                    error("Multi-cite sequence not supported yet.")
                end
                single_atom = vec(sal)
            end
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

function calc_survival(datas)
    data_dict = Dict{Float64,Vector{Float64}}()
    local num_cnts::Int
    for data::AbstractMatrix in datas
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
calc_survival(fnames::AbstractVector{T} where T<:AbstractString) =
    calc_survival(readdlm(fname, ',', Float64, skipstart=1) for fname in fnames)
calc_survival(fnames::AbstractString) = calc_survival([fnames])

function load_striped_mat(fname)
    matopen(fname) do fd
        params = read(fd, "ParamList")
        @assert size(params, 1) == 1
        return params[1, :], read(fd, "SingleAtomLogical") .!= 0
    end
end

_seqid_iter(params, maxseq::Integer) = 1:min(length(params), maxseq)
_seqid_iter(params, filter) = (i for i in 1:length(params) if filter(i))

function select_count(_params::AbstractVector{T}, logicals::AbstractArray{TL,3} where TL,
                      selector, filter=length(_params)) where T
    num_cnt::Int = 0
    data_dict = Dict{T,Vector{Int}}()
    for seq in _seqid_iter(_params, filter)
        param = _params[seq]
        counts = selector(@view logicals[:, :, seq])
        if num_cnt == 0
            num_cnt = length(counts)
            if num_cnt == 0
                throw(ArgumentError("No counts returned by the selector"))
            end
        elseif num_cnt != length(counts)
            throw(ArgumentError("Selector does not return the same number of elements"))
        end
        if haskey(data_dict, param)
            data_dict[param] .+= counts
        else
            data_dict[param] = counts
        end
    end
    params = sort(collect(keys(data_dict)))
    len = length(params)
    counts = zeros(Int, len, num_cnt)
    for i in 1:len
        counts[i, :] = data_dict[params[i]]
    end
    return CountData(params, counts)
end
