#

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

function select_single(loads, survives)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        for l in loads
            if l > 0 && logicals[l] == 0
                return [1, 0, 0]
            elseif l < 0 && logicals[-l] != 0
                return [1, 0, 0]
            elseif l == 0
                throw(ArgumentError("Load logical cannot be 0."))
            end
        end
        for s in survives
            if s > 0 && logicals[s] == 0
                return [1, 1, 0]
            elseif s < 0 && logicals[-s] != 0
                return [1, 1, 0]
            elseif s == 0
                throw(ArgumentError("Survival logical cannot be 0."))
            end
        end
        return [1, 1, 1]
    end
end
