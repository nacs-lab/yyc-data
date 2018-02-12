#!/usr/bin/julia

struct State{N}
    α::Float64
    γ::Float64
    ρ::Float64
    σ::Float64
    xs::Vector{NTuple{N,Float64}}
    ys::Vector{Float64}
end

Base.issorted(state::State) = issorted(state.ys)

function Base.sort!(state::State)
    issorted(state) && return
    perm = sortperm(state.ys)
    state.ys .= state.ys[perm]
    state.xs .= state.xs[perm]
    return
end

function iterate(state::State{N}, cb) where N
    sort!(state)
    # Centroid
    xo = ntuple(x->0.0, Val{N}())
    xs = state.xs
    n = length(xs) - 1
    for i in 1:n
        xo = xo .+ xs[i]
    end
    xo = xo ./ n
    # Reflection
    xr = xo .+ state.α .* (xo .- xs[n + 1])
    fr = cb(xr)
    if fr < state.ys[n]
        state.ys[n + 1] = fr
        xs[n + 1] = xr
        if state.ys[1] <= fr
            return
        end
        # Expansion
        xe = xo .+ state.γ .* (xr .- xo)
        fe = cb(xe)
        if fe < fr
            state.ys[n + 1] = fe
            xs[n + 1] = xe
        end
        return
    end
    # Contraction
    xc = xo .+ state.ρ .* (xs[n + 1] .- xo)
    fc = cb(xc)
    if fc < state.ys[n + 1]
        state.ys[n + 1] = fc
        xs[n + 1] = xc
        return
    end
    # Shrink
    for i in 2:(n + 1)
        x = xs[i]
        x′ = x[1] .+ state.σ .* (x .- x[1])
        f′ = cb(x′)
        state.ys[i] = f′
        xs[i] = x′
    end
end

# function test_cb(f, n, xs::Vector{NTuple{N,Float64}}) where N
#     xs = copy(xs)
#     ys = Float64[f(x...) for x in xs]
#     state = State{N}(1.0, 2.0, 0.5, 0.5, copy(xs), ys)
#     cb = function (x)
#         push!(xs, x)
#         return f(x...)
#     end
#     for i in 1:n
#         iterate(state, cb)
#     end
#     return state.xs[1], state.ys[1], xs
# end

# res = test_cb((x, y) -> (x - 3)^2 + 3 * (x - 3) * (y + 2) + 4 * (y + 2)^2,
#               100, [(1.0, 3.0), (2.4, 9.0), (-2.0, 2.0)])
# @show res
