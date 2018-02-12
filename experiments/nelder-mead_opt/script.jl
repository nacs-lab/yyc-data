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

function read_float64(prompt)
    while true
        print(prompt)
        str = readline()
        if isempty(str)
            throw(InterruptException())
        end
        try
            return parse(Float64, str)
        catch ex
            if isa(ex, InterruptException)
                rethrow()
            end
            println(ex)
            println("Please try again")
        end
    end
end

function print_states(io::IO, state::State{N}) where N
    for i in 1:N
        print(io, "x$i,")
    end
    println(io, "y")
    for i in 1:length(state.xs)
        x = state.xs[i]
        y = state.ys[i]
        for j in 1:N
            print(io, x[j], ",")
        end
        println(io, y)
    end
end

function interactive_opt(xs, ys, reverse=false)
    ccall(:jl_exit_on_sigint, Void, (Cint,), 0)
    state = State(1.0, 2.0, 0.5, 0.5, xs, ys)
    cb = function (x)
        rv = read_float64("Value for $x: ")
        return reverse ? -rv : rv
    end
    try
        while true
            iterate(state, cb)
            println("Iteration result:")
            print_states(STDOUT, state)
        end
    catch ex
        if !isa(ex, InterruptException)
            rethrow()
        end
        println("Optimization ends.")
    end
    sort!(state)
    println("\n\nFinal result:")
    print_states(STDOUT, state)
end

try
    @eval using Markdown
catch
    @eval using Base.Markdown
end

function usage()
    display(md"""
    Usage:

        script.jl <input_file>

    Interactively minimize a blackbox target function.

    The `input_file` should be a csv file that contains an optional one line header
    Each row of the csv file should contain the list of input variables followed by the
    result.
    The number of rows should be 1 larger than the number of input variables
    (i.e. the csv file should contain the same number of columns and rows).

    At each iteration, the script will ask for one or more evaluations of the target function
    at specific points and the results should be typed in as a floating point number.
    The result of the optimization will be printed out after each iteration which can be used
    to continue an optimization by placing it into the input csv file.

    Stop the optimization process with either `Ctrl-C` or `Ctrl-D`.
    The final result will be printed before the process exit.
    """)
end

if length(ARGS) < 1
    println("ERROR: no input csv file specified")
    usage()
    exit(1)
end
const fname = ARGS[1]

const data = try
    readdlm(fname, ',', Float64)
catch
    readdlm(fname, ',', Float64, skipstart=1)
end

if size(data, 1) != size(data, 2)
    println("ERROR: Incorrect input data dimension")
    usage()
    exit(1)
end
const xs = [ntuple(x->data[j, x], size(data, 2) - 1) for j in 1:size(data, 1)]
const ys = data[:, end]
interactive_opt(xs, ys, false)
