#!/usr/bin/julia -f

using SingleAtom

function run_test(fname)
    println("Testing $fname")
    include(fname)
end

run_test("utils.jl")
run_test("optical.jl")
run_test("atomic.jl")
run_test("system.jl")
run_test("propagate.jl")
