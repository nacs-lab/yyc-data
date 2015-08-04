#!/usr/bin/julia -f

using SingleAtom

@generated function test_scale{N,T}(nele, factors::NTuple{N,Complex{T}},
                                    nloop=10_000)
    splat_factors = [gensym(:splat) for i in 1:N]
    init_ex = quote
        println("nele: $nele, ndim: $($N)")
        ary = Matrix{Complex{T}}(nele, $N)
        soary = convert(StructOfArrays, ary)
        factor2 = exp(rand(T, nele) * im)
        sofactor2 = convert(StructOfArrays, factor2)
        $([:($(splat_factors[j]) = factors[$j]) for j in 1:N]...)
    end
    rand_ex = quote
        rand!(ary)
        unsafe_copy!(soary, ary)
    end
    var_names = [gensym(:ele) for i in 1:N]
    eff_factors = [gensym(:factor) for i in 1:N]
    aos_run = quote
        @inbounds for i in 1:nele
            $([:($(var_names[j]) = ary[i, $j]) for j in 1:N]...)
            factor2i = factor2[i]
            $([:($(eff_factors[j]) = factor2i * $(splat_factors[j]))
               for j in 1:N]...)
            $([:(ary[i, $j] = $(var_names[j]) * $(eff_factors[j]))
               for j in 1:N]...)
        end
    end
    soa_run = quote
        unsafe_copy!(soary, ary)
        @inbounds @simd for i in 1:nele
            factor2i = sofactor2[i]
            $([:($(eff_factors[j]) = factor2i * $(splat_factors[j]);
                 $(var_names[j]) = soary[i, $j];
                 soary[i, $j] = $(var_names[j]) * $(eff_factors[j]))
               for j in 1:N]...)
        end
        unsafe_copy!(ary, soary)
    end
    quote
        $init_ex
        println("AoS")
        $rand_ex
        @time for n in 1:nloop
            $aos_run
        end
        println("SoA")
        $rand_ex
        @time for n in 1:nloop
            $soa_run
        end
    end
end

get_factors(n) = (map(x->exp(im * x), linspace(0f0, 10f0, n))...)

# test_scale(64, get_factors(5))
# test_scale(64, get_factors(10))
test_scale(64, get_factors(24))
println()
# test_scale(128, get_factors(5))
# test_scale(128, get_factors(10))
test_scale(128, get_factors(24))
println()
# test_scale(256, get_factors(5))
# test_scale(256, get_factors(10))
test_scale(256, get_factors(24))
println()
# test_scale(512, get_factors(5))
# test_scale(512, get_factors(10))
test_scale(512, get_factors(24))
println()
# test_scale(1024, get_factors(5))
# test_scale(1024, get_factors(10))
test_scale(1024, get_factors(24))

# @code_llvm test_scale(128, get_factors(2), 100000)
