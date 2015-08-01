#!/usr/bin/julia -f

module TestUtils

using Base.Test

using SingleAtom
using SingleAtom.Utils

# Test Complex Structure of Arrays
for T in (Float32, Float64)
    aos_comp_1 = Complex{T}[(1:2048) * (1 + im);]
    aos_comp_2 = reshape(copy(aos_comp_1), (2, 1024))
    soa_comp_1 = convert(StructOfArrays, aos_comp_1)
    soa_comp_2 = convert(StructOfArrays, aos_comp_2)

    # Test alias
    @test isa(soa_comp_1, SoCVector{T})
    @test isa(soa_comp_2, SoCMatrix{T})

    @test aos_comp_1 == soa_comp_1
    @test aos_comp_2 == soa_comp_2

    aos_comp_1 = 10 ./ aos_comp_1
    aos_comp_2 = 10 ./ aos_comp_2

    @test aos_comp_1 != soa_comp_1
    @test aos_comp_2 != soa_comp_2

    # Test copy from AoS to SoA
    Base.unsafe_copy!(soa_comp_1, aos_comp_1)
    Base.unsafe_copy!(soa_comp_2, aos_comp_2)

    @test aos_comp_1 == soa_comp_1
    @test aos_comp_2 == soa_comp_2

    aos_comp_1 = 10 ./ aos_comp_1
    aos_comp_2 = 10 ./ aos_comp_2

    @test aos_comp_1 != soa_comp_1
    @test aos_comp_2 != soa_comp_2

    # Test copy from SoA to AoS
    Base.unsafe_copy!(aos_comp_1, soa_comp_1)
    Base.unsafe_copy!(aos_comp_2, soa_comp_2)

    @test aos_comp_1 == soa_comp_1
    @test aos_comp_2 == soa_comp_2
end

# Trig Cache
for T in (Float32, Float64)
    rang = T(1):T(1024)
    trig_cache = TrigCache(rang)
    @test isa(trig_cache, TrigCache{T})

    for i in 1:length(rang)
        @test trig_cache.sins[i] == sin(rang[i])
        @test trig_cache.coss[i] == cos(rang[i])
    end
end

# meta_expr

@test macroexpand(:(@meta_expr inline)) == Expr(:meta, :inline)
@test macroexpand(:(@meta_expr noinline)) == Expr(:meta, :noinline)

end
