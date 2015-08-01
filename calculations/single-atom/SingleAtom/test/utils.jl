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

# sum2average
# count == 1 should make the uncertainty undefined
@test isequal(sum2average(1.0, 1.0, 1), (1.0, NaN))
# rounding error in Σ(x²) shouldn't throw error
@test sum2average(2.0, 2.0 - 0.01, 2) == (1.0, 0.0)

@test_approx_eq [sum2average(6.0, 14.0, 3)...] [2.0, sqrt(1 / 3)]

# Vec3D

let
    e_x = Vec3D(1.0, 0.0, 0.0)
    e_y = Vec3D(0.0, 1.0, 0.0)
    e_z = Vec3D(0.0, 0.0, 1.0)

    # == and isequal
    @test e_x == e_x
    @test e_y == e_y
    @test e_z == e_z
    @test e_x === e_x
    @test e_y === e_y
    @test e_z === e_z
    @test isequal(e_x, e_x)
    @test isequal(e_y, e_y)
    @test isequal(e_z, e_z)

    @test e_x != e_y
    @test e_y != e_z
    @test e_x != e_z

    # ==, isqual and NaN, ±0.0
    @test Vec3D(NaN, 0.0, 0.0) != Vec3D(NaN, 0.0, 0.0)
    @test !(Vec3D(NaN, 0.0, 0.0) == Vec3D(NaN, 0.0, 0.0))
    @test isequal(Vec3D(NaN, 0.0, 0.0), Vec3D(NaN, 0.0, 0.0))
    @test Vec3D(1.0, -0.0, 0.0) == Vec3D(1.0, 0.0, -0.0)
    @test Vec3D(1.0, -0.0, 0.0) !== Vec3D(1.0, 0.0, -0.0)
end

end
