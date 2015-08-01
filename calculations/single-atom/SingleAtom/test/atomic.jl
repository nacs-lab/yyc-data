#!/usr/bin/julia -f

module TestAtomic

using Base.Test
using SingleAtom
using SingleAtom.Atomic

# Polarization and transition types
let
    k1 = Vec3D(1.0, 0.0, 0.0)
    amp_σ⁺ = Vec3D{Complex128}(0.0, 1.0, -1.0im) / √2
    amp_σ⁻ = Vec3D{Complex128}(0.0, 1.0, 1.0im) / √2
    amp_π = Vec3D{Complex128}(1.0, 0.0, 0.0)

    @test (k1, Trans_σ⁺) * amp_π == 0
    @test (k1, Trans_σ⁺) * amp_σ⁻ == 0
    @test_approx_eq (k1, Trans_σ⁺) * amp_σ⁺ 1.0

    @test (k1, Trans_π) * amp_σ⁻ == 0
    @test (k1, Trans_π) * amp_σ⁺ == 0
    @test_approx_eq (k1, Trans_π) * amp_π 1.0

    @test (k1, Trans_σ⁻) * amp_σ⁺ == 0
    @test (k1, Trans_σ⁻) * amp_π == 0
    @test_approx_eq (k1, Trans_σ⁻) * amp_σ⁻ 1.0
end

end
