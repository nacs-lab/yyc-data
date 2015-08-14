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
    @test_approx_eq abs((k1, Trans_σ⁺) * amp_σ⁺) 1.0

    @test (k1, Trans_π) * amp_σ⁻ == 0
    @test (k1, Trans_π) * amp_σ⁺ == 0
    @test_approx_eq abs((k1, Trans_π) * amp_π) 1.0

    @test (k1, Trans_σ⁻) * amp_σ⁺ == 0
    @test (k1, Trans_σ⁻) * amp_π == 0
    @test_approx_eq abs((k1, Trans_σ⁻) * amp_σ⁻) 1.0
end

let
    builder = AtomBuilder{Float32}()
    add_state!(builder, :(G, 0, 0), :G, 0)
    add_state!(builder, :(E, 1, -1), :E, -1)
    add_state!(builder, :(E, 1, 0), :E, 0)
    add_state!(builder, :(E, 1, 1), :E, 1)

    add_transition!(builder, :(G, 0, 0), :(E, 1, -1),
                    Transition{Trans_σ⁻}(0.5f0, 10f0, 1f0))
    add_transition!(builder, :(G, 0, 0), :(E, 1, 0),
                    Transition{Trans_π}(0.5f0, 10f0, 1f0))
    add_transition!(builder, :(G, 0, 0), :(E, 1, 1),
                    Transition{Trans_σ⁺}(0.5f0, 10f0, 1f0))

    atom = InternStates(builder)
    @test isbits(atom)
end

end
