#!/usr/bin/julia -f

module TestOptical

using SingleAtom
using SingleAtom.Optical
using Base.Test

let
    amp = Vec3D{Complex64}(1, 2, 3im)

    # Deterministic phase tracker
    drive1 = Drive{amp}(0.5, 1, 0, Inf)
    @test isa(drive1, Drive{amp,Float32})

    tracker1 = PhaseTracker(drive1)
    @test isa(tracker1, PhaseTracker{Float32})

    for i in 1:4
        init_phase!(tracker1)
        @test tracker1.phase == 0

        @test update_phase!(tracker1, 0.1f0) == (-0.1f0, exp(-0.1f0 * im))
        @test update_phase!(tracker1, 0.1f0) == (-0.2f0, exp(-0.2f0 * im))
    end

    # Random phase tracker
    drive2 = Drive{amp}(0.5, 1, NaN, 100)
    @test isa(drive2, Drive{amp,Float32})

    tracker2 = PhaseTracker(drive2)
    @test isa(tracker2, PhaseTracker{Float32})

    phases = [init_phase!(tracker2).phase for i in 1:100]
    max_phase = maximum(phases)
    min_phase = minimum(phases)
    @test 0 <= min_phase < max_phase < 2π

    for i in 1:4
        init_phase!(tracker2)

        has_rand::Bool = false
        for j in 1:10
            prev_phase = tracker2.phase
            @test tracker2.exp_t == exp(im * prev_phase)
            cur_phase, exp_t = update_phase!(tracker2, 0.1f0)
            if cur_phase != prev_phase - 0.1f0
                has_rand = true
            end
        end

        @test has_rand
    end
end

end
