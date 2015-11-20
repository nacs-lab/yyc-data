#!/usr/bin/julia -f

module TestDrive

using Base.Test

info("Test Drive")

import TwoLevels: AbstractDrive, getδ, getΩ, getϕ₀
import TwoLevels: PhaseTracker, update, getϕ
import TwoLevels: ConstDrive, LinearRampDrive, RampToDrive

info("  Test ConstDrive")
let
    drive = ConstDrive(1f0, 2f0)
    @test getδ(drive, 0f0, 1f0, 10f0) === 1f0
    @test getΩ(drive, 0f0, 1f0, 10f0) === 2f0
    @test getϕ₀(drive) === 0f0

    tracker = PhaseTracker{Float32}()
    @test update(tracker, drive, 0.1f0, 0.1f0, 1f0, 10f0) === nothing
    @test getδ(tracker, drive, 0.1f0, 1f0, 10f0) === 1f0
    @test getΩ(tracker, drive, 0.1f0, 1f0, 10f0) === 2f0
    @test getϕ(tracker) === 1f0 * 0.1f0
    @test update(tracker, drive, 0.1f0, 0.2f0, 1f0, 10f0) === nothing
    @test getδ(tracker, drive, 0.2f0, 1f0, 10f0) === 1f0
    @test getΩ(tracker, drive, 0.2f0, 1f0, 10f0) === 2f0
    @test getϕ(tracker) === 1f0 * 0.1f0 * 2
end

info("  Test LinearRampDrive")
let
    drive = LinearRampDrive(1f0, 2f0, -3f0, -4f0)
    @test getδ(drive, 1f0, 2f0, 10f0)::Float32 ≈ 1.5f0
    @test getΩ(drive, 1f0, 2f0, 10f0)::Float32 ≈ -3.5f0
    @test getϕ₀(drive) === 0f0

    tracker = PhaseTracker{Float32}()
    @test update(tracker, drive, 0.1f0, 0.1f0, 1f0, 10f0) === nothing
    @test getδ(tracker, drive, 0.1f0, 1f0, 10f0)::Float32 ≈ 1.1f0
    @test getΩ(tracker, drive, 0.1f0, 1f0, 10f0)::Float32 ≈ -3.1f0
    @test getϕ(tracker)::Float32 ≈ 0.105f0
    @test update(tracker, drive, 0.1f0, 0.2f0, 1f0, 10f0) === nothing
    @test getδ(tracker, drive, 0.2f0, 1f0, 10f0)::Float32 ≈ 1.2f0
    @test getΩ(tracker, drive, 0.2f0, 1f0, 10f0)::Float32 ≈ -3.2f0
    @test getϕ(tracker)::Float32 ≈ 0.22f0
end

info("  Test RampToDrive")
let
    drive = RampToDrive(-2f0, 4f0)
    @test getδ(drive, 1f0, 2f0, -1f0)::Float32 ≈ -1.5f0
    @test getΩ(drive, 1f0, 2f0, 3f0)::Float32 ≈ 3.5f0
    @test getϕ₀(drive) === 0f0

    tracker = PhaseTracker{Float32}()
    @test update(tracker, drive, 0.1f0, 0.1f0, 1f0, -1f0) === nothing
    @test getδ(tracker, drive, 0.1f0, 1f0, -1f0)::Float32 ≈ -1.1f0
    @test getΩ(tracker, drive, 0.1f0, 1f0, 3f0)::Float32 ≈ 3.1f0
    @test getϕ(tracker)::Float32 ≈ -0.105f0
    @test update(tracker, drive, 0.1f0, 0.2f0, 1f0, -1f0) === nothing
    @test getδ(tracker, drive, 0.2f0, 1f0, -1f0)::Float32 ≈ -1.2f0
    @test getΩ(tracker, drive, 0.2f0, 1f0, 3f0)::Float32 ≈ 3.2f0
    @test getϕ(tracker)::Float32 ≈ -0.22f0
end

end
