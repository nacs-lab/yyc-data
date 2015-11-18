#!/usr/bin/julia -f

module TwoLevels

include("drive.jl")
import .Drive: AbstractDrive, getδ, getΩ, getϕ₀
import .Drive: DriveTracker, update, getϕ
import .Drive: ConstDrive, LinearRampDrive, RampToDrive

include("measure.jl")
import .Measure: AbstractMeasure, measure_snapshot
import .Measure: DummyMeasure, MeasureWrapper

end
