#!/usr/bin/julia -f

module TwoLevels

include("drive.jl")
import .Drives: AbstractDrive, getδ, getΩ, getϕ₀
import .Drives: PhaseTracker, update, getϕ
import .Drives: ConstDrive, LinearRampDrive, RampToDrive, SinsDrive

include("measure.jl")
import .Measures: AbstractMeasure
import .Measures: MeasureWrapper, MeasureList
import .Measures: DummyMeasure, FullMeasure, SingleMeasure

# sequence
include("sequence.jl")
import .Sequences: Sequence, propagate

# builder
include("builder.jl")
import .Builders: SeqBuilder

end
