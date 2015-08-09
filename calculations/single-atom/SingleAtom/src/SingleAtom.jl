#!/usr/bin/julia -f

module SingleAtom

include("utils.jl")
using .Utils
export StructOfArrays, SoCArray, SoCVector, SoCMatrix
export sum2average
export Vec3D

include("optical.jl")
using .Optical
export Drive

include("atomic.jl")
using .Atomic
export TransitionType, Trans_σ⁺, Trans_σ⁻, Trans_π
export Transition
export AtomBuilder, add_state!, add_transition!
export InternStates

include("system.jl")
using .System
export AbstractPotential, HarmonicPotential, ZeroPotential
export get_potential, get_kinetic
export SystemBuilder, add_state!, add_transition!, add_potential!, add_drive!
export MotionSystem

include("propagator.jl")
using .Propagate
export SystemPropagator, AbstractMeasure, DummyMeasure, MonteCarloMeasure
export AbstractSetup, StaticSetup, propagate, Propagate
export SnapshotType, SnapshotX, SnapshotK
export DecayType, DecayNone, DecayLeft, DecayRight, DecayMiddle

include("measure.jl")
using .Measure
export WaveFuncMeasure, EnergyMeasure
export WaveFuncMonteCarloMeasure, EnergyMonteCarloMeasure

include("plot.jl")
using .Plotting
export plot_measure

end
