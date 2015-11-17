#!/usr/bin/julia -f

module TwoLevels

include("drive.jl")
import .Drive: AbstractDrive, getδ, getΩ, getϕ₀
import .Drive: DriveTracker, update, getϕ
import .Drive: ConstDrive, LinearRampDrive, RampToDrive

end
