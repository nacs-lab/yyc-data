#!/usr/bin/julia -f

module Atomic

using ..Utils
import ..Optical
import Base: *

export TransitionType, Trans_σ⁺, Trans_σ⁻, Trans_π

@enum TransitionType Trans_σ⁺ Trans_σ⁻ Trans_π

"""
The type of an optical dipole transition. σ⁺, σ⁻, or π
end JuliaLang/julia#12419
"""
TransitionType

"""
Overlap of a drive with a dipole transition
"""
function *(ax_trans::Tuple{Vec3D,TransitionType}, amp::Vec3D)
    # A = a * π + b * σ⁺ + c * σ⁻
    ax = ax_trans[1]
    trans = ax_trans[2]
    ax_scale = 1 / abs(ax)
    norm_ax = ax * ax_scale
    if trans == Trans_π
        # ax ⋅ π = 1
        # ax ⋅ σ⁺ = 0
        # ax ⋅ σ⁻ = 0
        return abs(norm_ax * amp)
    end
    # TODO figure out the sign of σ⁺ and σ⁻
    # ax × π = 0
    # ax × σ⁺ = i σ⁺
    # ax × σ⁻ = -i σ⁻
    trans_diff = norm_ax × amp # i (σ⁺ - σ⁻)
    # ax × (ax × π) = 0
    # ax × (ax × σ⁺) = -σ⁺
    # ax × (ax × σ⁻) = -σ⁻
    trans_amp = norm_ax × trans_diff # -(σ⁺ + σ⁻)
    σ_sign = ifelse(trans == Trans_σ⁺, -1im, 1im)
    return abs(trans_diff + σ_sign * trans_amp) / 2
end

*(amp::Vec3D, ax_trans::Tuple{Vec3D,TransitionType}) = ax_trans * amp

###
# Internal states of the atom

immutable Transition{T,Pol} # Pol::TransitionType
    k::T # wave vector
    Γ::T # Line width
    α::T # Dipole moment
end

end
