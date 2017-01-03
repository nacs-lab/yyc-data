#!/usr/bin/julia

module Format

function fp_digit(io::IO, a::AbstractFloat, d::Int)
    isfinite(a) || throw(ArgumentError("Invalid floating point number $a"))
    isneg = false
    if a < 0
        a = -a
        isneg = true
    end
    ten = oftype(a, 10)
    scaled = a * ten^d
    if !isfinite(scaled) || eps(scaled) > 1
        throw(ArgumentError("$a does not have enough precision for $d digits"))
    end
    scaled_int = round(Integer, scaled)
    isneg && write(io, '-')
    if d == 0
        show(io, scaled_int)
        return
    end

    ten_int_n = oftype(scaled_int, 10)^d
    int_part = scaled_int ÷ ten_int_n
    float_part = scaled_int % ten_int_n
    show(io, int_part)
    d == 0 && return
    write(io, '.')
    if float_part == 0
        for i in 1:d
            write(io, '0')
        end
    else
        float_digits = floor(Integer, log10(float_part)) + 1
        @assert float_digits <= d
        for i in 1:(d - float_digits)
            write(io, '0')
        end
        show(io, float_part)
    end
    return
end

end
