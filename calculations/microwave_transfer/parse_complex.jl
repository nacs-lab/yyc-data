#!/usr/bin/julia

function parse_complex(s::String)
    Base.@gc_preserve s begin
        ps = pointer(s)
        pout = Ref{Ptr{UInt8}}(0)
        vre = ccall(:strtod, Cdouble, (Ptr{UInt8}, Ptr{Ptr{UInt8}}), ps, pout)
        if ps == pout[]
            throw(ArgumentError("Cannot parse the real part."))
        end
        pm = unsafe_load(pout[])
        if pm == 0
            return Complex(vre)
        elseif pm != UInt8('+') && pm != UInt8('-')
            throw(ArgumentError("Cannot parse the real part."))
        end
        ps = pout[]
        vim = ccall(:strtod, Cdouble, (Ptr{UInt8}, Ptr{Ptr{UInt8}}), ps, pout)
        pi = unsafe_load(pout[])
        if ps == pout[] || pi != UInt8('i')
            throw(ArgumentError("Cannot parse the imaginary part."))
        end
        return vre + vim * im
    end
end
