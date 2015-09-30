#!/usr/bin/julia -f

function read_beamage_txt(fd::IO)
    len = 0
    data_lines = Vector{Int}[]
    for line in eachline(fd)
        isempty(line) && continue
        if !(';' in line)
            len == 0 && continue
            throw(ParseError("Non-data line following data lines"))
        end
        data = Int[parse(Int, field) for field in split(strip(line), ';',
                                                        keep=false)]
        line_len = length(data)
        len == 0 || line_len == len || throw(ParseError("Data length mismatch"))
        len = line_len
        push!(data_lines, data)
    end
    res = Matrix{Int}(len, length(data_lines))
    @inbounds for i in 1:length(data_lines)
        data = data_lines[i]
        @simd for j in 1:len
            res[j, i] = data[j]
        end
    end
    res
end

function read_beamage_txt(fname)
    open(read_beamage_txt, fname)
end
