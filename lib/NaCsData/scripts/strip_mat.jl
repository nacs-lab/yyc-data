#!/usr/bin/julia --startup-file=no

using MAT

iname = ARGS[1]
oname = ARGS[2]
if isdir(oname)
    oname = joinpath(oname, "$(splitext(basename(iname))[1]).mat")
end

const mf = matopen(iname)
const pl = read(mf, "ParamList")
const sa = read(mf, "Analysis")["SingleAtomLogical"]
const param_name = read(mf, "Scan")["ParamName"]
close(mf)

const out = matopen(oname, "w")
write(out, "ParamList", pl)
write(out, "SingleAtomLogical", sa)
write(out, "ParamName", param_name)
close(out)
