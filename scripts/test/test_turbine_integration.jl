using Pkg; Pkg.activate("/home/daniel/cycle")
using GasCycle
const G = GasCycle

# NPSS turbine map (TurbinePRmap): independents (NpMap, PRmap); outputs WpMap and effMap (PR is an
# input, not an output). So we read Wp and eff along the 100% speed line and check the design point.
tables = G.read_npss_map("/home/daniel/bru_maps/bru_maps_project/maps/turbine/turbine_argon.map")
println("Parsed: ", collect(keys(tables)))

wc, ef = tables["TB_Wp"], tables["TB_eff"]
Np = wc.speeds[1]                                  # alpha index 1 (2-level table)
println("Np_axis = ", round.(Np, digits=3))
i = findfirst(x -> isapprox(x, 1.0; atol=1e-6), Np)
isnothing(i) && error("100% speed line not found")

# linear interpolation along the (sorted) PR axis with end clamping
function lininterp(x, y, xq)
    xq <= x[1]   && return y[1]
    xq >= x[end] && return y[end]
    k = searchsortedlast(x, xq)
    y[k] + (xq - x[k]) / (x[k+1] - x[k]) * (y[k+1] - y[k])
end

PR_des = 1.69                                       # design total-to-static PR (Fig 8 marker)
Wp_des = lininterp(wc.coords[1][i], wc.values[1][i], PR_des)
eta_des = lininterp(ef.coords[1][i], ef.values[1][i], PR_des)
println("Query at design (NpMap=1.0, PRmap=", PR_des, "): WpMap=", round(Wp_des, digits=3),
        "  eff=", round(eta_des, digits=3))
println("   expected ~ WpMap 1.0 (design flow), eff 0.913")
