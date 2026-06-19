using Pkg; Pkg.activate("/home/daniel/cycle")
using GasCycle
const G = GasCycle
mappath = "/home/daniel/bru_maps/bru_maps_project/maps/compressor/compressor_argon.map"
tables = G.read_npss_map(mappath)
println("Parsed tables: ", collect(keys(tables)))
pm = G.to_performance_map(tables; flow="TB_Wc", pr="TB_PR", eff="TB_eff")
println("Nc_axis = ", round.(pm.Nc_axis, digits=3))
println("Wc_axis range = ", round(minimum(pm.Wc_axis),digits=3), " .. ", round(maximum(pm.Wc_axis),digits=3))
PR, eta = G.query(pm, 1.0, 1.0)   # design: Nc=1, Wc=1 (normalized)
println("Query at design (Nc=1.0, Wc=1.0):  PR=", round(PR,digits=3), "  eta=", round(eta,digits=3))
println("   expected ~ PR 1.90, eta 0.795")
