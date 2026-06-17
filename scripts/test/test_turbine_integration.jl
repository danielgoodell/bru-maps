using Pkg; Pkg.activate("/home/daniel/cycle")
using GasCycle
const G = GasCycle
tables = G.read_npss_map("/home/daniel/bru_maps/bru_maps_project/maps/turbine/turbine_argon.map")
println("Parsed: ", collect(keys(tables)))
pm = G.to_performance_map(tables; flow="S_map.WcMap", pr="S_map.PRmap", eff="S_map.effMap")
println("Nc_axis = ", round.(pm.Nc_axis, digits=3))
PR, eta = G.query(pm, 1.0, 1.0)
println("Query (Nc=1.0, Wc=1.0): PR=", round(PR,digits=3), " eta=", round(eta,digits=3), "  (expect ~1.69 / 0.913)")
