# M7Multiscale/visualization/BiofilmViewer/BiofilmViewer.jl
# Workspace 3: Biofilm, Cellular & Network World (Julia biomass, EPS/PGA/Curli, Droop Q, BARE graph,
# discrete low-occupancy extinction, and typed EPS interaction graph G_EPS(t) = (V, E_pp, E_pb, E_mineral)).

using .M7VisualizationBridge:
    extract_m7_extended_conservation_and_topology_summary

function extract_biofilm_viewport_summary(frame::M7SpatialFrame)
    jl_fields = [f for f in frame.fields if f.producer == :Julia]
    bg = frame.biological_graph
    return (
        workspace = :biological_world,
        producer = :Julia,
        field_names = [f.name for f in jl_fields],
        has_graph = (bg !== nothing),
        beta_0 = bg === nothing ? 0 : bg.beta_0,
        beta_1 = bg === nothing ? 0 : bg.beta_1,
        tarjan_bridge_count = bg === nothing ? 0 : bg.tarjan_bridge_count,
        topological_regime = bg === nothing ? :none : bg.topological_regime
    )
end

