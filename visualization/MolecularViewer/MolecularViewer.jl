# M7Multiscale/visualization/MolecularViewer/MolecularViewer.jl
# Workspace 4: Molecular & Macromolecular World (Sobolev SMACRO1 active-site + coarse-grained beads).

function extract_molecular_viewport_summary(frame::M7SpatialFrame)
    n_anchors = length(frame.molecular_anchors)
    total_full = sum((size(a.local_coords, 1) for a in frame.molecular_anchors); init=0)
    total_lod = sum((size(a.coarse_grain.coords, 1) for a in frame.molecular_anchors); init=0)
    return (
        workspace = :molecular_world,
        producer = :Sobolev,
        anchor_count = n_anchors,
        total_full_beads = total_full,
        total_lod_beads = total_lod,
        sobolev_relaxation_step = frame.clocks.sobolev_relaxation_step
    )
end
