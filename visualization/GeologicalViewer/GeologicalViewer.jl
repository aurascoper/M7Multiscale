# M7Multiscale/visualization/GeologicalViewer/GeologicalViewer.jl
# Workspace 1: Geological & Hydrodynamic World (OpenFOAM velocity, pressure, temperature, porosity).

function extract_geological_viewport_summary(frame::M7SpatialFrame)
    of_fields = [f for f in frame.fields if f.producer == :OpenFOAM]
    return (
        workspace = :geological_world,
        producer = :OpenFOAM,
        cell_count = frame.mesh.cell_count,
        topology = frame.mesh.topology,
        geometry_sha256 = frame.mesh.geometry_sha256,
        field_names = [f.name for f in of_fields],
        mean_porosity = sum(frame.mesh.porosity_fraction) / max(length(frame.mesh.porosity_fraction), 1)
    )
end
