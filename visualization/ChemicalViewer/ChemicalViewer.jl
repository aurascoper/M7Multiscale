# M7Multiscale/visualization/ChemicalViewer/ChemicalViewer.jl
# Workspace 2: Geochemical & Speciation World (PHREEQC pH, pe, ionic strength, species, mineral SI).

function extract_chemical_viewport_summary(frame::M7SpatialFrame)
    ph_fields = [f for f in frame.fields if f.producer == :PHREEQC]
    return (
        workspace = :chemical_world,
        producer = :PHREEQC,
        cell_count = frame.mesh.cell_count,
        field_names = [f.name for f in ph_fields],
        has_ph = any(f -> f.name == "pH", ph_fields),
        has_pe = any(f -> f.name == "pe", ph_fields)
    )
end
