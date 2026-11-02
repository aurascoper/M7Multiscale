# M7Multiscale/coupling/ConservativeMeshMapping.jl
# Geometry-aware conservative 3D bounding-box overlap mapping across OpenFOAM, PHREEQC, and Julia.

using .FlowCytometryObservationModel:
    build_geometric_overlap_matrix,
    map_openfoam_field_conservative,
    evaluate_mesh_mapping_reconstruction_error

using .M7VisualizationBridge:
    map_multiphysics_fields_to_m7_mesh
