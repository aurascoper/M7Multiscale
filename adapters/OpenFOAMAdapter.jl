# M7Multiscale/adapters/OpenFOAMAdapter.jl
# OpenFOAM spatial mesh, boundary patch, velocity/pressure/tracer field, and CFL stability adapter.

using .FlowCytometryObservationModel:
    OpenFOAMBoundaryPatch, OpenFOAMMeshExchange, OpenFOAMTransportSnap,
    serialize_openfoam_ascii_field, parse_openfoam_ascii_field,
    map_openfoam_field_conservative, audit_openfoam_transport_stability,
    build_geometric_overlap_matrix, evaluate_mesh_mapping_reconstruction_error,
    summarize_openfoam_mesh_to_snap
