# M7Multiscale/contracts/SpatialFrame.jl
# Canonical M7SpatialFrame v1.0.0 spatial coordinate and mesh contract.

using .M7VisualizationBridge:
    M7_SCHEMA_NAME, M7_SCHEMA_VERSION,
    VALID_M7_MESH_TOPOLOGIES, VALID_M7_CELL_ORDERINGS, VALID_M7_LOCAL_UNITS,
    M7CoordinateFrameSpec, M7NestedTransform, M7MeshDescriptor, M7SpatialFrame,
    compose_nested_transform, inverse_nested_transform,
    transform_local_to_world, transform_world_to_local,
    transform_instance_matrices_to_world, audit_nested_transform,
    serialize_m7_spatial_frame_json, validate_m7_spatial_frame_schema
