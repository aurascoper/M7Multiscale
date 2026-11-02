# M7Multiscale/validation/SpatialAlignment.jl
# Nested SO(3) chirality, unit scale, coordinate round-trip, and cross-panel spatial picking validation.

using .M7VisualizationBridge:
    audit_nested_transform,
    pick_m7_spatial_location,
    inspect_m7_microscope_scale
