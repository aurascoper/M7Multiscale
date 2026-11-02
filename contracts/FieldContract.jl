# M7Multiscale/contracts/FieldContract.jl
# Field descriptor contract with producer attribution, units, SI exponents, validity masks, and uncertainties.

using .M7VisualizationBridge:
    VALID_M7_PRODUCERS, VALID_M7_ASSOCIATIONS,
    M7FieldDescriptor, M7FieldDisplayProvenance,
    audit_m7_field_display_provenance
