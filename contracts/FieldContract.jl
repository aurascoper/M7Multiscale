# M7Multiscale/contracts/FieldContract.jl
# Field descriptor and frozen m7-field-bundle v1.0.0 contract with producer attribution, units,
# explicit basis, SI exponents, validity masks, uncertainties, and negative-fixture rejection gates.

using .M7VisualizationBridge:
    VALID_M7_PRODUCERS, VALID_M7_ASSOCIATIONS,
    VALID_M7_FIELD_BASES, VALID_M7_GRAPH_TYPES,
    M7_FIELD_BUNDLE_SCHEMA_VERSION,
    M7FieldDescriptor, M7FieldDisplayProvenance,
    M7FieldBundleV1,
    audit_m7_field_display_provenance,
    export_m7_field_bundle_v1,
    validate_m7_field_bundle_v1,
    verify_oklo_pool_mobile_cross_language_parity,
    evaluate_layer_down_solver_clearinghouse_ab

