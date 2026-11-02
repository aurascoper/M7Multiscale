# M7Multiscale/adapters/OkloScenarioAdapter.jl
# Oklo <-> M7 6-field hydrothermal voxel scenario and frozen-fixture importer (Milestone M7.6, Track A).

using .M7VisualizationBridge:
    VALID_OKLO_PULSE_MODES,
    OkloScenarioVoxelSpec,
    evaluate_oklo_conformance_voxel,
    import_oklo_frozen_fixture_to_m7_frame
