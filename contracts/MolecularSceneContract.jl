# M7Multiscale/contracts/MolecularSceneContract.jl
# Sobolev SMACRO1 binary frame, far-field coarse-graining, and OpenUSD PointInstancer contract.

using .M7VisualizationBridge:
    SMACRO1_MAGIC, SMACRO1_VERSION,
    M7AssetPrototype, M7CoarseGrainResult, M7MolecularAnchor,
    encode_smacro1_frame, decode_smacro1_frame,
    coarse_grain_far_field, create_m7_molecular_anchor,
    write_m7_usda_point_instancer
