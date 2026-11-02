# M7Multiscale/adapters/FlowCytometryAdapter.jl
# Julia <-> Flow Cytometry compartment observation-to-model alignment adapter (Milestone M7.4).

using .FlowCytometryObservationModel:
    CellStatePrediction, GateCalibration, FCMSampleSummary,
    MESFCalibrationCurve, BeadSizeCalibration,
    expected_cell_events, predicted_gate_probabilities, absolute_cells_per_ml,
    joint_loglikelihood, recover_latent_state_mle

using .M7VisualizationBridge:
    align_m7_frame_to_flow_cytometry
