# M7Multiscale/coupling/BiologicalFeedback.jl
# Coupled biofilm EPS/PGA/Curli swelling, effective diffusivity D_eff(x,t), Droop quota Q,
# discrete low-occupancy absorbing extinction, dynamic EPS cross-linking (N_pp, N_pb, N_mineral),
# V. cholerae cell-surface remodeling, and typed EPS interaction graph G_EPS(t) = (V, E_pp, E_pb, E_mineral).

using .MicrofluidicClassificationEngine:
    BiofilmMatrixComposition, ThinFilmSwellingMechanicsSpec,
    evaluate_biofilm_swelling_mechanics, HeterogeneousDiffusivitySpec,
    compute_effective_diffusivity, calibrate_and_audit_effective_diffusivity,
    solve_1d_reactive_transport_biofilm, audit_droop_mass_conservation,
    run_biofilm_ablation_matrix,
    DiscretePopulationTransitionContract, HybridPopulationCellState,
    step_hybrid_population_extinction!, step_multicell_hybrid_dispersal!,
    DynamicEPSCrosslinkState, step_dynamic_eps_crosslinks!,
    simulate_cell_surface_remodeling!,
    EPSGraphVertex, EPSGraphEdge, TypedEPSInteractionGraph,
    build_typed_eps_interaction_graph

using .M7VisualizationBridge:
    M7BiologicalGraphNode, M7BiologicalGraphEdge, M7BiologicalGraphLayer,
    compute_biological_graph_topology

