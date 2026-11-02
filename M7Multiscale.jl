module M7Multiscale

include(joinpath(@__DIR__, "..", "flow-cytometry-model-bridge", "src", "M7VisualizationBridge.jl"))

using .M7VisualizationBridge
using .FlowCytometryObservationModel
using .MicrofluidicClassificationEngine

include(joinpath(@__DIR__, "contracts", "SpatialFrame.jl"))
include(joinpath(@__DIR__, "contracts", "TimeFrame.jl"))
include(joinpath(@__DIR__, "contracts", "FieldContract.jl"))
include(joinpath(@__DIR__, "contracts", "MolecularSceneContract.jl"))
include(joinpath(@__DIR__, "contracts", "ProvenanceContract.jl"))

include(joinpath(@__DIR__, "adapters", "OpenFOAMAdapter.jl"))
include(joinpath(@__DIR__, "adapters", "PHREEQCAdapter.jl"))
include(joinpath(@__DIR__, "adapters", "FlowCytometryAdapter.jl"))
include(joinpath(@__DIR__, "adapters", "OkloScenarioAdapter.jl"))

include(joinpath(@__DIR__, "coupling", "ConservativeMeshMapping.jl"))
include(joinpath(@__DIR__, "coupling", "ChemicalStateMapping.jl"))
include(joinpath(@__DIR__, "coupling", "BiologicalFeedback.jl"))
include(joinpath(@__DIR__, "coupling", "MultirateClock.jl"))

include(joinpath(@__DIR__, "visualization", "GeologicalViewer", "GeologicalViewer.jl"))
include(joinpath(@__DIR__, "visualization", "ChemicalViewer", "ChemicalViewer.jl"))
include(joinpath(@__DIR__, "visualization", "BiofilmViewer", "BiofilmViewer.jl"))
include(joinpath(@__DIR__, "visualization", "MolecularViewer", "MolecularViewer.jl"))
include(joinpath(@__DIR__, "visualization", "Timeline4D", "Timeline4D.jl"))

include(joinpath(@__DIR__, "validation", "MassConservation.jl"))
include(joinpath(@__DIR__, "validation", "ChargeConservation.jl"))
include(joinpath(@__DIR__, "validation", "SpatialAlignment.jl"))
include(joinpath(@__DIR__, "validation", "TemporalAlignment.jl"))
include(joinpath(@__DIR__, "validation", "EvidenceClassification.jl"))

export M7VisualizationBridge,
       extract_geological_viewport_summary,
       extract_chemical_viewport_summary,
       extract_biofilm_viewport_summary,
       extract_molecular_viewport_summary,
       extract_timeline4d_summary,
       extract_m7_extended_conservation_and_topology_summary,
       VALID_ACHIEVEMENT_LEVELS,
       VALID_EVIDENCE_CERTIFICATES,
       ACHIEVEMENT_STRATIFICATION_CONTRACT,
       MENKE_2026_BENCHMARK_REFERENCE,
       MenkeConservativeCellInventory,
       ExtendedBiologicalConservativeInventory,
       DiscretePopulationTransitionContract,
       HybridPopulationCellState,
       DynamicEPSCrosslinkState,
       EPSGraphVertex,
       EPSGraphEdge,
       TypedEPSInteractionGraph,
       NaCaMontmorilloniteExperimentSpec,
       audit_epistemic_achievement_stratification,
       evaluate_menke_conservative_inventory,
       evaluate_extended_biological_inventory,
       audit_four_compartment_mass_transfer,
       build_extended_inventory_from_hybrid_state,
       step_hybrid_population_extinction!,
       step_multicell_hybrid_dispersal!,
       step_dynamic_eps_crosslinks!,
       simulate_cell_surface_remodeling!,
       build_typed_eps_interaction_graph,
       evaluate_eps_topology_hypotheses,
       evaluate_na_vs_ca_montmorillonite_protocol,
       VALID_M7_EXPLORER_PANELS,
       VALID_M7_THREE_SYSTEMS,
       M7ThreeSystemExplorerSpec,
       export_m7_frame_xdmf_csv,
       reload_m7_frame_from_xdmf_csv,
       reload_m7_metadata_from_openusd,
       build_m7_three_system_explorer_frames,
       build_m7_explorer_payload,
       verify_m7_complete_roundtrip,
       evaluate_m7_acceptance_ladder,
       audit_csga_dsrb_motif_and_synonymous_group,
       M6_M7_CONCEPTUAL_BRIDGES_VS_FORMAL_ANALOGIES,
       MYCELIAL_ATTRACTOR_STATE_SPACE_DIMENSIONS,
       INTRAHYPHAL_TRANSPORT_PHENOTYPES,
       audit_m6_m7_conceptual_bridges_vs_formal_analogies,
       evaluate_mycelial_attractor_basin_remodeling,
       evaluate_chemotropic_homing_and_mutant_ablation,
       compute_m6_m7_persistent_homology_and_transport_metrics,
       CATALYTIC_COUPLE_MOTIF_JOINS,
       audit_catalytic_couple_motifs_and_cantor_joins,
       audit_pair_of_pants_cobordism_and_cantor_branch_contrast,
       evaluate_openfoam_regime_and_boundary_handoff,
       audit_conserved_quantity_ledger_across_shuffle,
       audit_glmm_openfoam_conjecture_and_empowerment,
       verify_m7_topological_paraview_timeseries_suite,
       verify_m7_five_layer_reference_dataset,
       VALID_M7_WORLD_EVIDENCE_TIERS,
       VALID_M7_SCENARIO_CLASSIFICATIONS,
       M7WorldScenario,
       validate_m7_world_scenario,
       build_canonical_m7_world_scenarios,
       simulate_m7_planetary_feedback_loop,
       serialize_m7_world_scenario_json,
       evaluate_m7_world_scenario_acceptance_suite,
       FUN_NEUROLOGY_COMPARTMENTALIZED_NETWORK_SPEC,
       FOUR_REGIME_SURFACE_SCARS_AND_ARC_CAPSID_BRIDGE,
       ASR_PHONETIC_RECOVERY_DICTIONARY_V2,
       build_fun_neurology_compartmentalized_network_graph,
       evaluate_surface_scars_equifinality_and_viral_arc_bridge,
       M7GeometryTopologyGate,
       validate_m7_geometry_import_topology,
       SPICE_LOCK,
       evaluate_dual_path_spice_ephemeris_contract,
       audit_m7_v2_non_recursive_poller_and_fun_neurology_graph

end # module M7Multiscale



