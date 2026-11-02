using Test

include(joinpath(@__DIR__, "..", "M7Multiscale.jl"))
using .M7Multiscale
using .M7Multiscale.M7VisualizationBridge

@testset "M7Multiscale Repository Architecture & Contract Suite" begin
    refs = build_m7_reference_case_frames()

    geo = extract_geological_viewport_summary(refs.frame_2)
    @test geo.workspace == :geological_world
    @test geo.producer == :OpenFOAM
    @test geo.cell_count == 4

    chem = extract_chemical_viewport_summary(refs.frame_2)
    @test chem.workspace == :chemical_world
    @test chem.producer == :PHREEQC
    @test chem.has_ph && chem.has_pe

    bio = extract_biofilm_viewport_summary(refs.frame_2)
    @test bio.workspace == :biological_world
    @test bio.producer == :Julia
    @test bio.has_graph
    @test bio.beta_1 == 2

    mol = extract_molecular_viewport_summary(refs.frame_2)
    @test mol.workspace == :molecular_world
    @test mol.producer == :Sobolev
    @test mol.anchor_count == 1
    @test mol.total_full_beads == 8
    @test mol.total_lod_beads == 5

    t4d = extract_timeline4d_summary([refs.frame_1, refs.frame_2])
    @test t4d.workspace == :spatiotemporal_world
    @test t4d.frame_count == 2
    @test t4d.conservation_audit.all_checks_passed

    oklo = import_oklo_frozen_fixture_to_m7_frame()
    @test oklo.track_audit.separation_valid
    @test validate_m7_spatial_frame_schema(oklo.frame)
    @test validate_m7_spatial_frame_schema(serialize_m7_spatial_frame_json(oklo.frame))

    cap = evaluate_capillary_meniscus_habitat(CapillaryMeniscusHabitatSpec())
    @test cap.cell_fully_immersed
    @test cap.within_microct_5_to_150um_band

    @test_throws ArgumentError evaluate_hydrothermal_vent_plume_profile(180.0, 0.0)
    plume = evaluate_hydrothermal_vent_plume_profile(
        180.0, 0.0;
        magma_temp_degC = 1000.0,
        orifice_temp_degC = 350.0,
        plume_core_temp_degC = 2.05,
        ambient_ocean_temp_degC = 2.00
    )
    @test plume.zone == :buoyant_plume

    # Step 1-5 verification through M7Multiscale
    @test MENKE_2026_BENCHMARK_REFERENCE.closed_batch_max_rel_inventory_change == 8.9e-10
    ext_sum = extract_m7_extended_conservation_and_topology_summary(
        induce_extinction = true,
        cell_surface_remodeling_active = true
    )
    @test ext_sum.discrete_population_state.is_extinct
    @test ext_sum.discrete_population_state.discrete_occupancy_count == 0
    @test ext_sum.discrete_population_state.mass_conservation_audit.is_conserved
    @test ext_sum.eps_crosslink_state.remodeling_report.vps_matrix_intact
    @test ext_sum.eps_crosslink_state.remodeling_report.cells_released_from_intact_mesh
    @test ext_sum.hypothesis_certificates.H1 == :witness
    @test ext_sum.hypothesis_certificates.H2 == :witness
    @test ext_sum.hypothesis_certificates.H3 == :witness

    # Complete round-trip, 4-panel explorer payload, 3-system selector, V1-V6 ladder, and motif audit
    rt = verify_m7_complete_roundtrip()
    @test rt.roundtrip_passed
    @test rt.field_integrals_preserved
    @test rt.units_preserved
    @test rt.coordinates_preserved
    @test rt.three_clocks_preserved
    @test rt.provenance_sha256_preserved

    mktempdir() do tmpd
        exp_b = export_m7_frame_xdmf_csv(tmpd, refs.frame_2; prefix = "m7_standalone")
        rel_b = reload_m7_frame_from_xdmf_csv(exp_b.xdmf_path)
        @test rel_b.sha256_verified
        @test rel_b.frame.frame_sha256 == refs.frame_2.frame_sha256
        @test rel_b.frame.biological_graph !== nothing
        @test length(rel_b.frame.molecular_anchors) == 1
    end

    three_sys = build_m7_three_system_explorer_frames()
    @test !three_sys.integrated_into_existing_openfoam_phreeqc_solvers
    for sym in VALID_M7_THREE_SYSTEMS
        s = getfield(three_sys, sym)
        @test s.spec.epistemic_system_status == :modeled_system
        @test s.spec.validation_realization_status == :not_experimentally_validated_realization
        @test !s.spec.integrated_into_existing_openfoam_phreeqc_solvers
        f_names = Set(f.name for f in s.frame.fields)
        @test "surface_OCa_plus" in f_names
        @test "surface_O_minus" in f_names
    end

    payload = build_m7_explorer_payload(refs.frame_2)
    @test keys(payload.panels) == VALID_M7_EXPLORER_PANELS
    @test payload.panels.panel_b_geochemistry.has_surface_complexation

    ladder = evaluate_m7_acceptance_ladder()
    @test ladder.all_gates_passed
    @test ladder.gates.V1_geometry_mesh_spatial_frame_consistency
    @test ladder.gates.V2_conservative_openfoam_field_mapping
    @test ladder.gates.V3_phreeqc_species_units_chemical_state_integrity
    @test ladder.gates.V4_independent_clock_preservation
    @test ladder.gates.V5_difference_frame_and_visual_provenance
    @test ladder.gates.V6_scene_export_reload_and_reproducibility

    motif = audit_csga_dsrb_motif_and_synonymous_group()
    @test motif.all_reconciled
    @test motif.csga.core_length_aa == 12
    @test motif.csga.repeat_length_aa == 18
    @test motif.dsrb.n_X_internal_wildcards == 5
    @test motif.synonymous_group.log2_group_order == 38

    # Known-bad structure fail-loud inverse test:
    # 1. Corrupted SMACRO1 magic header or truncated payload must throw ArgumentError
    @test_throws ArgumentError decode_smacro1_frame(UInt8[0x00, 0x01, 0x02])
    bad_magic = Vector{UInt8}(undef, 32)
    fill!(bad_magic, 0x42)
    @test_throws ArgumentError decode_smacro1_frame(bad_magic)
    # 2. Tampered in-memory M7SpatialFrame field values or malformed JSON SHA-256 must throw ArgumentError
    tampered_frame = deepcopy(refs.frame_2)
    tampered_frame.fields[1].values[1, 1] += 1.0
    @test_throws ArgumentError validate_m7_spatial_frame_schema(tampered_frame)
    valid_json = serialize_m7_spatial_frame_json(refs.frame_2)
    tampered_json = replace(valid_json, "sha256:" => "sha256:corrupted_")
    @test_throws ArgumentError validate_m7_spatial_frame_schema(tampered_json)
    # 3. Missing hydrothermal boundary conditions must throw ArgumentError
    @test_throws ArgumentError evaluate_hydrothermal_vent_plume_profile(100.0, 0.0)
    # 4. Linear interpolation of Sobolev relaxation steps as physical time must be blocked
    f_sob_bad = M7SpatialFrame(
        simulation_id = refs.frame_1.simulation_id,
        frame_index = refs.frame_1.frame_index,
        clocks = M7MultirateClock(
            physical_time_s = refs.frame_1.clocks.physical_time_s,
            biological_time_s = refs.frame_1.clocks.biological_time_s,
            sobolev_relaxation_step = refs.frame_1.clocks.sobolev_relaxation_step + 50,
            snapshot_sync_id = "sync-bad-sob"
        ),
        mesh = refs.frame_1.mesh,
        fields = refs.frame_1.fields
    )
    @test_throws ArgumentError interpolate_m7_spatial_frames(refs.frame_1, f_sob_bad, refs.frame_1.clocks.physical_time_s)

    # M6-7 Persistent Homology (H_0, H_1 Betti curves), cycle rank, mean transport path length, and Conceptual Bridges vs Formal Analogies audit
    cb_audit = audit_m6_m7_conceptual_bridges_vs_formal_analogies()
    @test cb_audit.all_bridges_strictly_separated
    @test cb_audit.m6_internal_project_label_disclosed
    @test cb_audit.dimension_count == 4
    @test cb_audit.non_mechanistic_guardrails.binary_nucleotide_flags_not_folding_mechanics
    @test cb_audit.non_mechanistic_guardrails.synonymous_codons_not_independent_binary_flips
    @test cb_audit.non_mechanistic_guardrails.cycle_rank_not_3d_tertiary_structure
    @test cb_audit.non_mechanistic_guardrails.tree_automorphisms_not_biophysical_control

    ph_metrics = compute_m6_m7_persistent_homology_and_transport_metrics()
    @test ph_metrics.node_count == 5
    @test first(ph_metrics.H0_betti_curve) == 5
    @test last(ph_metrics.H0_betti_curve) == 1
    @test maximum(ph_metrics.H1_betti_curve) >= 2
    @test ph_metrics.active_network.beta_0 == 1
    @test ph_metrics.active_network.beta_1 == 2
    @test isfinite(ph_metrics.active_network.mean_euclidean_path_length)
    @test ph_metrics.ablation_regimes.disable_anastomosis.beta_1 == 0
    @test ph_metrics.ablation_regimes.disable_anastomosis.mean_euclidean_path_length >=
          ph_metrics.ablation_regimes.full_coupled.mean_euclidean_path_length

    # 4D Mycelial State-Space Attractor Remodeling, Chemotropic Homing, Mutant Ablations, and 5 Transport Phenotypes
    @test rt.mycelial_attractor_basin_converged
    att = ph_metrics.mycelial_attractor_remodeling
    @test att.trap_space_converged
    @test keys(att.state_space_dimensions) == (
        :topological_density,
        :edge_conductivity,
        :spatial_partitioning,
        :information_routing_spectral_entropy
    )
    @test att.divergent_remodeling_paths.high_centrality_cord_thickening_and_anastomosis.mean_diameter_increase_ratio > 1.0
    @test att.divergent_remodeling_paths.high_centrality_cord_thickening_and_anastomosis.anastomosis_cycles_retained == 2
    @test att.divergent_remodeling_paths.low_transport_link_pruning.mean_diameter_attenuation_ratio < 1.0
    @test att.quantitative_state_space_metrics.topological_density.pareto_front_trapped
    @test att.quantitative_state_space_metrics.edge_conductivity.fixed_point_stabilized
    @test att.quantitative_state_space_metrics.spatial_partitioning.phase_field_converged
    @test att.quantitative_state_space_metrics.information_routing_spectral_entropy.avoids_chaotic_dissipation

    hom = ph_metrics.chemotropic_homing_and_mutants
    @test hom.all_homing_and_mutant_gates_verified
    @test hom.chemotropic_homing_gate.self_fusion.regulatory_loci == (:so, :Prm1)
    @test hom.chemotropic_homing_gate.self_fusion.anastomosis_pore_open
    @test !hom.chemotropic_homing_gate.self_fusion.cell_lysis_triggered
    @test hom.chemotropic_homing_gate.self_fusion.cord_thickening_active
    @test hom.chemotropic_homing_gate.non_self_conflict.regulatory_loci == (:het,)
    @test hom.chemotropic_homing_gate.non_self_conflict.cell_lysis_triggered
    @test !hom.chemotropic_homing_gate.non_self_conflict.anastomosis_pore_open
    @test isapprox(hom.mutant_ablation_states.delta_Prm1.translocation_reduction_fraction, 0.50; atol = 1.0e-12)
    @test hom.mutant_ablation_states.delta_so.heterogeneous_media_translocation_flux == 0.0
    @test 0.0 <= hom.reticulation_index_alpha <= 1.0

    phen = ph_metrics.intrahyphal_transport_phenotypes
    @test keys(phen) == (
        :turgor_mass_flow,
        :trunk_channeling,
        :pulsatile_exchange,
        :loop_redundancy,
        :autophagic_recycling
    )

    cc_audit = audit_catalytic_couple_motifs_and_cantor_joins()
    @test cc_audit.all_checks_passed
    @test cc_audit.catalytic_couple_count == 6
    @test cc_audit.all_catalytic_couples_verified
    @test !cc_audit.catalytic_couples.couple_1_radiolytic_h2_sulfate_reduction.shared_flat_primary_key
    @test cc_audit.catalytic_couples.couple_1_radiolytic_h2_sulfate_reduction.catalytic_turnover_join_valid
    @test cc_audit.separation_purity_specification.within_radionuclidic_spec_at_calibration
    @test cc_audit.separation_purity_specification.within_elemental_purity_spec
    @test cc_audit.analytical_fraction_separation.elemental_fraction_icp_ms.instrument == :ICP_MS
    @test cc_audit.analytical_fraction_separation.structural_mineral_fraction_xrd_xrf.instrument == :XRD_XRF
    @test cc_audit.analytical_fraction_separation.cellular_uptake_fraction_flow_cytometry.instrument == :Flow_Cytometry
    @test cc_audit.cobordism_homology_primitive.all_cobordism_invariants_verified
    @test cc_audit.cobordism_homology_primitive.pair_of_pants_primitive.beta_1 == 2
    @test cc_audit.cobordism_homology_primitive.higher_genus_cobordism_decomposition.beta_1 == 9
    @test cc_audit.cobordism_homology_primitive.higher_genus_cobordism_decomposition.pairs_of_pants_count == 8
    @test cc_audit.cobordism_homology_primitive.cantor_tree_vs_cobordism_contrast.unclosed_tree_beta_1 == 0

    cc_equal_hl = audit_catalytic_couple_motifs_and_cantor_joins(product_half_life_h = 60.0, contaminant_half_life_h = 60.0)
    @test cc_equal_hl.all_checks_passed
    @test cc_equal_hl.separation_purity_specification.decay_corrected_validity_window_h == Inf

    cc_eop_bad = audit_catalytic_couple_motifs_and_cantor_joins(contaminant_initial_activity_Bq = 5.0e3, max_radionuclidic_impurity_ratio = 1.0e-3)
    @test !cc_eop_bad.all_checks_passed
    @test cc_eop_bad.separation_purity_specification.decay_corrected_validity_window_h == 0.0

    cob_closed = audit_pair_of_pants_cobordism_and_cantor_branch_contrast(genus = 2, boundary_components = 0)
    @test cob_closed.all_cobordism_invariants_verified
    @test cob_closed.higher_genus_cobordism_decomposition.pairs_of_pants_count == 2

    # OpenFOAM regime classification and conserved-flux boundary-condition handoff (dec_038)
    of_reg = evaluate_openfoam_regime_and_boundary_handoff(wall_temperature_K = 310.15, bulk_temperature_K = 298.15)
    @test of_reg.all_checks_passed
    @test of_reg.boundary_condition_handoff_ledger.mass_flux_conserved
    @test of_reg.boundary_condition_handoff_ledger.wall_heat_flux_W_m2 > 0.0

    # External Conserved-Quantity Ledger across permutations / shuffles (dec_030, dec_034, dec_035, dec_036, dec_037)
    shuf_led = audit_conserved_quantity_ledger_across_shuffle(refs.frame_2; permutation = [4, 2, 1, 3])
    @test shuf_led.ledger_verified
    @test shuf_led.adjacency_destroyed_by_shuffle
    @test shuf_led.max_conserved_integral_rel_error <= 1.0e-12
    @test shuf_led.max_spatial_remap_recovery_abs_error <= 1.0e-12
    @test shuf_led.cobordism_boundary_count_preserved

    # 5-Layer Known-Good Reference Dataset & Workbench HTML verification
    wb_html_path = normpath(joinpath(@__DIR__, "..", "workbench", "m7_multiscale_explorer.html"))
    @test isfile(wb_html_path)
    ref5 = verify_m7_five_layer_reference_dataset()
    @test ref5.all_five_layers_verified
    @test ref5.environments_present
    @test ref5.layer_1_cfd_openfoam.passed
    @test ref5.layer_1_cfd_openfoam.voxel_32cubed_cell_count == 32768
    @test ref5.layer_2_geochemistry_phreeqc.passed
    @test ref5.layer_2_geochemistry_phreeqc.arithmetic_vs_activity_pH_bias > 0.05
    @test ref5.layer_3_coupling_and_conserved_ledger.passed
    @test ref5.layer_4_m7_data_contract_and_clocks.passed
    @test ref5.layer_5_molecular_biology_structures.passed
    @test ref5.layer_5_molecular_biology_structures.experimental_pdb_2V4J_3M0J_1J58.is_experimental_structure
    @test ref5.layer_5_molecular_biology_structures.theoretical_model_csgA_P28307.is_theoretical_model

    conj_emp = audit_glmm_openfoam_conjecture_and_empowerment()
    @test conj_emp.all_verified
    @test conj_emp.dec_040_conjecture_comparator.reputation_rejected_as_category_error
    @test conj_emp.dec_041_upstream_and_downstream_composition.downstream_glmm_valid
    @test isapprox(conj_emp.dec_042_empowerment_and_conserved_option_space.empowerment_cantor_T26_bits, 6.0; atol = 1.0e-12)
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.verified
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.disk_logs_cross_verified
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.turn_count == 21
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.first_turn_index == 27
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.last_turn_index == 47
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.total_health_turn_count == 48
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.latest_health_turn_count == 50
    @test conj_emp.post_dec046_creole_hoffman_eeg_telemetry.prosody_trace_line_count == 60
    @test length(conj_emp.post_dec046_creole_hoffman_eeg_telemetry.post_000950_silent_turns) == 2

    pv_suite = verify_m7_topological_paraview_timeseries_suite()
    @test pv_suite.passed
    @test pv_suite.timesteps_s == [0.0, 300.0, 600.0, 900.0, 1200.0]
    @test pv_suite.global_beta_0_trajectory == [2, 1, 1, 1, 1]
    @test pv_suite.global_beta_1_trajectory == [0, 1, 2, 3, 2]
    @test pv_suite.glmm_max_decomposition_abs_error <= 1.0e-12
    @test pv_suite.pvd_valid
    @test pv_suite.xdmf_temporal_valid
    @test pv_suite.pvsm_valid
    @test pv_suite.criterial_tubes_valid
    @test pv_suite.all_pngs_1920x1080
    @test pv_suite.zero_em_en_dashes

    # Perturbation fail-loud test on GLMM residual and 0/beta_0
    m7_ref_dir = normpath(joinpath(@__DIR__, "..", "reference_dataset"))
    mktempdir() do tmp_root
        tmp_ref = joinpath(tmp_root, "ref_copy")
        cp(m7_ref_dir, tmp_ref; force = true)
        eps_600_path = joinpath(tmp_ref, "openfoam_case", "600", "glmm_residual_eps_ij")
        orig_eps = read(eps_600_path, String)
        write(eps_600_path, replace(orig_eps, "-0.006" => "-0.005"))
        pv_perturbed_glmm = verify_m7_topological_paraview_timeseries_suite(dataset_dir = tmp_ref, atol = 1.0e-12)
        @test !pv_perturbed_glmm.passed
        @test pv_perturbed_glmm.glmm_max_decomposition_abs_error > 1.0e-4

        write(eps_600_path, orig_eps)
        b0_0_path = joinpath(tmp_ref, "openfoam_case", "0", "beta_0")
        orig_b0 = read(b0_0_path, String)
        write(b0_0_path, replace(orig_b0, "\n2\n2\n" => "\n1\n1\n"))
        pv_perturbed_b0 = verify_m7_topological_paraview_timeseries_suite(dataset_dir = tmp_ref, atol = 1.0e-12)
        @test !pv_perturbed_b0.beta_progression_valid
        @test !pv_perturbed_b0.passed
    end

    # Planetary WorldScenario Contract & Closed Feedback Loop Verification
    ws_suite = evaluate_m7_world_scenario_acceptance_suite()
    @test ws_suite.all_passed
    @test ws_suite.epistemic_tiers_distinct
    @test ws_suite.terrestrial_gravity_reuse_blocked
    @test ws_suite.fake_solver_promotion_blocked
    @test ws_suite.feedback_budgets_closed
    @test ws_suite.unclosed_budget_rejected
    @test ws_suite.json_contract_complete
    @test ws_suite.html_loader_complete
    @test ws_suite.pvsm_zero_users_paths

    # Verify Enceladus uses H2 + HCO3- -> CH4 while Earth Oklo and Mars use H2 + SO4-2 -> HS-
    @test ws_suite.simulations.enceladus.donor_key == "H2"
    @test ws_suite.simulations.enceladus.acceptor_key == "HCO3-"
    @test ws_suite.simulations.enceladus.product_key == "CH4"
    @test ws_suite.simulations.earth_oklo.donor_key == "H2"
    @test ws_suite.simulations.earth_oklo.acceptor_key == "SO4-2"
    @test ws_suite.simulations.earth_oklo.product_key == "HS-"
    @test ws_suite.simulations.mars_noachian.donor_key == "H2"
    @test ws_suite.simulations.mars_noachian.acceptor_key == "SO4-2"
    @test ws_suite.simulations.mars_noachian.product_key == "HS-"

    # Verify M7SpatialFrame schema & field display provenance on simulated frames
    @test validate_m7_spatial_frame_schema(ws_suite.simulations.enceladus.spatial_frame)
    @test validate_m7_spatial_frame_schema(ws_suite.simulations.earth_oklo.spatial_frame)
    @test audit_m7_field_display_provenance(ws_suite.simulations.enceladus.spatial_frame).all_fields_provenance_complete
    @test audit_m7_field_display_provenance(ws_suite.simulations.earth_oklo.spatial_frame).all_fields_provenance_complete

    # M7 v2: Generalized Compartmentalized "Fun Neurology" Network Graph (Turns L110-L114),
    # 4-Regime Surface-Scars & Viral Arc Capsid + Basilisk/Vizard/SPICE.jl Bridge (Turns L75, L79, L102-L109, L115),
    # and Turn L112 Non-Recursive Poller Verification Suite
    fn_graph = build_fun_neurology_compartmentalized_network_graph(time_s = 6.0, drive_amplitude = 1.0)
    @test fn_graph.n_vertices == 5
    @test fn_graph.n_edges == 7
    @test fn_graph.beta_0 == 1
    @test fn_graph.beta_1 == 3
    @test fn_graph.beta_1_strictly_positive
    @test fn_graph.tau_fast_s == 1.2
    @test fn_graph.tau_slow_s == 8.0
    @test fn_graph.consolidation_rule_immutable
    @test fn_graph.recursion_depth == 0
    @test fn_graph.argument_is_dance_not_war
    @test_throws ArgumentError build_fun_neurology_compartmentalized_network_graph(attempt_recursive_rule_rewrite = true)

    surf_arc = evaluate_surface_scars_equifinality_and_viral_arc_bridge()
    @test surf_arc.all_four_surface_regimes_present
    @test surf_arc.viral_arc_capsid.triangulation_number_T == 4
    @test surf_arc.viral_arc_capsid.capsid_diameter_nm == 32.0
    @test surf_arc.viral_arc_capsid.stokes_creeping_flow_valid
    @test surf_arc.equifinality_and_surface_scars.equifinal_forward_converged
    @test surf_arc.equifinality_and_surface_scars.unique_backward_inversion_forbidden
    @test isapprox(surf_arc.basilisk_vizard_spice_bridge.spice_irradiance_W_m2, 1361.0; atol = 1.0e-9)
    @test_throws ArgumentError evaluate_surface_scars_equifinality_and_viral_arc_bridge(heliocentric_distance_au = -1.0)

    v2_audit = audit_m7_v2_non_recursive_poller_and_fun_neurology_graph()
    @test v2_audit.all_passed
    @test v2_audit.dialectic_log_line_count >= 115
    @test v2_audit.ground_truth_turns_verified.L79_viral_arc_openfoam_pymol_chimerax
    @test v2_audit.ground_truth_turns_verified.L109_four_regime_surface_scars
    @test v2_audit.ground_truth_turns_verified.L111_fun_neurology_tau_1_2s_8_0s
    @test v2_audit.ground_truth_turns_verified.L112_non_recursive_consolidation_guard
    @test v2_audit.ground_truth_turns_verified.L114_argument_dance_beta_1_nonzero
    @test v2_audit.ground_truth_turns_verified.L115_vizard_spice_julia_kernels
    @test v2_audit.poller_v2_verified
    @test v2_audit.master_event_contract_gates_verified
    @test v2_audit.substrate_swap_L29_to_L32_verified
    @test v2_audit.recursive_rewrite_blocked
    @test v2_audit.doc_valid
    @test v2_audit.html_v2_valid

    # Verified Claims Memo (2026-11-01), Single Boundary Converter convert_spice_state_to_basilisk_si, and SMBR Bridge
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.memo_date == "2026-11-01"
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_2_spice_planet_state_msg_payload.outputs ==
          (:spiceTimeOutMsg, :planetStateOutMsgs, :scStateOutMsgs, :attRefStateOutMsgs, :transRefStateOutMsgs)
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_4_vizard_live_modes.handshake_tokens ==
          ("PING", "SIM_UPDATE", "OK", "REQUEST_INPUT", "VIZARD_INPUT", "REQUEST_IMAGE_")
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_4_vizard_live_modes.clock_synch.c_header_field_status == :Incomplete
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_5_traceable_seven_budget_criteria.c5_atp_yield_uncertainty_interval_mol_per_mol_co2 == (0.5, 1.5)
    @test VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_5_traceable_seven_budget_criteria.c7_planetprofile_and_europa_status.europa_nutmeg_module_status == :Incomplete
    @test length(VERIFIED_BASILISK_SPICE_M7_MEMO_2026_11_01_SPEC.section_9_downloads_images_2026_11_01) == 5

    dcm_test = [cos(0.1) sin(0.1) 0.0; -sin(0.1) cos(0.1) 0.0; 0.0 0.0 1.0]
    c_un = convert_spice_state_to_basilisk_si([2.279392e8, 0.0, 0.0], [0.0, 24.07, 0.0], dcm_test; j2000_current_s = 600.0, planet_name = "MARS")
    @test c_un.spice_jl_unit_verdict == :Incomplete
    @test !c_un.units_verified_si
    @test c_un.planet_frame == "IAU_MARS"
    @test isapprox(c_un.PositionVector[1], 2.279392e11; rtol = 1.0e-12)
    @test isapprox(c_un.VelocityVector[2], 24070.0; rtol = 1.0e-12)
    @test isapprox(c_un.J20002Pfix, dcm_test; atol = 1.0e-12)
    c_chk = convert_spice_state_to_basilisk_si([2.279392e8, 0.0, 0.0], [0.0, 24.07, 0.0], dcm_test; j2000_current_s = 600.0, planet_name = "MARS", expected_distance_m_bounds = (2.0e11, 2.6e11))
    @test c_chk.spice_jl_unit_verdict == :Accepted
    @test c_chk.units_verified_si
    @test_throws ArgumentError convert_spice_state_to_basilisk_si([1.0, 2.0], [0.0, 1.0, 0.0])

    memo_res = evaluate_verified_basilisk_spice_memo_and_smbr_bridge(temperature_K = 298.15, atp_per_ch4_mol = 1.0)
    @test memo_res.all_verified
    @test memo_res.memos_identical
    @test memo_res.memos_zero_dashes
    @test memo_res.converter_valid
    @test memo_res.downloads_images_verified
    @test memo_res.smbr_acceptance_passed
    @test !memo_res.below_280K_underconstrained
    @test memo_res.within_calibrated_288_371K
    memo_cold_res = evaluate_verified_basilisk_spice_memo_and_smbr_bridge(temperature_K = 273.15, atp_per_ch4_mol = 0.5)
    @test memo_cold_res.below_280K_underconstrained
    @test_throws ArgumentError evaluate_verified_basilisk_spice_memo_and_smbr_bridge(atp_yield_interval_mol_per_mol_co2 = (1.0, 1.0))
    @test_throws ArgumentError evaluate_verified_basilisk_spice_memo_and_smbr_bridge(atp_per_ch4_mol = 2.0)

    matlab_alara_aud = audit_m7_matlab_simulink_3d_eda_and_alara_hippocratism()
    @test matlab_alara_aud.all_passed
    @test matlab_alara_aud.fixtures_present
    @test matlab_alara_aud.zero_dashes_fixtures
    @test matlab_alara_aud.receipt_checks_passed
    @test matlab_alara_aud.mat_v7_and_v73_hdf5_present
    @test matlab_alara_aud.simulink_slx_present
    @test matlab_alara_aud.png_dimensions_valid
    @test matlab_alara_aud.png_dimensions.width >= 1920
    @test matlab_alara_aud.png_dimensions.height >= 1080
    @test matlab_alara_aud.plugin_and_workspace_skills_verified
    @test matlab_alara_aud.health_physics_alara.alara_satisfied
    @test matlab_alara_aud.health_physics_alara.primum_non_nocere_enforced
    @test matlab_alara_aud.payload_controller.alara_satisfied
    @test matlab_alara_aud.melanin_neurospora_alara.alara_biological_protection_verified
    @test matlab_alara_aud.epistemic_hippocratism_gate.verdict == :Accepted
    @test matlab_alara_aud.m7_schema_validation.schema_valid
    @test matlab_alara_aud.selfish_gene_comfort_asr.reticulated_asr_verified
    @test matlab_alara_aud.tarjan_bridge_dag.is_acyclic_tarjan_bridge_chain
    @test matlab_alara_aud.matlab_neuralcompose_multipool_seam.collar_and_multipool_valid
    @test matlab_alara_aud.matlab_neuralcompose_multipool_seam.unified_substrate_swap_L31_valid
    @test matlab_alara_aud.matlab_neuralcompose_multipool_seam.multipool_biofilm_seam.biomass_conserved
    @test matlab_alara_aud.openmc_cross_section_swap_L32.valid
    @test matlab_alara_aud.neurospora_4block_8pool_state.all_4_blocks_valid
    @test matlab_alara_aud.neurospora_4block_8pool_state.n_pool === Int64(8)
    @test length(matlab_alara_aud.neurospora_4block_8pool_state.C_O2_i) == 8
    @test length(matlab_alara_aud.neurospora_4block_8pool_state.conidial_cytometry.FSC_A) == 8
    @test matlab_alara_aud.neurospora_4block_8pool_state.hormesis_quarantine_verdict == :Incomplete
end






