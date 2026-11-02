# M7Multiscale/validation/MassConservation.jl
# OpenFOAM conservative field mapping, Droop nutrient quota mass balance, Menke (2026) benchmark
# (M = theta * C + a_surf * Gamma with epsilon != theta in general), and proposed 4-compartment
# biological extension (M_total = theta * C + a_surf * Gamma + X * q + E).

using .M7VisualizationBridge:
    M7ConservationAuditReport,
    audit_m7_spatiotemporal_conservation,
    MENKE_2026_BENCHMARK_REFERENCE,
    MenkeConservativeCellInventory,
    ExtendedBiologicalConservativeInventory,
    evaluate_menke_conservative_inventory,
    evaluate_extended_biological_inventory,
    audit_four_compartment_mass_transfer,
    build_extended_inventory_from_hybrid_state

