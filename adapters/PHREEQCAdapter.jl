# M7Multiscale/adapters/PHREEQCAdapter.jl
# PHREEQC geochemical speciation, charge/element balance, activity coefficient, and mineral SI adapter.

using .FlowCytometryObservationModel:
    PHREEQCMineralPhase, PHREEQCSpeciesAdapter, PHREEQCSpeciationSnap,
    audit_phreeqc_speciation_balance, verify_phreeqc_database_compatibility,
    convert_molality_to_molarity, run_menke_2026_retardation_benchmark,
    serialize_phreeqc_selected_output, parse_phreeqc_selected_output,
    summarize_phreeqc_adapter_to_snap
