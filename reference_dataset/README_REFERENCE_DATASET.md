# M7 Known-Good Five-Layer Reference Dataset & Planetary WorldScenario Extension

M7 is an integration, comparison, and provenance-preserving 3D/4D process-engineering workbench that sits above independent numerical solvers and structural tools. It does not replace OpenFOAM, PHREEQC, ParaView, PyMOL, or UCSF ChimeraX.

This directory provides a self-contained, known-good reference dataset where every layer can be opened and verified independently before loading into the M7 browser workbench (`visualization/m7_multiscale_explorer.html` / `workbench/m7_multiscale_explorer.html`).

## 1. Layer 1 (CFD + Topological Time-Series + OpenFOAM/GLMM): OpenFOAM Fixture and ParaView 6.2.0 Suite

Directory: `openfoam_case/`
- **Epistemic & Mesh Resolution Distinction**:
  - **Primary OpenFOAM Fixture Mesh (`4-Cell Hexahedral OpenFOAM Fixture (20 vertices)`)**: The OpenFOAM reference dataset in `openfoam_case/` is a reproducible visualization and coupling fixture (`evidence_tier = EARTH_REFERENCE_FIXTURE / SYNTHETIC_FIXTURE`, generated deterministically by `m7_paraview_interactive_state.py`), not evidence of a completed external `porousSimpleFoam` binary solve. Its `constant/polyMesh` defines 20 vertices, 21 faces, and 4 hexahedral cells (`1.0e-9 m^3` each).
  - **Overlap Validation Grid (`32^3 Conservative Voxel-Overlap Validation Grid`)**: The `32^3 = 32,768` voxel grid belongs to the separate conservative geometric overlap mapping and mass-conservation validation path (`ConservativeMeshMapping.jl`), not the cell count of `openfoam_case/`.
- `0/`, `300/`, `600/`, `900/`, `1200/`: Five physical time-step directories (`t = 0, 300, 600, 900, 1200 s`, pinned to Sobolev relaxation steps `0, 100, 200, 300, 400`) on the 4-cell porous clay-biofilm channel (`mesh_00035_clay_pore`, right-handed Cartesian meters). Each time directory contains 12 OpenFOAM v2.0 ASCII field dictionaries:
  - Transport & Geochemistry: `U`, `p`, `porosity`, `tracer`, `pH` (activity-averaged $\mathrm{pH}_{\mathrm{avg}}(t) = -\log_{10}(\mathrm{mean}(10^{-\mathrm{pH}_i(t)}))$ progressing `7.2134 -> 7.1245 -> 7.0317 -> 7.0317 -> 7.0317`).
  - Topological Homology: `beta_0`, `beta_1` (tracking mycelial/EPS anastomosis loop closure from $\beta_0=2,\beta_1=0$ unclosed tree at $t=0\text{ s}$ -> $\beta_1=1$ primary `so`/`Prm1` loop at $t=300\text{ s}$ -> $\beta_1=2$ pair-of-pants $\Sigma_{0,3}$ bifurcation at $t=600\text{ s}$ -> $\beta_1=3$ hollow-shell $\Psi$ excursion at $t=900\text{ s}$ -> $\beta_1=2$ Pareto trap-space anastomosis attractor at $t=1200\text{ s}$).
  - Coupled OpenFOAM/GLMM & Empowerment: `glmm_observed_y`, `glmm_fixed_effect_beta0_plus_beta1_x`, `glmm_random_intercept_b_i`, `glmm_residual_eps_ij` ($y_{ij}(t) = (\beta_0 + \beta_1 x_{ij}(t)) + b_i + \varepsilon_{ij}(t)$ with $\beta_0=1.20, \beta_1=0.85, b = [-0.06, -0.02, +0.03, +0.05], \mathrm{ICC} = 0.98537$), and `empowerment_bits` ($\mathfrak{E}(s_t) = \max_{p(a_t|s_t)} I(A_t; S_{t+1}|s_t)$).
- `constant/polyMesh/{points,faces,owner,neighbour,boundary}` and `constant/transportProperties`: Full hexahedral mesh topology (20 vertices, 21 faces, 4 cells) and Newtonian/Darcy transport properties.
- `system/{controlDict,fvSchemes,fvSolution}`: Solver controls (`porousSimpleFoam`, `startTime = 0`, `endTime = 1200`, `writeInterval = 300`).
- `case.foam`: Empty reader entry-point marker for ParaView (`SkipZeroTime = 0` loads all 5 time steps `[0.0, 300.0, 600.0, 900.0, 1200.0]`).
- `m7_topological_timeseries.pvd` + `timeseries/m7_channel_t{0000..1200}s.vtu`: ParaView PVD/VTU 3D hexahedral time-series collection with cell and point fields.
- `m7_topology_tubes.pvd` + `timeseries/m7_tubes_t{0000..1200}s.vtp`: ParaView PVD/VTP 3D mycelial/EPS anastomosis tube network time series with criterial separation (`criterial_status_code = +1.0` green verified active, `0.0` gold exploratory, `-1.0` red criterial pruned).
- `m7_topological_timeseries.xdmf`: XDMF v3.0 `<Grid GridType="Collection" CollectionType="Temporal">` across `t = 0, 300, 600, 900, 1200 s`.
- `m7_topological_glmm_timeseries.pvsm`, `remap_pvsm_paths.py`, and `m7_paraview_interactive_state.py`: Portable ParaView 6.2.0 interactive state file (with zero hardcoded `/Users/` paths; all reader paths are relative to `openfoam_case/`), CLI portability remapper/verifier (`python3 remap_pvsm_paths.py --check-only`), and `pvpython` generator/verifier with Slice, Contour, StreamTracer, Glyph, and Tube filters plus `1920x1080` comparison renders in `renders/m7_paraview_t{0000..1200}s.png` and `renders/m7_paraview_topological_glmm_comparison_1920x1080.png`.
- `m7_openfoam_reference.vtk` and `m7_openfoam_reference.xdmf` + `m7_openfoam_reference.csv`: Standalone single-frame baseline VTK and XDMF v3.0 exports.

### How to open independently in ParaView 6.2.0:
```bash
cd reference_dataset/openfoam_case
# 1. Start with the relative PVD time series or OpenFOAM case marker:
open -a ParaView m7_topological_timeseries.pvd
open -a ParaView case.foam
# 2. Verify or remap .pvsm paths if relocating to a custom working directory:
python3 remap_pvsm_paths.py --check-only
open -a ParaView m7_topological_glmm_timeseries.pvsm
# 3. Or re-run the headless + interactive pvpython verifier:
/Applications/ParaView-6.2.0.app/Contents/bin/pvpython m7_paraview_interactive_state.py
```

## 2. Layer 2 (Geochemistry): PHREEQC Reference Speciation and Cell Mapping

Directory: `phreeqc/`
- `reference_speciation.pqi`: Standalone PHREEQC input script defining `SOLUTION 1` (Upstream, `pH = 7.35, pe = 4.10, T = 25.0 degC`) and `SOLUTION 2` (Downstream, `pH = 6.85, pe = 2.40, T = 25.0 degC`), `EQUILIBRIUM_PHASES` (`Calcite`, `Ce_Oxalate_10H2O`), and `SURFACE` (`Mont_OH`). Hydrothermal temperatures are strictly quarantined as explicit scenario inputs rather than solver defaults.
- `reference_selected_output.tsv` and `reference_speciation.pqo`: Selected-output table and run log demonstrating zero electrical charge imbalance (`0.00%`), exact elemental conservation (`Ca`, `Fe`, `C`, `Cl`), and activity-based pH averaging (`pH_activity_avg = -log10(mean(10.^(-pH))) = 7.03169896` vs arithmetic `mean(pH) = 7.10000000`).
- `cell_mapping.csv`: Explicit cell-to-OpenFOAM region mapping table.

## 3. Layer 3 & Layer 4 (Coupling, Conserved-Quantity Ledger, `M7WorldScenario` & M7 Data Contract)

Verified programmatically via `verify_m7_five_layer_reference_dataset()` and `evaluate_m7_world_scenario_acceptance_suite()` in `M7VisualizationBridge.jl`:
- **Planetary `M7WorldScenario` & Closed Feedback Loop**: Parameterizes `environment`, `fluid_system`, `geological_substrate`, `energy_sources`, `redox_chemistry`, `biological_system`, and `evidence_tier` across three canonical worlds:
  1. `earth_oklo_natural_reactor_reference` (`EARTH_REFERENCE_FIXTURE` / `:literature_derived`)
  2. `enceladus_tiger_stripe_ocean_floor` (`PLANETARY_HYPOTHETICAL` / `:planetary_hypothetical`)
  3. `mars_noachian_radiolytic_subsurface` (`PLANETARY_HYPOTHETICAL` / `:planetary_hypothetical`)
  and couples `Biofilm growth -> EPS production -> Kozeny-Carman pore restructuring -> Darcy transport redistribution -> redox/mineral changes -> new biological niches` gated by strict mass and free-energy budget closure (`elemental_mass_conserved`, `energy_budget_closed`).
- **Spatial consistency**: Right-handed Cartesian convention, `32^3` conservative voxel overlap mapping preserving integrated quantities, partial-overlap / out-of-domain validity mask checks, and reversible `M7NestedTransform` (`angstrom`/`nm`/`um` <-> `m`) transforms.
- **Conserved-quantity ledger across permutations (`audit_conserved_quantity_ledger_across_shuffle`)**: Verifies that an external ledger tracks volume-integrated mass, charge, elemental totals, and cobordism boundary counts across cell/row permutations.
- **Temporal & provenance consistency**: Keeps `physical_time_s`, `biological_time_s`, and `sobolev_relaxation_step` strictly distinct and verifies SHA-256 digests across all files and frames.

## 4. Layer 5 (Molecular Biology): PyMOL and UCSF ChimeraX Independent Inspection

Directory: `molecular/`
- `dsrB_2V4J_catalytic_cluster.pdb`: **Composite active-site excerpt** (`EXPDTA X-RAY DIFFRACTION`) containing selected active-site chain fragments from three separate source crystal structures: `DsrB` (`PDB 2V4J`, `UniProt P45575`, `Cys261+Cys264+Cys268` in `CGNCYTMC` with `[Fe4S4]` and siroheme), `OahA` (`PDB 3M0J`, `Arg21+Ser26` in `MGFRIMIFSFA` with `Mg2+` and oxalate), and `OxDC` (`PDB 1J58`, `His95+His97+Glu101` with `Mn2+` and oxalate). It is a coordinate inspection excerpt, not a single crystallographic co-assembly of all three proteins; inspect the full original structures directly from the PDB archive via `fetch 2v4j 3m0j 1j58` (PyMOL) or `open 2v4j 3m0j 1j58` (ChimeraX).
- `csgA_curli_repeat_model.pdb`: Badged explicitly as `[THEORETICAL_MODEL_FRAGMENT]` (`EXPDTA THEORETICAL MODEL`, `UniProt P28307`), containing only an 18-residue R1 repeat fragment (`Ser43..Gln60` `SELNIYQYGNGNAALQTQ` with 12-aa core `Gln49..Gln60` `QYGNGNAALQTQ`), strictly quarantined from experimental structures and contrasted with full-length experimental cryo-EM curli fibril assemblies `PDB 8ENQ` and `PDB 8ENR` (`fetch 8enq 8enr` / `open 8enq 8enr`).
- `inspect_active_sites.pml`: Ready-to-run PyMOL script with composite-excerpt disclosures and `fetch 2v4j 3m0j 1j58` / `fetch 8enq 8enr` commands.
- `inspect_active_sites.cxc`: Ready-to-run UCSF ChimeraX script with composite-excerpt disclosures and `open 2v4j 3m0j 1j58` / `open 8enq 8enr` commands.

### How to open independently in PyMOL and ChimeraX:
```bash
cd reference_dataset/molecular
pymol inspect_active_sites.pml
chimerax inspect_active_sites.cxc
```

## 5. Support Micromamba Environments & Browser Workbench

Directory: `environments/`
Note: `m7-analysis.yml` (Python, NumPy, SciPy, pandas, meshio, PyVista, VTK) and `m7-molecular.yml` (Biopython, Biotite, MDTraj, PDBFixer) are Python analysis and molecular-data support environments. They do not install OpenFOAM, PHREEQC, ParaView, PyMOL, or ChimeraX, which remain separately installed external solver and visualization applications.
```bash
micromamba create -f reference_dataset/environments/m7-analysis.yml -y
micromamba create -f reference_dataset/environments/m7-molecular.yml -y
```

### Opening the M7 Browser Workbench:
```bash
open visualization/m7_multiscale_explorer.html
# Or from workbench/:
open workbench/m7_multiscale_explorer.html
```
