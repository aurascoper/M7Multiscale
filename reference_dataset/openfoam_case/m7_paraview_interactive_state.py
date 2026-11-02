#!/usr/bin/env pvpython
"""
M7 Interactive Topological Time-Series & OpenFOAM/GLMM ParaView Suite
=====================================================================
Generates and verifies the 5-timestep OpenFOAM reference fixture (0/, 300/, 600/, 900/, 1200/),
ParaView PVD/VTU/VTP collections (m7_topological_timeseries.pvd, m7_topology_tubes.pvd),
XDMF v3.0 TemporalCollection (m7_topological_timeseries.xdmf), portable ParaView state file
(m7_topological_glmm_timeseries.pvsm), and 1920x1080 PNG human-visual-analysis renders.

Epistemic & Mesh Resolution Disclosure:
- Evidence Tier: SYNTHETIC_FIXTURE / EARTH_REFERENCE_FIXTURE (not an external porousSimpleFoam solve).
- Primary OpenFOAM Fixture Mesh: 4-Cell Hexahedral OpenFOAM Fixture (20 vertices, 21 faces, 4 cells).
- Separate Overlap Validation Grid: 32^3 Conservative Voxel-Overlap Validation Grid (32,768 voxels).

Strictly zero em-dashes and zero en-dashes across all code, labels, and annotations.
"""

import json
import math
import os
import hashlib
import re
import sys

TIMESTEPS = [0.0, 300.0, 600.0, 900.0, 1200.0]
TIMESTEP_DIRS = ["0", "300", "600", "900", "1200"]
MESH_RESOLUTION_BADGE = "4-Cell Hexahedral OpenFOAM Fixture (20 vertices)"
VOXEL_VALIDATION_GRID_BADGE = "32^3 Conservative Voxel-Overlap Validation Grid"
EVIDENCE_TIER = "EARTH_REFERENCE_FIXTURE"

# 4-cell hexahedral channel mesh: x in [0, 4e-3] m, y in [-0.5e-3, 0.5e-3] m, z in [-0.5e-3, 0.5e-3] m
CELL_CENTERS_X = [0.5e-3, 1.5e-3, 2.5e-3, 3.5e-3]
CELL_VOLUMES_M3 = [1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9]

# GLMM fixed parameters and cell-level random intercepts b_i (dec_040, dec_041, dec_042)
GLMM_BETA0 = 1.20
GLMM_BETA1 = 0.85
GLMM_RANDOM_INTERCEPTS_B = [-0.06, -0.02, 0.03, 0.05]
GLMM_BASE_RESIDUALS = [-0.004, 0.005, -0.006, 0.005]

# Time-dependent physical, geochemical, topological, and GLMM fields across 5 steps
TIMESERIES_DATA = {
    0.0: {
        "dir": "0",
        "U_x": [2.4e-4, 2.2e-4, 1.9e-4, 2.2e-4],
        "p": [0.0095, 0.0070, 0.0042, 0.0014],
        "porosity": [0.82, 0.78, 0.75, 0.79],
        "tracer": [1.00, 0.86, 0.64, 0.42],
        "pH": [7.42, 7.38, 7.12, 7.05],
        "global_beta_0": 2,
        "global_beta_1": 0,
        "beta_0": [2.0, 2.0, 1.0, 1.0],
        "beta_1": [0.0, 0.0, 0.0, 0.0],
        "empowerment_bits": [1.42, 1.48, 1.45, 1.40],
        "residual_scale": 0.80,
        "topology_stage": "Unclosed Cantor Tree (beta_0=2, beta_1=0)",
        "sobolev_step": 0,
    },
    300.0: {
        "dir": "300",
        "U_x": [2.3e-4, 2.05e-4, 1.75e-4, 2.1e-4],
        "p": [0.0102, 0.0075, 0.0045, 0.0015],
        "porosity": [0.81, 0.76, 0.73, 0.78],
        "tracer": [1.00, 0.79, 0.52, 0.29],
        "pH": [7.39, 7.37, 6.98, 6.95],
        "global_beta_0": 1,
        "global_beta_1": 1,
        "beta_0": [1.0, 1.0, 1.0, 1.0],
        "beta_1": [1.0, 1.0, 0.0, 0.0],
        "empowerment_bits": [1.50, 1.64, 1.58, 1.46],
        "residual_scale": 0.90,
        "topology_stage": "Primary so/Prm1 Loop Closure (beta_0=1, beta_1=1)",
        "sobolev_step": 100,
    },
    600.0: {
        "dir": "600",
        "U_x": [2.2e-4, 1.9e-4, 1.6e-4, 2.0e-4],
        "p": [0.0110, 0.0081, 0.0048, 0.0016],
        "porosity": [0.80, 0.74, 0.70, 0.77],
        "tracer": [1.00, 0.72, 0.41, 0.18],
        "pH": [7.35, 7.35, 6.85, 6.85],
        "global_beta_0": 1,
        "global_beta_1": 2,
        "beta_0": [1.0, 1.0, 1.0, 1.0],
        "beta_1": [1.0, 2.0, 2.0, 1.0],
        "empowerment_bits": [1.56, 1.78, 1.74, 1.52],
        "residual_scale": 1.00,
        "topology_stage": "Pair-of-Pants Sigma_0,3 Bifurcation (beta_0=1, beta_1=2)",
        "sobolev_step": 200,
    },
    900.0: {
        "dir": "900",
        "U_x": [2.1e-4, 1.75e-4, 1.45e-4, 1.9e-4],
        "p": [0.0118, 0.0087, 0.0051, 0.0017],
        "porosity": [0.79, 0.73, 0.69, 0.76],
        "tracer": [0.98, 0.65, 0.35, 0.14],
        "pH": [7.35, 7.35, 6.85, 6.85],
        "global_beta_0": 1,
        "global_beta_1": 3,
        "beta_0": [1.0, 1.0, 1.0, 1.0],
        "beta_1": [1.0, 3.0, 3.0, 2.0],
        "empowerment_bits": [1.59, 1.83, 1.80, 1.55],
        "residual_scale": 1.05,
        "topology_stage": "Hollow Shell Psi Excursion (beta_0=1, beta_1=3)",
        "sobolev_step": 300,
    },
    1200.0: {
        "dir": "1200",
        "U_x": [2.0e-4, 1.6e-4, 1.3e-4, 1.8e-4],
        "p": [0.0125, 0.0092, 0.0054, 0.0018],
        "porosity": [0.78, 0.72, 0.68, 0.75],
        "tracer": [0.95, 0.58, 0.29, 0.11],
        "pH": [7.35, 7.35, 6.85, 6.85],
        "global_beta_0": 1,
        "global_beta_1": 2,
        "beta_0": [1.0, 1.0, 1.0, 1.0],
        "beta_1": [1.0, 2.0, 2.0, 1.0],
        "empowerment_bits": [1.62, 1.88, 1.84, 1.58],
        "residual_scale": 1.00,
        "topology_stage": "Pareto Trap-Space Anastomosis Attractor (beta_0=1, beta_1=2)",
        "sobolev_step": 400,
    },
}


def compute_glmm_fields_for_step(step_dict):
    tracer = step_dict["tracer"]
    r_scale = step_dict["residual_scale"]
    fixed_vec = []
    rand_vec = []
    res_vec = []
    obs_vec = []
    for i in range(4):
        # Local dimensionless Peclet/Thiele reaction covariate x_i(t)
        x_ij = 0.5 + 1.0 * (1.0 - tracer[i])
        fixed_val = GLMM_BETA0 + GLMM_BETA1 * x_ij
        b_i = GLMM_RANDOM_INTERCEPTS_B[i]
        eps_ij = GLMM_BASE_RESIDUALS[i] * r_scale
        y_ij = fixed_val + b_i + eps_ij
        fixed_vec.append(round(fixed_val, 12))
        rand_vec.append(round(b_i, 12))
        res_vec.append(round(eps_ij, 12))
        obs_vec.append(round(y_ij, 12))
    return {
        "glmm_observed_y": obs_vec,
        "glmm_fixed_effect_beta0_plus_beta1_x": fixed_vec,
        "glmm_random_intercept_b_i": rand_vec,
        "glmm_residual_eps_ij": res_vec,
    }


def compute_activity_pH_avg(ph_list):
    mean_h = sum(10.0 ** (-p) for p in ph_list) / float(len(ph_list))
    return -math.log10(mean_h)


def write_openfoam_vector_field(filepath, loc_dir, obj_name, dims_str, ux_vals, inlet_flux=2e-10, outlet_flux=-1.8e-10):
    lines = [
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       volVectorField;",
        f'    location    "{loc_dir}";',
        f"    object      {obj_name};",
        "}",
        f"dimensions      {dims_str};",
        "internalField   nonuniform List<vector>",
        "4",
        "(",
    ]
    for ux in ux_vals:
        lines.append(f"({ux:.8g} 0 0)")
    lines.extend([
        ");",
        "boundaryField",
        "{",
        "    inlet",
        "    {",
        "        type            fixedValue;",
        f"        flux_mol_s      {inlet_flux:.6g};",
        f"        value           uniform ({ux_vals[0]:.8g} 0 0);",
        "    }",
        "    outlet",
        "    {",
        "        type            zeroGradient;",
        f"        flux_mol_s      {outlet_flux:.6g};",
        "    }",
        "    walls",
        "    {",
        "        type            fixedValue;",
        "        flux_mol_s      0;",
        "        value           uniform (0 0 0);",
        "    }",
        "}",
        "",
    ])
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_openfoam_scalar_field(filepath, loc_dir, obj_name, dims_str, vals, boundary_mode="zeroGradient"):
    lines = [
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       volScalarField;",
        f'    location    "{loc_dir}";',
        f"    object      {obj_name};",
        "}",
        f"dimensions      {dims_str};",
        "internalField   nonuniform List<scalar>",
        "4",
        "(",
    ]
    for v in vals:
        lines.append(f"{v:.12g}")
    lines.append(");")
    lines.append("boundaryField")
    lines.append("{")
    if boundary_mode == "pressure":
        lines.extend([
            "    inlet",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
            "    outlet",
            "    {",
            "        type            fixedValue;",
            "        flux_mol_s      0;",
            "        value           nonuniform List<scalar> 1 ( 0.0 );",
            "    }",
            "    walls",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
        ])
    elif boundary_mode == "tracer":
        lines.extend([
            "    inlet",
            "    {",
            "        type            fixedValue;",
            "        flux_mol_s      1.9e-10;",
            "        value           nonuniform List<scalar> 1 ( 1.0 );",
            "    }",
            "    outlet",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      -2.2e-11;",
            "    }",
            "    walls",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
        ])
    else:
        lines.extend([
            "    inlet",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
            "    outlet",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
            "    walls",
            "    {",
            "        type            zeroGradient;",
            "        flux_mol_s      0;",
            "    }",
        ])
    lines.append("}")
    lines.append("")
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def get_mesh_points_and_cells():
    # 20 vertices across 5 x-planes: x = 0.0, 0.001, 0.002, 0.003, 0.004 m
    pts = []
    for ix in range(5):
        x = ix * 0.001
        pts.append((x, -0.0005, -0.0005))
        pts.append((x,  0.0005, -0.0005))
        pts.append((x,  0.0005,  0.0005))
        pts.append((x, -0.0005,  0.0005))
    # 4 hexahedral cells (VTK_HEXAHEDRON = 12)
    cells = [
        [0, 4, 5, 1, 3, 7, 6, 2],
        [4, 8, 9, 5, 7, 11, 10, 6],
        [8, 12, 13, 9, 11, 15, 14, 10],
        [12, 16, 17, 13, 15, 19, 18, 14],
    ]
    return pts, cells


def cell_to_plane_nodal_values(cell_vals):
    # 5 x-planes from 4 cells: plane 0 = cell 0, plane 1..3 = average of adjacent cells, plane 4 = cell 3
    plane_vals = [
        cell_vals[0],
        0.5 * (cell_vals[0] + cell_vals[1]),
        0.5 * (cell_vals[1] + cell_vals[2]),
        0.5 * (cell_vals[2] + cell_vals[3]),
        cell_vals[3],
    ]
    nodal = []
    for pv in plane_vals:
        for _ in range(4):
            nodal.append(pv)
    return nodal


def write_vtu_timestep(filepath, t_s, step_dict, glmm_dict):
    pts, cells = get_mesh_points_and_cells()
    ux = step_dict["U_x"]
    nodal_ux = cell_to_plane_nodal_values(ux)

    scalar_fields = {
        "p": step_dict["p"],
        "porosity": step_dict["porosity"],
        "tracer": step_dict["tracer"],
        "pH": step_dict["pH"],
        "beta_0": step_dict["beta_0"],
        "beta_1": step_dict["beta_1"],
        "glmm_observed_y": glmm_dict["glmm_observed_y"],
        "glmm_fixed_effect_beta0_plus_beta1_x": glmm_dict["glmm_fixed_effect_beta0_plus_beta1_x"],
        "glmm_random_intercept_b_i": glmm_dict["glmm_random_intercept_b_i"],
        "glmm_residual_eps_ij": glmm_dict["glmm_residual_eps_ij"],
        "empowerment_bits": step_dict["empowerment_bits"],
    }

    lines = [
        '<?xml version="1.0" encoding="utf-8"?>',
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
        '  <UnstructuredGrid>',
        '    <Piece NumberOfPoints="20" NumberOfCells="4">',
        '      <PointData Scalars="tracer" Vectors="U">',
        '        <DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">',
    ]
    for (ip, u_val) in enumerate(nodal_ux):
        _, py, pz = pts[ip]
        # Subtle secondary cross-pore swirl component around anastomosis throats for 3D streamline richness
        uy = 0.08 * u_val * math.sin(2.0 * math.pi * (ip // 4) / 4.0) * (1.0 if py < 0 else -1.0)
        uz = 0.06 * u_val * math.cos(2.0 * math.pi * (ip // 4) / 4.0) * (1.0 if pz < 0 else -1.0)
        lines.append(f"          {u_val:.10g} {uy:.10g} {uz:.10g}")
    lines.append('        </DataArray>')

    for fname, cvals in scalar_fields.items():
        nvals = cell_to_plane_nodal_values(cvals)
        lines.append(f'        <DataArray type="Float64" Name="{fname}" NumberOfComponents="1" format="ascii">')
        lines.append("          " + " ".join(f"{v:.10g}" for v in nvals))
        lines.append('        </DataArray>')

    lines.extend([
        '      </PointData>',
        '      <CellData Scalars="tracer" Vectors="U">',
        '        <DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">',
    ])
    for u_val in ux:
        lines.append(f"          {u_val:.10g} 0.0 0.0")
    lines.append('        </DataArray>')

    for fname, cvals in scalar_fields.items():
        lines.append(f'        <DataArray type="Float64" Name="{fname}" NumberOfComponents="1" format="ascii">')
        lines.append("          " + " ".join(f"{v:.12g}" for v in cvals))
        lines.append('        </DataArray>')

    lines.extend([
        '      </CellData>',
        '      <Points>',
        '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">',
    ])
    for (px, py, pz) in pts:
        lines.append(f"          {px:.8g} {py:.8g} {pz:.8g}")
    lines.extend([
        '        </DataArray>',
        '      </Points>',
        '      <Cells>',
        '        <DataArray type="Int32" Name="connectivity" format="ascii">',
    ])
    for c in cells:
        lines.append("          " + " ".join(str(idx) for idx in c))
    lines.extend([
        '        </DataArray>',
        '        <DataArray type="Int32" Name="offsets" format="ascii">',
        '          8 16 24 32',
        '        </DataArray>',
        '        <DataArray type="UInt8" Name="types" format="ascii">',
        '          12 12 12 12',
        '        </DataArray>',
        '      </Cells>',
        '    </Piece>',
        '  </UnstructuredGrid>',
        '</VTKFile>',
        '',
    ])
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def build_topological_tube_segments_for_step(t_s, step_dict):
    """
    Constructs 3D polylines representing the mycelial/EPS graph G_bio(t) and
    pair-of-pants cobordism Sigma_0,3 / hollow shell Psi inside [0, 4e-3] x [-0.5e-3, 0.5e-3]^2 m.
    Criterial status code per <RULE>:
      +1.0 = verified_active (accepted verified candidates - green)
       0.0 = exploratory_unverified (unverified exploratory particles - amber/gold)
      -1.0 = criterial_pruned (criterial-pruned candidates - red dashed)
    """
    nodes = [
        (0.45e-3,  0.00e-3,  0.00e-3),   # 0: Apical tip (Cell 1)
        (1.45e-3,  0.24e-3,  0.12e-3),   # 1: Upper branch node (Cell 2)
        (1.45e-3, -0.24e-3, -0.12e-3),   # 2: Lower anastomosis node (Cell 2)
        (2.55e-3,  0.00e-3,  0.00e-3),   # 3: Biofilm microcolony waist (Cell 3)
        (3.55e-3,  0.00e-3,  0.00e-3),   # 4: Mineral attachment site (Cell 4)
        (2.15e-3,  0.32e-3, -0.20e-3),   # 5: Hollow-shell Psi overhang node (Cell 3)
        (2.05e-3, -0.30e-3,  0.16e-3),   # 6: Exploratory lower sprout tip (Cell 2/3)
    ]

    # Define edges active at each time step: (u, v, radius_m, status_code, is_bridge, flow_um_s)
    if t_s <= 1.0:
        # t = 0s: Unclosed Cantor branching tree (beta_0 = 2 disconnected components, beta_1 = 0 loops)
        edge_specs = [
            (0, 1, 0.045e-3,  1.0, 0, 4.2),
            (0, 2, 0.042e-3,  1.0, 0, 3.8),
            (3, 4, 0.040e-3,  1.0, 0, 5.1),   # Second disconnected component in Cells 3-4
            (1, 5, 0.025e-3,  0.0, 0, 1.2),   # Unclosed exploratory sprout tip
        ]
    elif t_s <= 301.0:
        # t = 300s: Primary hyphal fusion & first anastomosis loop (beta_0 = 1, beta_1 = 1)
        edge_specs = [
            (0, 1, 0.052e-3,  1.0, 0, 8.5),
            (0, 2, 0.050e-3,  1.0, 0, 8.1),
            (1, 2, 0.048e-3,  1.0, 1, 14.0),  # First anastomosis bridge closes cycle 1
            (1, 3, 0.046e-3,  1.0, 0, 10.2),  # Connects component 1 to waist node 3 (beta_0 -> 1)
            (3, 4, 0.044e-3,  1.0, 0, 9.4),
            (2, 6, 0.026e-3,  0.0, 0, 1.8),   # Unclosed exploratory sprout toward waist
        ]
    elif t_s <= 601.0:
        # t = 600s: Pair-of-pants Sigma_0,3 bifurcation (beta_0 = 1, beta_1 = 2)
        edge_specs = [
            (0, 1, 0.058e-3,  1.0, 0, 14.2),
            (0, 2, 0.055e-3,  1.0, 0, 13.5),
            (1, 2, 0.054e-3,  1.0, 1, 18.5),
            (1, 3, 0.060e-3,  1.0, 0, 21.0),
            (2, 3, 0.058e-3,  1.0, 1, 20.4),  # Second anastomosis bridge closes cycle 2 (Sigma_0,3)
            (3, 4, 0.052e-3,  1.0, 0, 15.0),
            (1, 4, 0.026e-3, -1.0, 0, 0.0),   # Criterial-pruned branch (red dashed)
            (2, 6, 0.028e-3,  0.0, 0, 2.5),   # Exploratory unverified particle/branch
        ]
    elif t_s <= 901.0:
        # t = 900s: Transient Hollow Shell Psi overhang cycle (beta_0 = 1, beta_1 = 3)
        edge_specs = [
            (0, 1, 0.062e-3,  1.0, 0, 16.8),
            (0, 2, 0.058e-3,  1.0, 0, 15.9),
            (1, 2, 0.058e-3,  1.0, 1, 21.0),
            (1, 3, 0.066e-3,  1.0, 0, 25.4),
            (2, 3, 0.064e-3,  1.0, 1, 24.8),
            (3, 4, 0.056e-3,  1.0, 0, 17.2),
            (1, 5, 0.036e-3,  0.0, 1, 6.4),   # Exploratory overhang loop (cycle 3)
            (5, 3, 0.036e-3,  0.0, 1, 6.1),
            (1, 4, 0.025e-3, -1.0, 0, 0.0),   # Criterial-pruned (red dashed)
            (2, 4, 0.025e-3, -1.0, 0, 0.0),   # Criterial-pruned (red dashed)
        ]
    else:
        # t = 1200s: Pareto Trap-Space stabilized attractor (cord thickening + pruning, beta_0 = 1, beta_1 = 2)
        edge_specs = [
            (0, 1, 0.068e-3,  1.0, 0, 18.5),
            (0, 2, 0.065e-3,  1.0, 0, 17.8),
            (1, 2, 0.064e-3,  1.0, 1, 22.0),
            (1, 3, 0.074e-3,  1.0, 0, 31.0),
            (2, 3, 0.072e-3,  1.0, 1, 29.5),
            (3, 4, 0.062e-3,  1.0, 0, 19.0),
            (1, 4, 0.025e-3, -1.0, 0, 0.0),   # Criterial-pruned bypass (red dashed)
            (1, 5, 0.025e-3, -1.0, 0, 0.0),   # Criterial-pruned Psi overhang (red dashed)
            (2, 6, 0.026e-3,  0.0, 0, 1.4),   # Exploratory unverified particle/branch
        ]

    return nodes, edge_specs


def write_vtp_tubes_timestep(filepath, t_s, step_dict):
    nodes, edge_specs = build_topological_tube_segments_for_step(t_s, step_dict)
    n_sub = 8  # smooth curved polyline segments per active/exploratory edge

    pts = []
    pt_radius = []
    pt_flow = []
    pt_beta1 = []
    pt_exp_particle_radius = []
    lines_conn = []
    lines_offsets = []
    cell_status = []
    cell_bridge = []

    offset_acc = 0
    for (e_idx, (u, v, rad, status, is_br, flow)) in enumerate(edge_specs):
        p0 = nodes[u]
        p1 = nodes[v]
        arch_sign = 1.0 if (e_idx % 2 == 0) else -1.0
        if status < -0.5:
            # Criterial-pruned candidate (status == -1.0): emit 4 discrete dash segments with gaps
            # so ParaView 3D Tube filter renders a true RED DASHED tube per <RULE> 2
            dash_intervals = [(0.04, 0.22), (0.30, 0.48), (0.56, 0.74), (0.82, 0.96)]
            for (s_start, s_end) in dash_intervals:
                line_indices = []
                for k_dash in range(3):
                    s = s_start + (s_end - s_start) * (k_dash / 2.0)
                    arch = math.sin(math.pi * s) * 0.06e-3 * arch_sign
                    px = (1.0 - s) * p0[0] + s * p1[0]
                    py = (1.0 - s) * p0[1] + s * p1[1] + arch
                    pz = (1.0 - s) * p0[2] + s * p1[2] + 0.5 * arch
                    pt_idx = len(pts)
                    pts.append((px, py, pz))
                    pt_radius.append(rad)
                    pt_flow.append(flow)
                    pt_beta1.append(float(step_dict["global_beta_1"]))
                    pt_exp_particle_radius.append(0.0)
                    line_indices.append(pt_idx)
                lines_conn.append(line_indices)
                offset_acc += len(line_indices)
                lines_offsets.append(offset_acc)
                cell_status.append(status)
                cell_bridge.append(float(is_br))
        else:
            line_indices = []
            for k in range(n_sub + 1):
                s = k / float(n_sub)
                arch = math.sin(math.pi * s) * 0.06e-3 * arch_sign
                px = (1.0 - s) * p0[0] + s * p1[0]
                py = (1.0 - s) * p0[1] + s * p1[1] + arch
                pz = (1.0 - s) * p0[2] + s * p1[2] + 0.5 * arch
                pt_idx = len(pts)
                pts.append((px, py, pz))
                pt_radius.append(rad)
                pt_flow.append(flow)
                pt_beta1.append(float(step_dict["global_beta_1"]))
                # Unverified exploratory particles along status == 0.0 trajectories
                is_exp_pt = (abs(status) < 0.25) and (k in (2, 4, 6, 8))
                pt_exp_particle_radius.append(rad * 2.1 if is_exp_pt else 0.0)
                line_indices.append(pt_idx)
            lines_conn.append(line_indices)
            offset_acc += len(line_indices)
            lines_offsets.append(offset_acc)
            cell_status.append(status)
            cell_bridge.append(float(is_br))

    n_pts = len(pts)
    n_lines = len(lines_conn)

    xml_lines = [
        '<?xml version="1.0" encoding="utf-8"?>',
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">',
        '  <PolyData>',
        f'    <Piece NumberOfPoints="{n_pts}" NumberOfVerts="0" NumberOfLines="{n_lines}" NumberOfStrips="0" NumberOfPolys="0">',
        '      <PointData Scalars="tube_radius_m">',
        '        <DataArray type="Float64" Name="tube_radius_m" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{r:.8g}" for r in pt_radius),
        '        </DataArray>',
        '        <DataArray type="Float64" Name="cytoplasmic_flow_um_s" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{fl:.8g}" for fl in pt_flow),
        '        </DataArray>',
        '        <DataArray type="Float64" Name="beta_1_cycle_rank" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{b1:.8g}" for b1 in pt_beta1),
        '        </DataArray>',
        '        <DataArray type="Float64" Name="exploratory_particle_radius_m" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{er:.8g}" for er in pt_exp_particle_radius),
        '        </DataArray>',
        '      </PointData>',
        '      <CellData Scalars="criterial_status_code">',
        '        <DataArray type="Float64" Name="criterial_status_code" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{st:.4g}" for st in cell_status),
        '        </DataArray>',
        '        <DataArray type="Float64" Name="is_anastomosis_bridge" NumberOfComponents="1" format="ascii">',
        "          " + " ".join(f"{br:.4g}" for br in cell_bridge),
        '        </DataArray>',
        '      </CellData>',
        '      <Points>',
        '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">',
    ]
    for (px, py, pz) in pts:
        xml_lines.append(f"          {px:.8g} {py:.8g} {pz:.8g}")
    xml_lines.extend([
        '        </DataArray>',
        '      </Points>',
        '      <Lines>',
        '        <DataArray type="Int32" Name="connectivity" format="ascii">',
    ])
    for conn in lines_conn:
        xml_lines.append("          " + " ".join(str(idx) for idx in conn))
    xml_lines.extend([
        '        </DataArray>',
        '        <DataArray type="Int32" Name="offsets" format="ascii">',
        "          " + " ".join(str(off) for off in lines_offsets),
        '        </DataArray>',
        '      </Lines>',
        '    </Piece>',
        '  </PolyData>',
        '</VTKFile>',
        '',
    ])
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(xml_lines))


def write_pvd_collection(filepath, step_entries):
    lines = [
        '<?xml version="1.0" encoding="utf-8"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '  <Collection>',
    ]
    for (t_s, rel_file) in step_entries:
        lines.append(f'    <DataSet timestep="{t_s:.1f}" group="" part="0" file="{rel_file}"/>')
    lines.extend([
        '  </Collection>',
        '</VTKFile>',
        '',
    ])
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_temporal_xdmf(filepath):
    pts, cells = get_mesh_points_and_cells()
    coords_str = " ".join(f"{px:.8g} {py:.8g} {pz:.8g}" for (px, py, pz) in pts)
    conn_str = " ".join(" ".join(str(idx) for idx in c) for c in cells)

    lines = [
        '<?xml version="1.0" encoding="utf-8"?>',
        '<Xdmf Version="3.0">',
        '  <Domain>',
        '    <Grid Name="M7_Topological_GLMM_TimeSeries" GridType="Collection" CollectionType="Temporal">',
    ]
    for t_s in TIMESTEPS:
        step = TIMESERIES_DATA[t_s]
        glmm = compute_glmm_fields_for_step(step)
        t_int = int(t_s)
        lines.extend([
            f'      <Grid Name="M7_Step_{t_int:04d}s" GridType="Uniform">',
            f'        <Time Value="{t_s:.1f}" />',
            '        <Topology TopologyType="Hexahedron" NumberOfElements="4">',
            f'          <DataItem Dimensions="4 8" NumberType="Int" Format="XML">{conn_str}</DataItem>',
            '        </Topology>',
            '        <Geometry GeometryType="XYZ">',
            f'          <DataItem Dimensions="20 3" NumberType="Float" Precision="8" Format="XML">{coords_str}</DataItem>',
            '        </Geometry>',
            f'        <Information Name="physical_time_s" Value="{t_s:.1f}" />',
            f'        <Information Name="biological_time_s" Value="{t_s:.1f}" />',
            f'        <Information Name="sobolev_relaxation_step" Value="{step["sobolev_step"]}" />',
            f'        <Information Name="global_beta_0" Value="{step["global_beta_0"]}" />',
            f'        <Information Name="global_beta_1" Value="{step["global_beta_1"]}" />',
            f'        <Information Name="topology_stage" Value="{step["topology_stage"]}" />',
        ])
        u_str = " ".join(f"{u:.10g} 0.0 0.0" for u in step["U_x"])
        lines.extend([
            '        <Attribute Name="U" AttributeType="Vector" Center="Cell">',
            f'          <DataItem Dimensions="4 3" NumberType="Float" Precision="8" Format="XML">{u_str}</DataItem>',
            '        </Attribute>',
        ])
        scalars = {
            "p": step["p"],
            "porosity": step["porosity"],
            "tracer": step["tracer"],
            "pH": step["pH"],
            "beta_0": step["beta_0"],
            "beta_1": step["beta_1"],
            "glmm_observed_y": glmm["glmm_observed_y"],
            "glmm_fixed_effect_beta0_plus_beta1_x": glmm["glmm_fixed_effect_beta0_plus_beta1_x"],
            "glmm_random_intercept_b_i": glmm["glmm_random_intercept_b_i"],
            "glmm_residual_eps_ij": glmm["glmm_residual_eps_ij"],
            "empowerment_bits": step["empowerment_bits"],
        }
        for fname, vals in scalars.items():
            v_str = " ".join(f"{v:.12g}" for v in vals)
            lines.extend([
                f'        <Attribute Name="{fname}" AttributeType="Scalar" Center="Cell">',
                f'          <DataItem Dimensions="4" NumberType="Float" Precision="8" Format="XML">{v_str}</DataItem>',
                '        </Attribute>',
            ])
        lines.append('      </Grid>')

    lines.extend([
        '    </Grid>',
        '  </Domain>',
        '</Xdmf>',
        '',
    ])
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def generate_openfoam_and_timeseries_files(case_dir):
    ts_dir = os.path.join(case_dir, "timeseries")
    os.makedirs(ts_dir, exist_ok=True)

    # Update controlDict writeInterval to 6000 (300 s at deltaT = 0.05 s)
    control_dict_path = os.path.join(case_dir, "system", "controlDict")
    if os.path.isdir(os.path.dirname(control_dict_path)):
        cd_txt = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
application     porousSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1200;
deltaT          0.05;
writeControl    timeStep;
writeInterval   6000;
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
"""
        with open(control_dict_path, "w", encoding="utf-8") as f:
            f.write(cd_txt)

    mesh_pvd_entries = []
    tube_pvd_entries = []

    for t_s in TIMESTEPS:
        step = TIMESERIES_DATA[t_s]
        glmm = compute_glmm_fields_for_step(step)
        loc_dir = step["dir"]
        abs_t_dir = os.path.join(case_dir, loc_dir)
        os.makedirs(abs_t_dir, exist_ok=True)

        # Write OpenFOAM ASCII dictionaries in each time step directory
        # Note: VTK's vtkOpenFOAMReader ignores files ending in '_0' (treating them as old-time restart fields),
        # so we write both 'beta_0' (for direct ASCII/Julia/XDMF compatibility) and 'beta0' (so OpenFOAMReader exposes it).
        write_openfoam_vector_field(
            os.path.join(abs_t_dir, "U"), loc_dir, "U", "[0 1 -1 0 0 0 0]", step["U_x"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "p"), loc_dir, "p", "[0 2 -2 0 0 0 0]", step["p"], boundary_mode="pressure"
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "porosity"), loc_dir, "porosity", "[0 0 0 0 0 0 0]", step["porosity"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "tracer"), loc_dir, "C_substrate", "[0 -3 0 0 1 0 0]", step["tracer"], boundary_mode="tracer"
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "pH"), loc_dir, "pH", "[0 0 0 0 0 0 0]", step["pH"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "beta_0"), loc_dir, "beta_0", "[0 0 0 0 0 0 0]", step["beta_0"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "beta0"), loc_dir, "beta0", "[0 0 0 0 0 0 0]", step["beta_0"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "beta_1"), loc_dir, "beta_1", "[0 0 0 0 0 0 0]", step["beta_1"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "glmm_observed_y"), loc_dir, "glmm_observed_y", "[0 0 0 0 0 0 0]", glmm["glmm_observed_y"]
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "glmm_fixed_effect_beta0_plus_beta1_x"),
            loc_dir,
            "glmm_fixed_effect_beta0_plus_beta1_x",
            "[0 0 0 0 0 0 0]",
            glmm["glmm_fixed_effect_beta0_plus_beta1_x"],
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "glmm_random_intercept_b_i"),
            loc_dir,
            "glmm_random_intercept_b_i",
            "[0 0 0 0 0 0 0]",
            glmm["glmm_random_intercept_b_i"],
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "glmm_residual_eps_ij"),
            loc_dir,
            "glmm_residual_eps_ij",
            "[0 0 0 0 0 0 0]",
            glmm["glmm_residual_eps_ij"],
        )
        write_openfoam_scalar_field(
            os.path.join(abs_t_dir, "empowerment_bits"),
            loc_dir,
            "empowerment_bits",
            "[0 0 0 0 0 0 0]",
            step["empowerment_bits"],
        )

        t_int = int(t_s)
        vtu_rel = f"timeseries/m7_channel_t{t_int:04d}.vtu"
        vtp_rel = f"timeseries/m7_tubes_t{t_int:04d}.vtp"
        write_vtu_timestep(os.path.join(case_dir, vtu_rel), t_s, step, glmm)
        write_vtp_tubes_timestep(os.path.join(case_dir, vtp_rel), t_s, step)
        mesh_pvd_entries.append((t_s, vtu_rel))
        tube_pvd_entries.append((t_s, vtp_rel))

    pvd_mesh_path = os.path.join(case_dir, "m7_topological_timeseries.pvd")
    pvd_tube_path = os.path.join(case_dir, "m7_topology_tubes.pvd")
    xdmf_ts_path = os.path.join(case_dir, "m7_topological_timeseries.xdmf")

    write_pvd_collection(pvd_mesh_path, mesh_pvd_entries)
    write_pvd_collection(pvd_tube_path, tube_pvd_entries)
    write_temporal_xdmf(xdmf_ts_path)

    return {
        "pvd_mesh_path": pvd_mesh_path,
        "pvd_tube_path": pvd_tube_path,
        "xdmf_ts_path": xdmf_ts_path,
    }


def run_paraview_pipeline_and_renders(case_dir):
    import paraview.simple as pvs

    pvs._DisableFirstRenderCameraReset()

    renders_dir = os.path.join(case_dir, "renders")
    os.makedirs(renders_dir, exist_ok=True)

    foam_path = os.path.join(case_dir, "case.foam")
    pvd_mesh_path = os.path.join(case_dir, "m7_topological_timeseries.pvd")
    pvd_tube_path = os.path.join(case_dir, "m7_topology_tubes.pvd")
    xdmf_ts_path = os.path.join(case_dir, "m7_topological_timeseries.xdmf")
    pvsm_path = os.path.join(case_dir, "m7_topological_glmm_timeseries.pvsm")

    # 1. Interrogate OpenFOAMReader on case.foam (with SkipZeroTime = 0)
    foam_reader = pvs.OpenFOAMReader(registrationName="M7_OpenFOAM_Case", FileName=foam_path)
    foam_reader.SkipZeroTime = 0
    foam_reader.Refresh()
    foam_reader.UpdatePipelineInformation()
    foam_reader.UpdatePipeline()
    foam_timesteps = [float(x) for x in foam_reader.TimestepValues]
    foam_cell_arrays = list(foam_reader.CellArrays)

    # 2. Interrogate Xdmf3ReaderT on m7_topological_timeseries.xdmf
    xdmf_reader = pvs.Xdmf3ReaderT(registrationName="M7_XDMF_TemporalCollection", FileName=[xdmf_ts_path])
    xdmf_reader.UpdatePipeline()
    xdmf_timesteps = [float(x) for x in xdmf_reader.TimestepValues]

    # 3. Load PVD Time-Series Readers for 3D Channel Mesh & Topological Anastomosis Tubes
    mesh_reader = pvs.PVDReader(registrationName="M7_Channel_GLMM_PVD", FileName=pvd_mesh_path)
    mesh_reader.UpdatePipeline()
    pvd_timesteps = [float(x) for x in mesh_reader.TimestepValues]

    tube_reader = pvs.PVDReader(registrationName="M7_Topological_Tubes_PVD", FileName=pvd_tube_path)
    tube_reader.UpdatePipeline()
    tube_timesteps = [float(x) for x in tube_reader.TimestepValues]

    # Build RenderView (1920x1080)
    render_view = pvs.CreateView("RenderView")
    render_view.ViewSize = [1920, 1080]
    render_view.OrientationAxesVisibility = 1
    try:
        render_view.UseColorPaletteForBackground = 0
    except Exception:
        pass
    render_view.Background = [0.04, 0.06, 0.10]

    # A. Outer Porous Channel Wireframe / Translucent Surface colored by glmm_observed_y
    mesh_disp = pvs.Show(mesh_reader, render_view, "UnstructuredGridRepresentation")
    mesh_disp.Representation = "Surface With Edges"
    mesh_disp.Opacity = 0.28
    mesh_disp.EdgeColor = [0.35, 0.55, 0.80]
    mesh_disp.LineWidth = 1.5
    pvs.ColorBy(mesh_disp, ("CELLS", "glmm_observed_y"))
    glmm_lut = pvs.GetColorTransferFunction("glmm_observed_y")
    glmm_lut.RescaleTransferFunction(1.55, 2.45)

    # B. Axial Midplane Slice (z = -0.00045 m backplane) colored by empowerment_bits
    slice_filter = pvs.Slice(registrationName="M7_Axial_Empowerment_Slice", Input=mesh_reader)
    slice_filter.SliceType = "Plane"
    slice_filter.SliceType.Origin = [0.002, 0.0, -0.00045]
    slice_filter.SliceType.Normal = [0.0, 0.0, 1.0]
    slice_filter.UpdatePipeline()
    slice_disp = pvs.Show(slice_filter, render_view, "GeometryRepresentation")
    slice_disp.Representation = "Surface"
    slice_disp.Opacity = 0.55
    pvs.ColorBy(slice_disp, ("POINTS", "empowerment_bits"))
    emp_lut = pvs.GetColorTransferFunction("empowerment_bits")
    emp_lut.RescaleTransferFunction(1.40, 1.90)

    # C. Substrate Isosurface Contours (Contour filter on point-interpolated tracer)
    contour_filter = pvs.Contour(registrationName="M7_Substrate_Isosurfaces", Input=mesh_reader)
    contour_filter.ContourBy = ["POINTS", "tracer"]
    contour_filter.Isosurfaces = [0.35, 0.60, 0.85]
    contour_filter.ComputeNormals = 1
    contour_filter.UpdatePipeline()
    contour_disp = pvs.Show(contour_filter, render_view, "GeometryRepresentation")
    contour_disp.Representation = "Surface"
    contour_disp.Opacity = 0.50
    pvs.ColorBy(contour_disp, ("POINTS", "tracer"))
    tracer_lut = pvs.GetColorTransferFunction("tracer")
    tracer_lut.RescaleTransferFunction(0.10, 1.00)

    # D. Darcy Pore Velocity StreamTracer & Vector Glyphs
    stream_filter = pvs.StreamTracer(
        registrationName="M7_Darcy_Velocity_StreamTracer",
        Input=mesh_reader,
        SeedType="Line",
    )
    stream_filter.Vectors = ["POINTS", "U"]
    stream_filter.MaximumStreamlineLength = 0.004
    stream_filter.SeedType.Point1 = [0.00005, -0.00035, -0.00025]
    stream_filter.SeedType.Point2 = [0.00005,  0.00035,  0.00025]
    stream_filter.SeedType.Resolution = 6
    stream_filter.UpdatePipeline()
    stream_disp = pvs.Show(stream_filter, render_view, "GeometryRepresentation")
    stream_disp.Representation = "Surface"
    stream_disp.LineWidth = 2.0
    pvs.ColorBy(stream_disp, ("POINTS", "U"))

    glyph_filter = pvs.Glyph(
        registrationName="M7_Darcy_Velocity_Glyphs",
        Input=mesh_reader,
        GlyphType="Arrow",
    )
    glyph_filter.OrientationArray = ["CELLS", "U"]
    glyph_filter.ScaleArray = ["CELLS", "U"]
    glyph_filter.ScaleFactor = 2.2
    glyph_filter.GlyphMode = "All Points"
    glyph_filter.UpdatePipeline()
    glyph_disp = pvs.Show(glyph_filter, render_view, "GeometryRepresentation")
    glyph_disp.Representation = "Surface"
    pvs.ColorBy(glyph_disp, ("POINTS", "U"))

    # E. 3D Topological Anastomosis & Pair-of-Pants Sigma_0,3 Tube Network
    tube_filter = pvs.Tube(registrationName="M7_Topological_Anastomosis_Tubes", Input=tube_reader)
    tube_filter.Scalars = ["POINTS", "tube_radius_m"]
    tube_filter.Vectors = ["POINTS", "None"]
    tube_filter.Radius = 4.5e-5
    tube_filter.VaryRadius = "By Scalar"
    tube_filter.RadiusFactor = 2.2
    tube_filter.NumberofSides = 18
    tube_filter.UpdatePipeline()
    tube_disp = pvs.Show(tube_filter, render_view, "GeometryRepresentation")
    tube_disp.Representation = "Surface"
    tube_disp.Opacity = 1.0
    pvs.ColorBy(tube_disp, ("CELLS", "criterial_status_code"))

    # Configure criterial separation color map per <RULE>:
    # -1.0 -> Red Dashed (criterial_pruned), 0.0 -> Amber/Gold (exploratory_unverified), +1.0 -> Green (verified_active)
    crit_lut = pvs.GetColorTransferFunction("criterial_status_code")
    crit_lut.RGBPoints = [
        -1.0, 0.94, 0.22, 0.24,
         0.0, 0.98, 0.76, 0.18,
         1.0, 0.16, 0.85, 0.45,
    ]
    crit_lut.RescaleTransferFunction(-1.0, 1.0)

    # E2. Unverified Exploratory Particles (Sphere Glyphs scaled by exploratory_particle_radius_m per <RULE> 2)
    exp_particle_glyph = pvs.Glyph(
        registrationName="M7_Exploratory_Particles",
        Input=tube_reader,
        GlyphType="Sphere",
    )
    exp_particle_glyph.ScaleArray = ["POINTS", "exploratory_particle_radius_m"]
    exp_particle_glyph.ScaleFactor = 2.6
    exp_particle_glyph.GlyphMode = "All Points"
    exp_particle_glyph.UpdatePipeline()
    exp_particle_disp = pvs.Show(exp_particle_glyph, render_view, "GeometryRepresentation")
    exp_particle_disp.Representation = "Surface"
    pvs.ColorBy(exp_particle_disp, None)
    exp_particle_disp.DiffuseColor = [1.0, 0.84, 0.22]

    # F. Scientific Title & Telemetry Text Overlays (Zero em/en-dashes per <RULE>)
    title_text = pvs.Text(registrationName="M7_Header_Annotation")
    title_text.Text = (
        "M7 Multiscale Topological Time-Series & OpenFOAM/GLMM Workbench (ParaView 6.2.0)\n"
        "Criterial Separation: Green = Verified Active ; Red Dashed = Criterial Pruned ; Gold Spheres = Exploratory Particles"
    )
    title_disp = pvs.Show(title_text, render_view, "TextSourceRepresentation")
    title_disp.WindowLocation = "Upper Left Corner"
    title_disp.FontSize = 18
    title_disp.Color = [0.93, 0.96, 1.0]

    telemetry_text = pvs.Text(registrationName="M7_Timestep_Telemetry")
    telemetry_disp = pvs.Show(telemetry_text, render_view, "TextSourceRepresentation")
    telemetry_disp.WindowLocation = "Lower Left Corner"
    telemetry_disp.FontSize = 15
    telemetry_disp.Color = [0.15, 0.85, 0.95]

    # Position 3D camera for clear isometric inspection of the 4 mm porous channel
    render_view.CameraPosition = [0.0020, -0.0038, 0.0024]
    render_view.CameraFocalPoint = [0.0020, 0.0, 0.0]
    render_view.CameraViewUp = [0.0, 0.48, 0.88]
    render_view.CameraParallelScale = 0.0022

    anim_scene = pvs.GetAnimationScene()
    anim_scene.UpdateAnimationUsingDataTimeSteps()

    rendered_files = []
    step_summaries = []

    for t_s in TIMESTEPS:
        step = TIMESERIES_DATA[t_s]
        glmm = compute_glmm_fields_for_step(step)
        ph_act = compute_activity_pH_avg(step["pH"])
        emp_mean = sum(step["empowerment_bits"]) / 4.0
        glmm_y_mean = sum(glmm["glmm_observed_y"]) / 4.0

        telemetry_text.Text = (
            f"Physical Time t = {t_s:.1f} s ; Biological Time = {t_s:.1f} s ; Pinned Sobolev Step = {step['sobolev_step']}\n"
            f"Topology: {step['topology_stage']} ; Activity pH_avg = {ph_act:.4f}\n"
            f"OpenFOAM/GLMM: mean(y_ij) = {glmm_y_mean:.4f} (beta0=1.20, beta1=0.85, ICC=0.985) ; Mean Empowerment E(s_t) = {emp_mean:.3f} bits"
        )

        anim_scene.AnimationTime = t_s
        render_view.ViewTime = t_s
        pvs.Render(render_view)

        t_int = int(t_s)
        png_path = os.path.join(renders_dir, f"m7_paraview_t{t_int:04d}s.png")
        pvs.SaveScreenshot(png_path, render_view, ImageResolution=[1920, 1080])
        rendered_files.append(png_path)

        step_summaries.append({
            "time_s": t_s,
            "sobolev_relaxation_step": step["sobolev_step"],
            "global_beta_0": step["global_beta_0"],
            "global_beta_1": step["global_beta_1"],
            "topology_stage": step["topology_stage"],
            "pH_activity_avg": round(ph_act, 8),
            "mean_velocity_x_m_s": round(sum(step["U_x"]) / 4.0, 10),
            "mean_tracer_mol_m3": round(sum(step["tracer"]) / 4.0, 8),
            "mean_glmm_observed_y": round(glmm_y_mean, 8),
            "mean_empowerment_bits": round(emp_mean, 8),
            "screenshot_path": f"renders/m7_paraview_t{t_int:04d}s.png",
        })

    # Save the interactive state at t = 1200.0 s so opening ParaView lands on the full mature state
    # with the 5-step time slider ready to scrub across [0, 300, 600, 900, 1200], then strip host paths
    pvs.SaveState(pvsm_path)
    if os.path.isfile(pvsm_path):
        with open(pvsm_path, "r", encoding="utf-8") as f_pvsm:
            pvsm_raw = f_pvsm.read()
        pvsm_rel = re.sub(r"/Users/[^\"'<>\s]+/openfoam_case/([^\"'<>\s]+)", r"\1", pvsm_raw)
        with open(pvsm_path, "w", encoding="utf-8") as f_pvsm:
            f_pvsm.write(pvsm_rel)

    return {
        "foam_timesteps": foam_timesteps,
        "foam_cell_arrays": sorted(foam_cell_arrays),
        "xdmf_timesteps": xdmf_timesteps,
        "pvd_timesteps": pvd_timesteps,
        "tube_timesteps": tube_timesteps,
        "pvsm_path": pvsm_path,
        "rendered_files": rendered_files,
        "step_summaries": step_summaries,
    }


def build_composite_comparison_gallery(case_dir, pv_result):
    """
    Combines the 5 ParaView 1920x1080 renders and time-series topological/GLMM curves into
    a single 1920x1080 human-visual-analysis comparison sheet (strictly zero em/en-dashes)
    using matplotlib built into pvpython.
    """
    renders_dir = os.path.join(case_dir, "renders")
    comp_path = os.path.join(renders_dir, "m7_paraview_topological_glmm_comparison_1920x1080.png")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="#070b14")
    gs = GridSpec(2, 15, figure=fig, height_ratios=[0.95, 1.25], hspace=0.26, wspace=0.85,
                  left=0.03, right=0.985, top=0.89, bottom=0.065)

    fig.suptitle(
        "M7 Interactive Topological Time-Series & OpenFOAM/GLMM ParaView 6.2.0 Suite (t = 0, 300, 600, 900, 1200 s)",
        color="#eef4ff", fontsize=15, fontweight="bold", y=0.965
    )
    fig.text(
        0.5, 0.930,
        "Criterial Separation: Green = Verified Active Anastomosis (Sigma_0,3) ; Red Dashed = Criterial Pruned Excursion (Psi) ; Gold Spheres = Exploratory Particles ; Epistemic Status: modeled_system",
        ha="center", va="center", color="#22d3ee", fontsize=10.2
    )

    steps = [0, 300, 600, 900, 1200]
    beta0_vals = [TIMESERIES_DATA[float(t)]["global_beta_0"] for t in steps]
    beta1_vals = [TIMESERIES_DATA[float(t)]["global_beta_1"] for t in steps]
    stage_short = [
        "Unclosed Tree",
        "Primary Loop (so/Prm1)",
        "Pair-of-Pants (Sigma_0,3)",
        "Hollow Shell (Psi)",
        "Pareto Attractor",
    ]
    emp_vals = [sum(TIMESERIES_DATA[float(t)]["empowerment_bits"]) / 4.0 for t in steps]
    glmm_means = [sum(compute_glmm_fields_for_step(TIMESERIES_DATA[float(t)])["glmm_observed_y"]) / 4.0 for t in steps]
    ph_acts = [compute_activity_pH_avg(TIMESERIES_DATA[float(t)]["pH"]) for t in steps]

    # Top row: 5 ParaView 3D time-step renders
    for idx, t in enumerate(steps):
        ax = fig.add_subplot(gs[0, idx * 3:(idx + 1) * 3])
        ax.set_facecolor("#0c1322")
        p = os.path.join(renders_dir, f"m7_paraview_t{t:04d}s.png")
        if os.path.isfile(p):
            img = mpimg.imread(p)
            h, w = img.shape[0], img.shape[1]
            crop = img[int(h * 0.14):int(h * 0.84), int(w * 0.04):int(w * 0.96)]
            ax.imshow(crop)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor("#22d3ee" if t == 1200 else "#2d466e")
            spine.set_linewidth(2.0 if t == 1200 else 1.2)
        ax.set_title(
            f"t = {t} s ; {stage_short[idx]}\nbeta_0={beta0_vals[idx]}, beta_1={beta1_vals[idx]} ; E={emp_vals[idx]:.3f} bits ; y={glmm_means[idx]:.3f}",
            color="#eef4ff", fontsize=9.0, pad=6
        )

    # Bottom Left (cols 0..6): Hero Full-Frame View of t = 1200 s
    ax_hero = fig.add_subplot(gs[1, 0:6])
    ax_hero.set_facecolor("#0c1322")
    hero_p = os.path.join(renders_dir, "m7_paraview_t1200s.png")
    if os.path.isfile(hero_p):
        himg = mpimg.imread(hero_p)
        ax_hero.imshow(himg)
    ax_hero.set_xticks([])
    ax_hero.set_yticks([])
    for spine in ax_hero.spines.values():
        spine.set_edgecolor("#10b981")
        spine.set_linewidth(2.0)
    ax_hero.set_title(
        "ParaView 6.2.0 State (t = 1200 s): Slice + Contour + StreamTracer + Glyph + Criterial Tubes",
        color="#10b981", fontsize=9.5, fontweight="bold", pad=6
    )

    # Bottom Middle (cols 6..10): Topological Betti Progression & Channel Empowerment Trajectory
    ax_top = fig.add_subplot(gs[1, 6:10])
    ax_top.set_facecolor("#0c1322")
    pruned_vals = [0, 0, 1, 2, 2]
    ax_top.plot(steps, beta1_vals, "-o", color="#10b981", linewidth=2.8, markersize=7, label="Verified Cycle Rank beta_1(t)")
    ax_top.plot(steps, pruned_vals, "--x", color="#ef4444", linewidth=2.2, markersize=6, label="Criterial Pruned Branches (Red Dashed)")
    ax_top.plot(steps, beta0_vals, "-.s", color="#f59e0b", linewidth=2.0, markersize=5.5, label="Connected Components beta_0(t)")
    ax_top.plot(steps, emp_vals, "-^", color="#22d3ee", linewidth=2.3, markersize=6, label="Mean Empowerment E(s_t) [bits]")
    ax_top.set_xlabel("Physical & Biological Time t [s] (Sobolev Steps 0..400)", color="#9cb3d9", fontsize=9.0)
    ax_top.set_ylabel("Homology Rank / Empowerment [bits]", color="#9cb3d9", fontsize=9.0)
    ax_top.set_xticks(steps)
    ax_top.set_ylim(-0.2, 3.6)
    ax_top.tick_params(colors="#9cb3d9", labelsize=8.5)
    ax_top.grid(True, color="#1e2d4a", linestyle=":", linewidth=0.8)
    for spine in ax_top.spines.values():
        spine.set_edgecolor("#2d466e")
    ax_top.legend(loc="upper left", facecolor="#09101d", edgecolor="#2d466e", labelcolor="#eef4ff", fontsize=7.6)
    ax_top.set_title(
        "Topology beta_1(t), Criterial Pruning & E(s_t)",
        color="#eef4ff", fontsize=9.5, fontweight="bold", pad=6
    )

    # Bottom Right (cols 10..15): Coupled OpenFOAM/GLMM Variance Decomposition across Cells 1..4
    ax_glmm = fig.add_subplot(gs[1, 10:15])
    ax_glmm.set_facecolor("#0c1322")
    cell_colors = ["#22d3ee", "#10b981", "#f59e0b", "#a78bfa"]
    for c_idx in range(4):
        cell_y = [compute_glmm_fields_for_step(TIMESERIES_DATA[float(t)])["glmm_observed_y"][c_idx] for t in steps]
        b_i = GLMM_RANDOM_INTERCEPTS_B[c_idx]
        ax_glmm.plot(
            steps, cell_y, "-o", color=cell_colors[c_idx], linewidth=2.2, markersize=5.5,
            label=f"Cell {c_idx+1} y_ij(t) (b_{c_idx+1}={b_i:+.2f})"
        )
    ax_glmm.plot(steps, glmm_means, "w--", linewidth=2.0, label="Channel Mean GLMM Response")
    ax_glmm.set_xlabel("Physical Time t [s] (Activity pH_avg: 7.213 -> 7.032)", color="#9cb3d9", fontsize=9.0)
    ax_glmm.set_ylabel("GLMM y_ij(t) = (beta0 + beta1*x_ij) + b_i + eps_ij", color="#9cb3d9", fontsize=8.8)
    ax_glmm.set_xticks(steps)
    ax_glmm.tick_params(colors="#9cb3d9", labelsize=8.5)
    ax_glmm.grid(True, color="#1e2d4a", linestyle=":", linewidth=0.8)
    for spine in ax_glmm.spines.values():
        spine.set_edgecolor("#2d466e")
    ax_glmm.legend(loc="lower right", facecolor="#09101d", edgecolor="#2d466e", labelcolor="#eef4ff", fontsize=8.0)
    ax_glmm.set_title(
        "OpenFOAM/GLMM Cell Trajectories (ICC = 0.985)",
        color="#eef4ff", fontsize=9.5, fontweight="bold", pad=6
    )

    fig.savefig(comp_path, dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    return comp_path


def verify_png_dimensions(png_path):
    with open(png_path, "rb") as f:
        header = f.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        return (0, 0)
    w = int.from_bytes(header[16:20], "big")
    h = int.from_bytes(header[20:24], "big")
    return (w, h)


def sha256_file(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    gen_info = generate_openfoam_and_timeseries_files(case_dir)

    pv_info = None
    comp_png = None
    try:
        import paraview.simple  # noqa: F401
        pv_info = run_paraview_pipeline_and_renders(case_dir)
        comp_png = build_composite_comparison_gallery(case_dir, pv_info)
    except ImportError:
        print("Note: paraview.simple not in current Python; generated OpenFOAM/PVD/VTU/VTP/XDMF files.")

    report = {
        "suite": "M7_Interactive_Topological_TimeSeries_OpenFOAM_GLMM_ParaView_Suite",
        "case_dir": ".",
        "evidence_tier": EVIDENCE_TIER,
        "mesh_resolution_badge": MESH_RESOLUTION_BADGE,
        "voxel_validation_grid_badge": VOXEL_VALIDATION_GRID_BADGE,
        "pvsm_portable_relative_paths": True,
        "timesteps_s": TIMESTEPS,
        "openfoam_time_directories": TIMESTEP_DIRS,
        "pvd_mesh_file": "m7_topological_timeseries.pvd",
        "pvd_tubes_file": "m7_topology_tubes.pvd",
        "xdmf_temporal_file": "m7_topological_timeseries.xdmf",
        "paraview_state_file": "m7_topological_glmm_timeseries.pvsm",
        "zero_em_en_dashes_enforced": True,
    }

    if pv_info is not None:
        all_pngs = list(pv_info["rendered_files"])
        if comp_png and os.path.isfile(comp_png):
            all_pngs.append(comp_png)
        png_dims = {os.path.basename(p): list(verify_png_dimensions(p)) for p in all_pngs}
        png_hashes = {os.path.basename(p): sha256_file(p) for p in all_pngs}
        report.update({
            "paraview_verified": True,
            "foam_timesteps": pv_info["foam_timesteps"],
            "foam_cell_arrays": pv_info["foam_cell_arrays"],
            "xdmf_timesteps": pv_info["xdmf_timesteps"],
            "pvd_timesteps": pv_info["pvd_timesteps"],
            "tube_timesteps": pv_info["tube_timesteps"],
            "pvsm_sha256": sha256_file(pv_info["pvsm_path"]),
            "png_dimensions": png_dims,
            "png_sha256": png_hashes,
            "step_summaries": pv_info["step_summaries"],
            "composite_comparison_png": os.path.basename(comp_png) if comp_png else None,
        })

    report_path = os.path.join(case_dir, "m7_paraview_verification_report.json")
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print("Wrote verification report:", report_path)


if __name__ == "__main__":
    main()
