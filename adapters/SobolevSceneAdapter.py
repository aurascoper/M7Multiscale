"""M7 multiscale scene adapter bridging SobolevMacromolecule to M7SpatialFrame v1.0.0.

This module keeps SobolevMacromolecule headless while connecting its SMACRO1
binary coordinate frames, (N, 4, 4) instance matrices, far-field coarse-graining,
and OpenUSD PointInstancer exports to the M7 continuum coordinate/time/provenance
architecture (OpenFOAM + PHREEQC + Julia BiofilmMultiphysics).

Key invariants enforced:
1. Nested spatial transform:
       x_world = T_cell * T_molecule * x_local
   with explicit unit conversion (angstrom, nm, um, m -> m), proper rotation
   validation (SO(3), det(R) = +1), and exact inverse mapping.
2. Three distinct temporal processes (physical_time_s, biological_time_s,
   sobolev_relaxation_step) that are never conflated.
3. Four-scale microscope epistemic disclosure badges distinguishing direct solver
   outputs, experimental imaging reconstructions, external structural database
   templates, coarse-grained approximations, and conceptual/synthetic elements.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

import sobolev_visualization as sv


M7_SCHEMA_NAME = "M7SpatialFrame"
M7_SCHEMA_VERSION = "1.0.0"
M7_ADAPTER_VERSION = "1.0.0"

UNIT_TO_METERS: dict[str, float] = {
    "angstrom": 1.0e-10,
    "A": 1.0e-10,
    "nm": 1.0e-9,
    "nanometer": 1.0e-9,
    "um": 1.0e-6,
    "micrometer": 1.0e-6,
    "mm": 1.0e-3,
    "m": 1.0,
    "meter": 1.0,
}

VALID_PRODUCERS = ("OpenFOAM", "PHREEQC", "Julia", "Sobolev")
VALID_ASSOCIATIONS = ("cell", "point", "face", "molecular_instance", "graph_edge")
VALID_TOPOLOGIES = (
    "structured_cartesian",
    "unstructured",
    "anisotropic_hexahedral",
    "cut_cell_porous",
)
VALID_EVIDENCE_CLASSES = (
    "measured",
    "source_reported",
    "calculated",
    "inferred",
    "assumed",
    "synthetic",
    "unknown",
)
VALID_MICROSCOPE_SCALES = (
    "planetary",
    "geological",
    "biofilm",
    "molecular",
)
VALID_DISCLOSURE_CLASSES = (
    "direct_numerical_solver",
    "experimental_imaging_reconstruction",
    "external_structural_database",
    "coarse_grained_approximation",
    "conceptual_or_synthetic",
)


def sha256_prefixed(payload: bytes | str) -> str:
    """Return a canonical 'sha256:<64-lowercase-hex>' digest over exact bytes."""
    raw = payload.encode("utf-8") if isinstance(payload, str) else bytes(payload)
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def verify_sha256_prefixed(digest: str, label: str = "digest") -> str:
    """Validate that a digest string has the canonical 'sha256:<64hex>' format."""
    if not isinstance(digest, str) or not digest.startswith("sha256:"):
        raise ValueError(f"{label} must start with 'sha256:'")
    hex_part = digest[7:]
    if len(hex_part) != 64 or any(c not in "0123456789abcdef" for c in hex_part):
        raise ValueError(f"{label} must contain 64 lowercase hexadecimal characters after 'sha256:'")
    return digest


@dataclass(frozen=True)
class M7MultirateClockState:
    """Three-clock temporal state for M7 4D frames.

    Physical time (OpenFOAM/PHREEQC), biological time (Julia growth/remodeling),
    and Sobolev structural relaxation step (integer optimization step) are stored
    on separate axes and must never be conflated.
    """

    physical_time_s: float
    biological_time_s: float
    sobolev_relaxation_step: int
    snapshot_sync_id: str
    interpolation_valid: bool = True
    parent_snapshot_id: str | None = None

    def __post_init__(self) -> None:
        if isinstance(self.physical_time_s, bool) or not np.isfinite(self.physical_time_s) or self.physical_time_s < 0.0:
            raise ValueError("physical_time_s must be a finite non-negative float")
        if isinstance(self.biological_time_s, bool) or not np.isfinite(self.biological_time_s) or self.biological_time_s < 0.0:
            raise ValueError("biological_time_s must be a finite non-negative float")
        if isinstance(self.sobolev_relaxation_step, bool) or not isinstance(self.sobolev_relaxation_step, (int, np.integer)):
            raise ValueError("sobolev_relaxation_step must be an integer step counter, not a float time")
        if int(self.sobolev_relaxation_step) < 0:
            raise ValueError("sobolev_relaxation_step must be non-negative")
        if not isinstance(self.snapshot_sync_id, str) or not self.snapshot_sync_id.strip():
            raise ValueError("snapshot_sync_id must be a non-empty string")
        if self.parent_snapshot_id is not None:
            if not isinstance(self.parent_snapshot_id, str) or not self.parent_snapshot_id.strip():
                raise ValueError("parent_snapshot_id must be None or a non-empty string")
            if self.parent_snapshot_id == self.snapshot_sync_id:
                raise ValueError("parent_snapshot_id cannot equal snapshot_sync_id")

    def to_dict(self) -> dict[str, Any]:
        return {
            "physical_time_s": float(self.physical_time_s),
            "biological_time_s": float(self.biological_time_s),
            "sobolev_relaxation_step": int(self.sobolev_relaxation_step),
            "snapshot_sync_id": self.snapshot_sync_id,
            "interpolation_valid": bool(self.interpolation_valid),
            "parent_snapshot_id": self.parent_snapshot_id,
        }


def validate_temporal_interpolation(
    clock_a: M7MultirateClockState,
    clock_b: M7MultirateClockState,
    alpha: float,
    *,
    geometry_hash_a: str,
    geometry_hash_b: str,
    topology_changed: bool = False,
) -> dict[str, Any]:
    """Check whether interpolating between two 4D snapshots is physically valid.

    Refuses to interpolate when:
    - either frame has interpolation_valid == False;
    - mesh geometry hashes differ;
    - discrete biological/molecular graph topology changed;
    - physical_time_s did not advance while sobolev_relaxation_step changed
      (preventing conflation of optimization iterations with physical time).
    """

    if not (0.0 <= float(alpha) <= 1.0) or not np.isfinite(alpha):
        raise ValueError("interpolation parameter alpha must lie in [0.0, 1.0]")

    reasons: list[str] = []
    if not clock_a.interpolation_valid or not clock_b.interpolation_valid:
        reasons.append("interpolation_valid_flag_false")
    if geometry_hash_a != geometry_hash_b:
        reasons.append("mesh_geometry_hash_mismatch")
    if topology_changed:
        reasons.append("discrete_graph_topology_discontinuity")

    dt_phys = float(clock_b.physical_time_s - clock_a.physical_time_s)
    dt_bio = float(clock_b.biological_time_s - clock_a.biological_time_s)
    dstep_sob = int(clock_b.sobolev_relaxation_step - clock_a.sobolev_relaxation_step)

    if dt_phys < 0.0 or dt_bio < 0.0:
        reasons.append("non_monotonic_physical_or_biological_clock")
    if abs(dt_phys) <= 1.0e-15 and abs(dt_bio) <= 1.0e-15 and dstep_sob != 0:
        reasons.append("sobolev_relaxation_steps_cannot_be_interpolated_as_physical_time")

    allowed = len(reasons) == 0
    interp_phys = float(clock_a.physical_time_s + alpha * dt_phys) if allowed else None
    interp_bio = float(clock_a.biological_time_s + alpha * dt_bio) if allowed else None
    pinned_sob_step = (
        int(clock_a.sobolev_relaxation_step)
        if alpha < 0.5
        else int(clock_b.sobolev_relaxation_step)
    )

    return {
        "interpolation_allowed": allowed,
        "blocking_reasons": reasons,
        "physical_time_s": interp_phys,
        "biological_time_s": interp_bio,
        "pinned_sobolev_relaxation_step": pinned_sob_step,
        "alpha": float(alpha),
    }


def _validate_so3_rotation(rot: np.ndarray, name: str = "rotation", atol: float = 1.0e-5) -> np.ndarray:
    arr = np.asarray(rot, dtype=np.float64)
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        raise ValueError(f"{name} must be a finite (3, 3) matrix")
    ident = arr.T @ arr
    if not np.allclose(ident, np.eye(3), atol=atol):
        raise ValueError(f"{name} must be orthogonal (R^T R = I)")
    det = float(np.linalg.det(arr))
    if not np.isclose(det, 1.0, atol=atol):
        raise ValueError(f"{name} must be a proper right-handed rotation with det(R) = +1, got {det:.6f}")
    return arr


@dataclass(frozen=True)
class M7NestedTransform:
    """Nested continuum-to-molecular transform: x_world = T_cell * T_molecule * x_local.

    T_cell places the continuum cell centroid and local axes in the right-handed
    world coordinate frame (meters).
    T_molecule scales local molecular coordinates (e.g. Angstroms or nanometers)
    into meters, applies the molecular orientation R_molecule, and offsets the
    molecular anchor within the cell.
    """

    cell_id: int
    T_cell: np.ndarray
    T_molecule: np.ndarray
    local_unit: str = "angstrom"

    def __post_init__(self) -> None:
        if isinstance(self.cell_id, bool) or int(self.cell_id) < 1:
            raise ValueError("cell_id must be >= 1 (1-based M7 continuum cell index)")
        if self.local_unit not in UNIT_TO_METERS:
            raise ValueError(f"unsupported local_unit '{self.local_unit}', expected one of {sorted(UNIT_TO_METERS)}")

        t_cell = np.asarray(self.T_cell, dtype=np.float64)
        t_mol = np.asarray(self.T_molecule, dtype=np.float64)
        if t_cell.shape != (4, 4) or not np.isfinite(t_cell).all():
            raise ValueError("T_cell must be a finite (4, 4) homogeneous matrix")
        if t_mol.shape != (4, 4) or not np.isfinite(t_mol).all():
            raise ValueError("T_molecule must be a finite (4, 4) homogeneous matrix")
        if not np.allclose(t_cell[3, :], [0.0, 0.0, 0.0, 1.0], atol=1.0e-9):
            raise ValueError("T_cell bottom row must be [0, 0, 0, 1]")
        if not np.allclose(t_mol[3, :], [0.0, 0.0, 0.0, 1.0], atol=1.0e-9):
            raise ValueError("T_molecule bottom row must be [0, 0, 0, 1]")

        _validate_so3_rotation(t_cell[:3, :3], name="T_cell[:3, :3]")

        expected_scale = UNIT_TO_METERS[self.local_unit]
        mol_upper = t_mol[:3, :3]
        if expected_scale <= 0.0:
            raise ValueError("unit scale must be positive")
        normalized_rot = mol_upper / expected_scale
        _validate_so3_rotation(normalized_rot, name="T_molecule[:3, :3] / unit_scale")

    @classmethod
    def from_cell_and_molecule(
        cls,
        cell_id: int,
        cell_center_m: Sequence[float] | np.ndarray,
        *,
        cell_rotation: Sequence[Sequence[float]] | np.ndarray | None = None,
        molecule_offset_m: Sequence[float] | np.ndarray = (0.0, 0.0, 0.0),
        molecule_rotation: Sequence[Sequence[float]] | np.ndarray | None = None,
        local_unit: str = "angstrom",
    ) -> "M7NestedTransform":
        if local_unit not in UNIT_TO_METERS:
            raise ValueError(f"unsupported local_unit '{local_unit}'")
        c_pos = np.asarray(cell_center_m, dtype=np.float64)
        if c_pos.shape != (3,) or not np.isfinite(c_pos).all():
            raise ValueError("cell_center_m must be a finite 3-vector")
        m_off = np.asarray(molecule_offset_m, dtype=np.float64)
        if m_off.shape != (3,) or not np.isfinite(m_off).all():
            raise ValueError("molecule_offset_m must be a finite 3-vector")

        r_cell = np.eye(3, dtype=np.float64) if cell_rotation is None else _validate_so3_rotation(np.asarray(cell_rotation, dtype=np.float64), "cell_rotation")
        r_mol = np.eye(3, dtype=np.float64) if molecule_rotation is None else _validate_so3_rotation(np.asarray(molecule_rotation, dtype=np.float64), "molecule_rotation")

        scale = UNIT_TO_METERS[local_unit]
        t_cell = np.eye(4, dtype=np.float64)
        t_cell[:3, :3] = r_cell
        t_cell[:3, 3] = c_pos

        t_mol = np.eye(4, dtype=np.float64)
        t_mol[:3, :3] = scale * r_mol
        t_mol[:3, 3] = m_off

        return cls(
            cell_id=int(cell_id),
            T_cell=t_cell,
            T_molecule=t_mol,
            local_unit=local_unit,
        )

    @property
    def unit_scale_m(self) -> float:
        return UNIT_TO_METERS[self.local_unit]

    @property
    def T_world_from_local(self) -> np.ndarray:
        """Return the composite 4x4 matrix T_cell @ T_molecule."""
        return np.asarray(self.T_cell, dtype=np.float64) @ np.asarray(self.T_molecule, dtype=np.float64)

    @property
    def T_local_from_world(self) -> np.ndarray:
        """Return the analytical inverse 4x4 matrix T_molecule^(-1) @ T_cell^(-1)."""
        t_cell = np.asarray(self.T_cell, dtype=np.float64)
        t_mol = np.asarray(self.T_molecule, dtype=np.float64)
        r_cell_t = t_cell[:3, :3].T
        inv_cell = np.eye(4, dtype=np.float64)
        inv_cell[:3, :3] = r_cell_t
        inv_cell[:3, 3] = -(r_cell_t @ t_cell[:3, 3])

        s = self.unit_scale_m
        r_mol = t_mol[:3, :3] / s
        r_mol_t_scaled = r_mol.T / s
        inv_mol = np.eye(4, dtype=np.float64)
        inv_mol[:3, :3] = r_mol_t_scaled
        inv_mol[:3, 3] = -(r_mol_t_scaled @ t_mol[:3, 3])
        return inv_mol @ inv_cell

    def transform_local_to_world(self, coords_local: np.ndarray | Sequence[Sequence[float]]) -> np.ndarray:
        """Map local molecular coordinates (N, 3) to world coordinates (N, 3) in meters."""
        arr = sv.validate_coords(coords_local, name="coords_local").astype(np.float64, copy=False)
        t_cell = np.asarray(self.T_cell, dtype=np.float64)
        t_mol = np.asarray(self.T_molecule, dtype=np.float64)
        in_cell = (arr @ t_mol[:3, :3].T) + t_mol[:3, 3][None, :]
        return (in_cell @ t_cell[:3, :3].T) + t_cell[:3, 3][None, :]

    def transform_world_to_local(self, coords_world: np.ndarray | Sequence[Sequence[float]]) -> np.ndarray:
        """Map world coordinates (N, 3) in meters back to local molecular coordinates (N, 3).

        Subtracts continuum cell translation and molecular offset in meters before
        dividing by the Angstrom/nm unit scale to avoid catastrophic cancellation.
        """
        arr = sv.validate_coords(coords_world, name="coords_world").astype(np.float64, copy=False)
        t_cell = np.asarray(self.T_cell, dtype=np.float64)
        t_mol = np.asarray(self.T_molecule, dtype=np.float64)
        s = self.unit_scale_m
        r_mol = t_mol[:3, :3] / s
        in_cell = (arr - t_cell[:3, 3][None, :]) @ t_cell[:3, :3]
        return ((in_cell - t_mol[:3, 3][None, :]) @ r_mol) / s

    def transform_instance_matrices_to_world(self, local_matrices: np.ndarray) -> np.ndarray:
        """Transform renderer (N, 4, 4) local instance matrices into world (N, 4, 4) matrices."""
        mats = np.asarray(local_matrices, dtype=np.float64)
        if mats.ndim != 3 or mats.shape[1:] != (4, 4):
            raise ValueError(f"local_matrices must have shape (N, 4, 4), got {mats.shape}")
        comp = self.T_world_from_local
        return np.einsum("ij,njk->nik", comp, mats)

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell_id": int(self.cell_id),
            "local_unit": self.local_unit,
            "unit_scale_m": float(self.unit_scale_m),
            "T_cell": np.asarray(self.T_cell, dtype=np.float64).tolist(),
            "T_molecule": np.asarray(self.T_molecule, dtype=np.float64).tolist(),
            "T_world_from_local": self.T_world_from_local.tolist(),
        }


@dataclass(frozen=True)
class M7ScaleProvenanceBadge:
    """Epistemic disclosure badge for a property at a specific microscope zoom scale."""

    scale: str
    layer_name: str
    disclosure_class: str
    evidence_class: str
    producer: str
    source_repo_or_db: str
    source_revision_or_accession: str
    notes: str = ""

    def __post_init__(self) -> None:
        if self.scale not in VALID_MICROSCOPE_SCALES:
            raise ValueError(f"scale must be one of {VALID_MICROSCOPE_SCALES}")
        if self.disclosure_class not in VALID_DISCLOSURE_CLASSES:
            raise ValueError(f"disclosure_class must be one of {VALID_DISCLOSURE_CLASSES}")
        if self.evidence_class not in VALID_EVIDENCE_CLASSES:
            raise ValueError(f"evidence_class must be one of {VALID_EVIDENCE_CLASSES}")
        if not self.layer_name.strip():
            raise ValueError("layer_name must be non-empty")
        if not self.source_repo_or_db.strip() or not self.source_revision_or_accession.strip():
            raise ValueError("source_repo_or_db and source_revision_or_accession must be non-empty (fail-closed pinning)")

    def to_dict(self) -> dict[str, str]:
        return {
            "scale": self.scale,
            "layer_name": self.layer_name,
            "disclosure_class": self.disclosure_class,
            "evidence_class": self.evidence_class,
            "producer": self.producer,
            "source_repo_or_db": self.source_repo_or_db,
            "source_revision_or_accession": self.source_revision_or_accession,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class M7MolecularAnchorBundle:
    """Coupled Sobolev molecular scene anchored to a continuum cell in M7SpatialFrame."""

    anchor_id: str
    nested_transform: M7NestedTransform
    smacro1_payload: bytes
    smacro1_sha256: str
    frame_index: int
    sobolev_relaxation_step: int
    local_coords: np.ndarray
    world_coords_m: np.ndarray
    asset_ids: np.ndarray
    prototypes: tuple[sv.AssetPrototype, ...]
    coarse_grain_result: sv.CoarseGrainResult
    coarse_world_coords_m: np.ndarray
    world_instance_matrices: np.ndarray
    disclosure_class: str = "coarse_grained_approximation"
    evidence_class: str = "synthetic"
    source_commit_sha: str = "a6e4917cdba8354421cd8ec2842a8047b698b448"

    def to_manifest_dict(self, binary_file_name: str = "sobolev_frame.smacro1") -> dict[str, Any]:
        verify_sha256_prefixed(self.smacro1_sha256, "smacro1_sha256")
        return {
            "anchor_id": self.anchor_id,
            "cell_id": int(self.nested_transform.cell_id),
            "binary_file": binary_file_name,
            "smacro1_sha256": self.smacro1_sha256,
            "frame_index": int(self.frame_index),
            "sobolev_relaxation_step": int(self.sobolev_relaxation_step),
            "n_nodes_full": int(self.local_coords.shape[0]),
            "n_nodes_lod": int(self.coarse_grain_result.coords.shape[0]),
            "n_coarse_groups": int(np.sum(self.coarse_grain_result.is_coarse)),
            "source_indices": [list(g) for g in self.coarse_grain_result.source_indices],
            "is_coarse": [bool(x) for x in self.coarse_grain_result.is_coarse],
            "nested_transform": self.nested_transform.to_dict(),
            "prototypes": [p.to_dict() for p in self.prototypes],
            "disclosure_class": self.disclosure_class,
            "evidence_class": self.evidence_class,
            "source_commit_sha": self.source_commit_sha,
        }


class M7SobolevSceneAdapter:
    """Adapter bridging SobolevMacromolecule coordinate frames to M7SpatialFrame v1.0.0."""

    def __init__(
        self,
        *,
        source_repo: str = "https://github.com/aurascoper/SobolevMacromolecule",
        source_commit_sha: str = "a6e4917cdba8354421cd8ec2842a8047b698b448",
    ) -> None:
        if not source_repo.strip():
            raise ValueError("source_repo must be non-empty")
        sha_clean = source_commit_sha.strip().lower()
        if len(sha_clean) != 40 or any(c not in "0123456789abcdef" for c in sha_clean):
            raise ValueError("source_commit_sha must be an immutable 40-character lowercase hex commit SHA")
        self.source_repo = source_repo.strip()
        self.source_commit_sha = sha_clean

    def create_molecular_anchor(
        self,
        anchor_id: str,
        coords_local: np.ndarray | Sequence[Sequence[float]],
        asset_ids: Sequence[int] | np.ndarray,
        prototypes: Sequence[sv.AssetPrototype],
        nested_transform: M7NestedTransform,
        *,
        frame_index: int = 0,
        sobolev_relaxation_step: int = 0,
        active_center_local: Sequence[float] | np.ndarray = (0.0, 0.0, 0.0),
        active_radius_local: float = 25.0,
        coarse_group_size: int = 5,
        local_scales: float | Sequence[float] | np.ndarray = 1.0,
        local_rotations: np.ndarray | None = None,
        disclosure_class: str = "coarse_grained_approximation",
        evidence_class: str = "synthetic",
    ) -> M7MolecularAnchorBundle:
        """Encode a Sobolev molecular coordinate set into SMACRO1 and place it in M7 world space."""
        if not anchor_id or not anchor_id.strip():
            raise ValueError("anchor_id must be non-empty")
        if int(sobolev_relaxation_step) < 0:
            raise ValueError("sobolev_relaxation_step must be non-negative")
        if disclosure_class not in VALID_DISCLOSURE_CLASSES:
            raise ValueError(f"disclosure_class must be one of {VALID_DISCLOSURE_CLASSES}")
        if evidence_class not in VALID_EVIDENCE_CLASSES:
            raise ValueError(f"evidence_class must be one of {VALID_EVIDENCE_CLASSES}")

        arr_local = sv.validate_coords(coords_local, name="coords_local").astype(np.float32, copy=True)
        ids_arr = np.asarray(asset_ids, dtype=np.uint32)
        if ids_arr.shape != (arr_local.shape[0],):
            raise ValueError(f"asset_ids shape {ids_arr.shape} does not match {(arr_local.shape[0],)}")
        if not prototypes:
            raise ValueError("at least one AssetPrototype is required")
        proto_tuple = tuple(prototypes)
        proto_ids = set()
        for proto in proto_tuple:
            proto.validate()
            proto_ids.add(int(proto.asset_id))
        missing = sorted(set(int(x) for x in ids_arr) - proto_ids)
        if missing:
            raise ValueError(f"missing AssetPrototype definitions for asset_ids {missing}")

        # Note: in SMACRO1, time_seconds records the relaxation step index as float metadata
        # while M7MultirateClockState keeps sobolev_relaxation_step strictly integer-typed.
        payload = sv.encode_coordinate_frame(
            arr_local,
            frame_index=int(frame_index),
            time_seconds=float(sobolev_relaxation_step),
            asset_ids=ids_arr,
        )
        payload_hash = sha256_prefixed(payload)

        cg = sv.coarse_grain_far_field(
            arr_local,
            active_center=active_center_local,
            active_radius=float(active_radius_local),
            group_size=int(coarse_group_size),
        )

        world_coords = nested_transform.transform_local_to_world(arr_local)
        coarse_world_coords = (
            nested_transform.transform_local_to_world(cg.coords)
            if cg.coords.shape[0] > 0
            else np.empty((0, 3), dtype=np.float64)
        )
        local_mats = sv.instance_matrices(
            arr_local,
            scales=local_scales,
            rotations=local_rotations,
            dtype=np.float64,
        )
        world_mats = nested_transform.transform_instance_matrices_to_world(local_mats)

        return M7MolecularAnchorBundle(
            anchor_id=anchor_id.strip(),
            nested_transform=nested_transform,
            smacro1_payload=payload,
            smacro1_sha256=payload_hash,
            frame_index=int(frame_index),
            sobolev_relaxation_step=int(sobolev_relaxation_step),
            local_coords=arr_local,
            world_coords_m=world_coords,
            asset_ids=ids_arr,
            prototypes=proto_tuple,
            coarse_grain_result=cg,
            coarse_world_coords_m=coarse_world_coords,
            world_instance_matrices=world_mats,
            disclosure_class=disclosure_class,
            evidence_class=evidence_class,
            source_commit_sha=self.source_commit_sha,
        )

    def decode_and_verify_anchor(
        self,
        smacro1_payload: bytes | bytearray | memoryview,
        expected_sha256: str,
        nested_transform: M7NestedTransform,
        *,
        active_center_local: Sequence[float] | np.ndarray = (0.0, 0.0, 0.0),
        active_radius_local: float = 25.0,
        coarse_group_size: int = 5,
    ) -> dict[str, Any]:
        """Verify SMACRO1 digest, decode coordinates, and project to M7 world space."""
        verify_sha256_prefixed(expected_sha256, "expected_sha256")
        actual_sha256 = sha256_prefixed(bytes(smacro1_payload))
        if actual_sha256 != expected_sha256:
            raise ValueError(f"SMACRO1 digest mismatch: expected {expected_sha256}, got {actual_sha256}")

        decoded = sv.decode_coordinate_frame(smacro1_payload)
        world_coords = nested_transform.transform_local_to_world(decoded.coords)
        roundtrip_local = nested_transform.transform_world_to_local(world_coords)
        max_roundtrip_err = float(np.max(np.abs(roundtrip_local - decoded.coords))) if decoded.n_nodes > 0 else 0.0

        cg = sv.coarse_grain_far_field(
            decoded.coords,
            active_center=active_center_local,
            active_radius=float(active_radius_local),
            group_size=int(coarse_group_size),
        )
        coarse_world = (
            nested_transform.transform_local_to_world(cg.coords)
            if cg.coords.shape[0] > 0
            else np.empty((0, 3), dtype=np.float64)
        )

        return {
            "frame": decoded,
            "smacro1_sha256": actual_sha256,
            "world_coords_m": world_coords,
            "coarse_grain": cg,
            "coarse_world_coords_m": coarse_world,
            "max_roundtrip_local_error": max_roundtrip_err,
        }

    def export_m7_usda_point_instancer(
        self,
        path: str | Path,
        anchor: M7MolecularAnchorBundle,
        clock: M7MultirateClockState,
        *,
        simulation_id: str = "m7-reference-case",
        use_world_coordinates: bool = True,
        scene_name: str = "M7SobolevScene",
    ) -> str:
        """Write an OpenUSD ASCII PointInstancer scene with M7 provenance metadata."""
        coords_to_write = anchor.world_coords_m if use_world_coordinates else anchor.local_coords
        temp_path = Path(path)
        sv.write_usda_point_instancer(
            temp_path,
            coords_to_write,
            anchor.asset_ids,
            anchor.prototypes,
            scene_name=scene_name,
        )
        base_usda = temp_path.read_text(encoding="utf-8")
        coord_space = "world_m" if use_world_coordinates else anchor.nested_transform.local_unit
        meta_block = "\n".join(
            [
                "#usda 1.0",
                "(",
                f'    defaultPrim = "{sv._usd_identifier(scene_name)}"',
                '    doc = "M7MultiscaleExplorer OpenUSD PointInstancer Handoff"',
                "    customLayerData = {",
                f'        string m7_schema = "{M7_SCHEMA_NAME}"',
                f'        string m7_schema_version = "{M7_SCHEMA_VERSION}"',
                f'        string simulation_id = "{simulation_id}"',
                f'        string anchor_id = "{anchor.anchor_id}"',
                f"        int cell_id = {int(anchor.nested_transform.cell_id)}",
                f"        double physical_time_s = {float(clock.physical_time_s)}",
                f"        double biological_time_s = {float(clock.biological_time_s)}",
                f"        int sobolev_relaxation_step = {int(clock.sobolev_relaxation_step)}",
                f'        string snapshot_sync_id = "{clock.snapshot_sync_id}"',
                f'        string coordinate_space = "{coord_space}"',
                f'        string smacro1_sha256 = "{anchor.smacro1_sha256}"',
                f'        string disclosure_class = "{anchor.disclosure_class}"',
                f'        string evidence_class = "{anchor.evidence_class}"',
                f'        string source_commit_sha = "{self.source_commit_sha}"',
                "    }",
                ")",
            ]
        )
        # Replace the simple 4-line header generated by sv.write_usda_point_instancer
        header_end = base_usda.find("\n)\n")
        if header_end != -1:
            updated_usda = meta_block + base_usda[header_end + 3 :]
        else:
            updated_usda = meta_block + "\n" + base_usda
        temp_path.write_text(updated_usda, encoding="utf-8")
        return sha256_prefixed(updated_usda)


def build_m7_spatial_frame_dict(
    *,
    simulation_id: str,
    frame_index: int,
    clock: M7MultirateClockState,
    cell_count: int,
    geometry_asset: str,
    geometry_sha256: str,
    fields: Sequence[Mapping[str, Any]],
    molecular_anchors: Sequence[M7MolecularAnchorBundle] = (),
    scale_disclosures: Sequence[M7ScaleProvenanceBadge] = (),
    topology: str = "unstructured",
    cell_ordering: str = "explicit_id_ascending",
    coordinate_origin_m: Sequence[float] = (0.0, 0.0, 0.0),
    evidence_class: str = "synthetic",
    solver_versions: Mapping[str, str] | None = None,
    source_hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Build and validate a canonical M7SpatialFrame v1.0.0 dictionary."""
    if not simulation_id or not simulation_id.strip():
        raise ValueError("simulation_id must be non-empty")
    if int(frame_index) < 0:
        raise ValueError("frame_index must be non-negative")
    if int(cell_count) <= 0:
        raise ValueError("cell_count must be positive")
    if topology not in VALID_TOPOLOGIES:
        raise ValueError(f"topology must be one of {VALID_TOPOLOGIES}")
    if evidence_class not in VALID_EVIDENCE_CLASSES:
        raise ValueError(f"evidence_class must be one of {VALID_EVIDENCE_CLASSES}")
    verify_sha256_prefixed(geometry_sha256, "geometry_sha256")

    origin = [float(x) for x in coordinate_origin_m]
    if len(origin) != 3 or not all(np.isfinite(origin)):
        raise ValueError("coordinate_origin_m must be a finite 3-vector")

    normalized_fields: list[dict[str, Any]] = []
    seen_field_names: set[str] = set()
    for f in fields:
        fname = str(f.get("name", "")).strip()
        if not fname:
            raise ValueError("field name must be non-empty")
        if fname in seen_field_names:
            raise ValueError(f"duplicate field name '{fname}'")
        seen_field_names.add(fname)

        producer = str(f.get("producer", "")).strip()
        if producer not in VALID_PRODUCERS:
            raise ValueError(f"field '{fname}' producer '{producer}' must be one of {VALID_PRODUCERS}")
        assoc = str(f.get("association", "cell")).strip()
        if assoc not in VALID_ASSOCIATIONS:
            raise ValueError(f"field '{fname}' association '{assoc}' must be one of {VALID_ASSOCIATIONS}")
        comps = int(f.get("components", 1))
        if comps not in (1, 3):
            raise ValueError(f"field '{fname}' components must be 1 or 3, got {comps}")
        unit = str(f.get("unit", "")).strip()
        if not unit:
            raise ValueError(f"field '{fname}' unit must be non-empty")
        f_sha = str(f.get("sha256", "")).strip()
        verify_sha256_prefixed(f_sha, f"field '{fname}' sha256")

        validity_mask = f.get("validity_mask")
        if validity_mask is not None:
            if len(validity_mask) != int(cell_count):
                raise ValueError(f"field '{fname}' validity_mask length {len(validity_mask)} != cell_count {cell_count}")
            validity_mask = [bool(v) for v in validity_mask]

        uncertainty = f.get("uncertainty")
        if uncertainty is not None:
            if len(uncertainty) != int(cell_count):
                raise ValueError(f"field '{fname}' uncertainty length {len(uncertainty)} != cell_count {cell_count}")
            uncertainty = [float(u) for u in uncertainty]
            if any((not np.isfinite(u)) or u < 0.0 for u in uncertainty):
                raise ValueError(f"field '{fname}' uncertainty values must be finite and non-negative")

        normalized_fields.append(
            {
                "name": fname,
                "producer": producer,
                "association": assoc,
                "components": comps,
                "unit": unit,
                "file": str(f.get("file", f"{fname}.bin")),
                "sha256": f_sha,
                "validity_mask": validity_mask,
                "uncertainty": uncertainty,
            }
        )

    anchor_dicts = [
        a.to_manifest_dict(binary_file_name=f"{a.anchor_id}.smacro1")
        for a in molecular_anchors
    ]
    badge_dicts = [b.to_dict() for b in scale_disclosures]

    s_hashes = dict(source_hashes or {})
    for h_key, h_val in s_hashes.items():
        verify_sha256_prefixed(h_val, f"source_hashes['{h_key}']")

    frame_dict: dict[str, Any] = {
        "schema": M7_SCHEMA_NAME,
        "schema_version": M7_SCHEMA_VERSION,
        "simulation_id": simulation_id.strip(),
        "frame_index": int(frame_index),
        "time": {
            "value": float(clock.physical_time_s),
            "unit": "s",
        },
        "clocks": clock.to_dict(),
        "coordinate_frame": {
            "name": "cartesian_right_handed",
            "unit": "m",
            "origin": origin,
        },
        "mesh": {
            "topology": topology,
            "geometry_asset": geometry_asset,
            "geometry_sha256": geometry_sha256,
            "cell_count": int(cell_count),
            "cell_ordering": cell_ordering,
        },
        "fields": normalized_fields,
        "molecular_anchors": anchor_dicts,
        "microscope_disclosures": badge_dicts,
        "provenance": {
            "evidence_class": evidence_class,
            "solver_versions": dict(solver_versions or {}),
            "source_hashes": s_hashes,
        },
    }
    canonical_bytes = json.dumps(frame_dict, sort_keys=True, separators=(",", ":")).encode("utf-8")
    frame_dict["frame_sha256"] = sha256_prefixed(canonical_bytes)
    return frame_dict


def validate_m7_spatial_frame_dict(
    frame_dict: dict[str, Any],
    *,
    verify_canonical_frame_sha256: bool = True,
) -> bool:
    """Validate an M7SpatialFrame v1.0.0 dictionary against schema, 3-clock, and SHA-256 rules."""

    if not isinstance(frame_dict, dict):
        raise ValueError("frame_dict must be a dictionary")
    if frame_dict.get("schema") != M7_SCHEMA_NAME:
        raise ValueError(f"schema must be '{M7_SCHEMA_NAME}'")
    if frame_dict.get("schema_version") != M7_SCHEMA_VERSION:
        raise ValueError(f"schema_version must be '{M7_SCHEMA_VERSION}'")
    if not str(frame_dict.get("simulation_id", "")).strip():
        raise ValueError("simulation_id must be non-empty")
    if int(frame_dict.get("frame_index", -1)) < 0:
        raise ValueError("frame_index must be >= 0")

    clocks = frame_dict.get("clocks")
    if not isinstance(clocks, dict):
        raise ValueError("clocks section is required")
    M7MultirateClockState(
        physical_time_s=clocks.get("physical_time_s", -1.0),
        biological_time_s=clocks.get("biological_time_s", -1.0),
        sobolev_relaxation_step=clocks.get("sobolev_relaxation_step", -1),
        snapshot_sync_id=str(clocks.get("snapshot_sync_id", "")),
        interpolation_valid=bool(clocks.get("interpolation_valid", True)),
        parent_snapshot_id=clocks.get("parent_snapshot_id"),
    )

    mesh = frame_dict.get("mesh")
    if not isinstance(mesh, dict):
        raise ValueError("mesh section is required")
    if mesh.get("topology") not in VALID_TOPOLOGIES:
        raise ValueError(f"mesh.topology must be one of {VALID_TOPOLOGIES}")
    cell_count = int(mesh.get("cell_count", 0))
    if cell_count < 1:
        raise ValueError("mesh.cell_count must be >= 1")
    verify_sha256_prefixed(str(mesh.get("geometry_sha256", "")), "mesh.geometry_sha256")

    fields = frame_dict.get("fields")
    if not isinstance(fields, list) or len(fields) == 0:
        raise ValueError("fields must be a non-empty list")
    seen_names: set[str] = set()
    for f in fields:
        fname = str(f.get("name", "")).strip()
        if not fname or fname in seen_names:
            raise ValueError(f"invalid or duplicate field name '{fname}'")
        seen_names.add(fname)
        if f.get("producer") not in VALID_PRODUCERS:
            raise ValueError(f"field '{fname}' has invalid producer '{f.get('producer')}'")
        if f.get("association") not in VALID_ASSOCIATIONS:
            raise ValueError(f"field '{fname}' has invalid association '{f.get('association')}'")
        verify_sha256_prefixed(str(f.get("sha256", "")), f"field '{fname}' sha256")
        vmask = f.get("validity_mask")
        if vmask is not None and len(vmask) != cell_count:
            raise ValueError(f"field '{fname}' validity_mask length != cell_count")

    for anchor in frame_dict.get("molecular_anchors", []):
        cid = int(anchor.get("cell_id", 0))
        if cid < 1 or cid > cell_count:
            raise ValueError(f"molecular anchor cell_id {cid} outside 1..{cell_count}")
        verify_sha256_prefixed(str(anchor.get("smacro1_sha256", "")), "anchor.smacro1_sha256")

    prov = frame_dict.get("provenance")
    if not isinstance(prov, dict) or prov.get("evidence_class") not in VALID_EVIDENCE_CLASSES:
        raise ValueError("provenance.evidence_class must be valid")

    f_sha = str(frame_dict.get("frame_sha256", ""))
    verify_sha256_prefixed(f_sha, "frame_sha256")

    if verify_canonical_frame_sha256:
        copy_without_sha = {k: v for k, v in frame_dict.items() if k != "frame_sha256"}
        expected_sha = sha256_prefixed(
            json.dumps(copy_without_sha, sort_keys=True, separators=(",", ":")).encode("utf-8")
        )
        if f_sha != expected_sha:
            raise ValueError(
                f"frame_sha256 mismatch (tampered payload): recorded {f_sha} != recomputed {expected_sha}"
            )

    return True


__all__ = [
    "M7_ADAPTER_VERSION",
    "M7_SCHEMA_NAME",
    "M7_SCHEMA_VERSION",
    "M7MolecularAnchorBundle",
    "M7MultirateClockState",
    "M7NestedTransform",
    "M7ScaleProvenanceBadge",
    "M7SobolevSceneAdapter",
    "UNIT_TO_METERS",
    "VALID_ASSOCIATIONS",
    "VALID_DISCLOSURE_CLASSES",
    "VALID_EVIDENCE_CLASSES",
    "VALID_MICROSCOPE_SCALES",
    "VALID_PRODUCERS",
    "VALID_TOPOLOGIES",
    "build_m7_spatial_frame_dict",
    "sha256_prefixed",
    "validate_m7_spatial_frame_dict",
    "validate_temporal_interpolation",
    "verify_sha256_prefixed",
]

