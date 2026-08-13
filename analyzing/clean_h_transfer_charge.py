#!/usr/bin/env python3
"""
Estimate the effective charge carried in clean H-transfer events.

This script is designed to run after the existing H-transfer detector.  It
reads h_rearrangement_events.csv, finds the nearest pre- and post-transfer
frames at which the transferred H is covalently associated with the expected
O, and reports both (1) the charge estimate from those two nearest frames and
(2) the estimate from four-frame averages immediately outside the anchors.
The two estimates are written together so that their sampling sensitivity can
be assessed directly.

The density partition is the same nearest-nucleus ("Cody") real-space method
used in density_integrate.py: every grid point within MAX_DENSITY_DISTANCE_A of
an atom is assigned to the nearest atom.  Atomic charges are Z - N_e.

Fragment bookkeeping uses fixed atom identities.  At the pre-transfer frame,
each non-transferred H assigned to the donor or acceptor is either:
  * excluded from both frames if it had already met the H-ejection definition;
  * included in the same parent fragment at both frames otherwise, even if it
    ejects between the two frames.
The transferred H is the only atom moved from donor to acceptor.  Thus H
ejection does not create a trivial change in fragment nuclear composition.

For corrected donor/acceptor fragment charges Q_D and Q_A,

    q_H(donor)    = Q_D(before) - Q_D(after)
    q_H(acceptor) = Q_A(after)  - Q_A(before)

Agreement between these two values is the main internal sanity check.

Run from the directory that contains the tdl* folders and the event CSV:

    python clean_h_transfer_charge_analysis.py

Only NumPy and the Python standard library are required.
"""

import csv
import faulthandler
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# If a compiled numerical library or the operating system terminates Python,
# print the active Python stack instead of only the words "Segmentation fault".
faulthandler.enable(all_threads=True)

import numpy as np


# ============================================================================
# USER-TUNABLE PARAMETERS
# ============================================================================

# Input layout
BASE_DIR = Path(".")
EVENTS_CSV_NAME = "h_rearrangement_events.csv"
XYZ_NAME = "trajectory.xyz"
DENSITY_SUBDIR = ""                 # e.g. "density" if dens*.bov is elsewhere
BOV_PREFIX = "dens"
BOV_DIGITS = 5

# Which time reported by the previous script should define the transfer time.
# Usually start_time_fs is the first frame captured by the new O.  Alternatives
# in the existing CSV are stable_time_fs and confirm_time_fs.
TRANSFER_TIME_COLUMN = "start_time_fs"

# Geometry criterion for the transferred O-H only (Angstrom).  Other O-H bond
# lengths do not affect anchor selection.
TRANSFERRED_OH_BOND_CUTOFF_A = 1.5

# H-ejection definition, matched to h_ejection_by_run.py.
# Clean: first frame with nearest-O distance > 3.0 A before either relevant
# atom reaches the CAP.  CAP special case: at the first frame at which H or its
# nearest O reaches the CAP, count ejection only if nearest-O distance > 2.5 A.
CAP_BOUND_A = 11.0
H_EJECTION_CLEAN_DISTANCE_A = 3.0
H_EJECTION_CAP_DISTANCE_A = 2.5

# A target H that was already classified as ejected and was later recaptured
# is not a clean direct H-transfer event.  Keep this True for the main analysis.
REJECT_IF_TRANSFERRED_H_WAS_EJECTED = True

# Density partitioning close to the CAP is difficult to interpret.  Reject an
# anchor if the transferred H, donor O, or acceptor O is already in the CAP.
REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS = True

# Both single-frame charge samples must be at/after this time, when the
# principal ionization stage is assumed to be over.
MIN_ANALYSIS_TIME_FS = 17.5

# Search no farther than this from the transfer time for a clean anchor.  Set
# to None to search all available trajectory frames.
MAX_ANCHOR_SEARCH_FS = None

# A clean event may contain H ejection, but no second H-transfer event may
# begin from the selected pre-transfer frame through the selected
# post-transfer frame.
REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW = True

# Time matching.  Charge sample times should normally coincide with saved XYZ
# and density frames.  Increase only if output times have small rounding noise.
XYZ_TIME_MATCH_TOL_FS = 1.0e-4
DENSITY_TIME_MATCH_TOL_FS = 1.0e-6
EVENT_TIME_COMPARE_TOL_FS = 1.0e-6

# Density-file timing, copied from density_integrate(2).py.  This controls the
# mapping time -> densNNNNN.bov.
SIMULATION_TIME_STEP_FS = 0.001
TD_DENSITY_OUTPUT_FREQUENCY = 500
DENSITY_OUTPUT_DT_FS = SIMULATION_TIME_STEP_FS * TD_DENSITY_OUTPUT_FREQUENCY

# Four-frame comparison.  With the current 0.5 fs density interval, the before
# samples are anchor-2.0, -1.5, -1.0, and -0.5 fs, while the after samples are
# anchor+0.5, +1.0, +1.5, and +2.0 fs.  The anchors themselves are reported
# separately as the nearest-frame result.
FOUR_FRAME_COUNT = 4
FOUR_FRAME_SPACING_FS = DENSITY_OUTPUT_DT_FS
REQUIRE_TARGET_OH_GEOMETRY_AT_FOUR_FRAME_SAMPLES = True

# Nearest-nucleus density partition (same physical method as the uploaded
# density integration script).
MAX_DENSITY_DISTANCE_A = 10.0
DENSITY_BINARY_DTYPE = "<f4"         # little-endian float32

# Hard upper bound on the number of grid points handled at once.  The original
# vectorized implementation formed a (points x atoms x 3) temporary array,
# which could exhaust memory and be reported by some HPC systems as a bare
# segmentation fault.  The streaming implementation below loops over atoms and
# therefore uses O(points), rather than O(points * atoms), temporary memory.
DENSITY_POINT_CHUNK_SIZE = 100_000

# Print each density file immediately before it is integrated.  Keeping this on
# is useful on batch systems, whose buffered stdout can otherwise make a later
# failure look as if it happened at program startup.
PRINT_DENSITY_PROGRESS = True

# Nuclear charges needed for water.  Add entries here for other elements.
NUCLEAR_CHARGE = {"H": 1.0, "O": 8.0}

# Sanity-check tolerances (elementary charge, e)
DONOR_ACCEPTOR_MISMATCH_TOL_E = 1.25
SYSTEM_CHARGE_CHANGE_TOL_E = 1.25

# Outputs
OUT_EVENTS_CSV = "clean_h_transfer_charges.csv"
OUT_FRAMES_CSV = "clean_h_transfer_charge_frames.csv"
OUT_SUMMARY_CSV = "clean_h_transfer_charge_summary.csv"
OUT_EJECTIONS_CSV = "clean_h_transfer_detected_ejections.csv"


TIME_RE = re.compile(r"time\[fs\]\s*=\s*([0-9Ee+\-.]+)")
ITER_RE = re.compile(r"iter\s*=\s*([0-9]+)")


# ============================================================================
# XYZ trajectory and geometry analysis
# ============================================================================

def parse_xyz_time(comment, iter_index):
    match = TIME_RE.search(comment)
    if match:
        return float(match.group(1))
    if iter_index is not None:
        return iter_index * SIMULATION_TIME_STEP_FS
    raise ValueError("XYZ frame has neither time[fs] nor iter")


def read_xyz_trajectory(path):
    """Read trajectory.xyz while preserving the fixed atom ordering."""
    times = []
    iterations = []
    all_coords = []
    reference_species = None

    with path.open("r") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.strip():
                continue
            natoms = int(line.strip())
            comment = handle.readline()
            if not comment:
                raise RuntimeError("EOF after XYZ atom count in %s" % path)

            iter_match = ITER_RE.search(comment)
            iter_index = int(iter_match.group(1)) if iter_match else None
            time_fs = parse_xyz_time(comment, iter_index)

            species = []
            coords = np.empty((natoms, 3), dtype=float)
            for atom_i in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise RuntimeError("EOF inside XYZ frame in %s" % path)
                parts = atom_line.split()
                if len(parts) < 4:
                    raise RuntimeError("Bad XYZ atom line: %s" % atom_line.strip())
                species.append(parts[0].capitalize())
                coords[atom_i] = [float(parts[1]), float(parts[2]), float(parts[3])]

            if reference_species is None:
                reference_species = species
            elif species != reference_species:
                raise ValueError("Atom identities/order change between frames in %s" % path)

            times.append(time_fs)
            iterations.append(iter_index)
            all_coords.append(coords)

    if not all_coords:
        raise ValueError("No XYZ frames found in %s" % path)

    times_array = np.asarray(times, dtype=float)
    if np.any(np.diff(times_array) <= 0.0):
        raise ValueError("XYZ times are not strictly increasing in %s" % path)

    return Trajectory(
        path=path,
        times=times_array,
        iterations=iterations,
        species=reference_species,
        coords=np.stack(all_coords, axis=0),
    )


class Trajectory:
    def __init__(self, path, times, iterations, species, coords):
        self.path = path
        self.times = times
        self.iterations = iterations
        self.species = species
        self.coords = coords
        self.oxygen_indices = np.asarray(
            [i for i, symbol in enumerate(species) if symbol == "O"], dtype=int
        )
        self.hydrogen_indices = np.asarray(
            [i for i, symbol in enumerate(species) if symbol == "H"], dtype=int
        )
        if self.oxygen_indices.size == 0 or self.hydrogen_indices.size == 0:
            raise ValueError("Trajectory must contain both O and H atoms: %s" % path)

        self.nearest_o = np.full(
            (len(times), len(species)), -1, dtype=int
        )
        self.nearest_o_distance = np.full(
            (len(times), len(species)), np.inf, dtype=float
        )
        self._calculate_nearest_oxygen()
        self.ejection_time = {}
        self.ejection_frame_index = {}
        self.ejection_type = {}
        self.ejection_nearest_o = {}
        self.ejection_distance = {}
        self.ejection_h_in_cap = {}
        self.ejection_o_in_cap = {}
        self.ejection_tracking_stop_time = {}
        self.ejection_tracking_stop_reason = {}
        self._detect_h_ejections()

    def _calculate_nearest_oxygen(self):
        oxygen_coords = self.coords[:, self.oxygen_indices, :]
        for h_index in self.hydrogen_indices:
            displacement = oxygen_coords - self.coords[:, h_index, :][:, None, :]
            distances = np.sqrt(np.sum(displacement * displacement, axis=2))
            nearest_local = np.argmin(distances, axis=1)
            rows = np.arange(len(self.times))
            self.nearest_o[:, h_index] = self.oxygen_indices[nearest_local]
            self.nearest_o_distance[:, h_index] = distances[rows, nearest_local]

    def _detect_h_ejections(self):
        """Apply the exact first-crossing/CAP logic of h_ejection_by_run.py."""
        for h_index in self.hydrogen_indices:
            h_index = int(h_index)
            self.ejection_time[h_index] = None
            self.ejection_frame_index[h_index] = None
            self.ejection_type[h_index] = None
            self.ejection_nearest_o[h_index] = None
            self.ejection_distance[h_index] = None
            self.ejection_h_in_cap[h_index] = None
            self.ejection_o_in_cap[h_index] = None
            self.ejection_tracking_stop_time[h_index] = None
            self.ejection_tracking_stop_reason[h_index] = None

            for frame_i in range(len(self.times)):
                nearest_o = int(self.nearest_o[frame_i, h_index])
                nearest_distance = float(self.nearest_o_distance[frame_i, h_index])
                h_in_cap = self.atom_is_in_cap(frame_i, h_index)
                o_in_cap = self.atom_is_in_cap(frame_i, nearest_o)

                # The reference script stops tracking at the first CAP touch,
                # whether or not the 2.5 A fallback criterion is satisfied.
                if h_in_cap or o_in_cap:
                    self.ejection_tracking_stop_time[h_index] = float(self.times[frame_i])
                    if nearest_distance > H_EJECTION_CAP_DISTANCE_A:
                        self._record_ejection(
                            h_index, frame_i, "CAP", nearest_o,
                            nearest_distance, h_in_cap, o_in_cap,
                        )
                        self.ejection_tracking_stop_reason[h_index] = "CAP_EJECTION"
                    else:
                        self.ejection_tracking_stop_reason[h_index] = "CAP_NO_EJECTION"
                    break

                if nearest_distance > H_EJECTION_CLEAN_DISTANCE_A:
                    self._record_ejection(
                        h_index, frame_i, "CLEAN", nearest_o,
                        nearest_distance, h_in_cap, o_in_cap,
                    )
                    self.ejection_tracking_stop_time[h_index] = float(self.times[frame_i])
                    self.ejection_tracking_stop_reason[h_index] = "CLEAN_EJECTION"
                    break

    def _record_ejection(
        self, h_index, frame_i, event_type, nearest_o, nearest_distance,
        h_in_cap, o_in_cap,
    ):
        self.ejection_time[h_index] = float(self.times[frame_i])
        self.ejection_frame_index[h_index] = int(frame_i)
        self.ejection_type[h_index] = event_type
        self.ejection_nearest_o[h_index] = int(nearest_o)
        self.ejection_distance[h_index] = float(nearest_distance)
        self.ejection_h_in_cap[h_index] = bool(h_in_cap)
        self.ejection_o_in_cap[h_index] = bool(o_in_cap)

    def atom_is_in_cap(self, frame_i, atom_i):
        x, y, z = self.coords[frame_i, atom_i]
        return (
            abs(float(x)) >= CAP_BOUND_A
            or abs(float(y)) >= CAP_BOUND_A
            or abs(float(z)) >= CAP_BOUND_A
        )

    def closest_frame_index(self, target_time_fs, tolerance_fs=XYZ_TIME_MATCH_TOL_FS):
        insertion = int(np.searchsorted(self.times, target_time_fs))
        candidates = []
        if insertion < len(self.times):
            candidates.append(insertion)
        if insertion > 0:
            candidates.append(insertion - 1)
        if not candidates:
            raise ValueError("No XYZ frames available")
        best = min(candidates, key=lambda i: abs(self.times[i] - target_time_fs))
        mismatch = abs(self.times[best] - target_time_fs)
        if mismatch > tolerance_fs:
            raise ValueError(
                "No XYZ frame at %.9f fs; closest is %.9f fs (difference %.9g fs)"
                % (target_time_fs, self.times[best], mismatch)
            )
        return best

    def h_is_ejected_by(self, h_index, time_fs):
        eject_time = self.ejection_time.get(h_index)
        return eject_time is not None and eject_time <= time_fs + EVENT_TIME_COMPARE_TOL_FS

    def geometry_is_clean_for_side(self, frame_i, target_h, target_o):
        if self.times[frame_i] + EVENT_TIME_COMPARE_TOL_FS < MIN_ANALYSIS_TIME_FS:
            return False
        if REJECT_IF_TRANSFERRED_H_WAS_EJECTED and self.h_is_ejected_by(
            target_h, self.times[frame_i]
        ):
            return False
        if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
            if self.atom_is_in_cap(frame_i, target_h) or self.atom_is_in_cap(frame_i, target_o):
                return False
        return (
            int(self.nearest_o[frame_i, target_h]) == target_o
            and self.nearest_o_distance[frame_i, target_h]
            <= TRANSFERRED_OH_BOND_CUTOFF_A
        )

    def find_anchor(self, transfer_time_fs, target_h, target_o, direction):
        if direction == "before":
            indices = np.where(self.times <= transfer_time_fs + EVENT_TIME_COMPARE_TOL_FS)[0][::-1]
        elif direction == "after":
            indices = np.where(self.times >= transfer_time_fs - EVENT_TIME_COMPARE_TOL_FS)[0]
        else:
            raise ValueError("direction must be before or after")

        for frame_i in indices:
            if MAX_ANCHOR_SEARCH_FS is not None:
                if abs(self.times[frame_i] - transfer_time_fs) > MAX_ANCHOR_SEARCH_FS:
                    continue
            if self.geometry_is_clean_for_side(frame_i, target_h, target_o):
                return int(frame_i)
        return None

    def build_fixed_fragment_ledger(
        self, before_frame_i, after_frame_i, target_h, donor_o, acceptor_o,
    ):
        """Build before/after atom lists with only target H changing parent."""
        before_time = float(self.times[before_frame_i])
        after_time = float(self.times[after_frame_i])

        donor_fixed_h = []
        acceptor_fixed_h = []
        excluded_ejected_before_donor = []
        excluded_ejected_before_acceptor = []

        for h_index_value in self.hydrogen_indices:
            h_index = int(h_index_value)
            if h_index == target_h:
                continue
            if self.h_is_ejected_by(h_index, before_time):
                ejection_parent = self.ejection_nearest_o.get(h_index)
                if ejection_parent == donor_o:
                    excluded_ejected_before_donor.append(h_index)
                elif ejection_parent == acceptor_o:
                    excluded_ejected_before_acceptor.append(h_index)
                continue

            owner_before = int(self.nearest_o[before_frame_i, h_index])
            if owner_before == donor_o:
                donor_fixed_h.append(h_index)
            elif owner_before == acceptor_o:
                acceptor_fixed_h.append(h_index)

        ejected_between_donor = [
            h for h in donor_fixed_h
            if self.ejection_time.get(h) is not None
            and before_time < self.ejection_time[h] <= after_time + EVENT_TIME_COMPARE_TOL_FS
        ]
        ejected_between_acceptor = [
            h for h in acceptor_fixed_h
            if self.ejection_time.get(h) is not None
            and before_time < self.ejection_time[h] <= after_time + EVENT_TIME_COMPARE_TOL_FS
        ]

        return {
            "donor_before_h": sorted(donor_fixed_h + [target_h]),
            "donor_after_h": sorted(donor_fixed_h),
            "acceptor_before_h": sorted(acceptor_fixed_h),
            "acceptor_after_h": sorted(acceptor_fixed_h + [target_h]),
            "excluded_ejected_before_donor": sorted(excluded_ejected_before_donor),
            "excluded_ejected_before_acceptor": sorted(excluded_ejected_before_acceptor),
            "ejected_between_donor": sorted(ejected_between_donor),
            "ejected_between_acceptor": sorted(ejected_between_acceptor),
        }


# ============================================================================
# BOV/DAT density integration
# ============================================================================

def bov_path_for_time(run_dir, target_time_fs):
    raw_index = target_time_fs / DENSITY_OUTPUT_DT_FS
    index = int(round(raw_index))
    actual_time = index * DENSITY_OUTPUT_DT_FS
    if abs(actual_time - target_time_fs) > DENSITY_TIME_MATCH_TOL_FS:
        raise ValueError(
            "Requested %.9f fs is not a density-output time (nearest %.9f fs)"
            % (target_time_fs, actual_time)
        )
    density_dir = run_dir / DENSITY_SUBDIR if DENSITY_SUBDIR else run_dir
    filename = "%s%0*d.bov" % (BOV_PREFIX, BOV_DIGITS, index)
    return density_dir / filename, actual_time, index


def read_bov_header(path):
    params = {}
    first_line = ""
    with path.open("r") as handle:
        for line_i, line in enumerate(handle):
            if line_i == 0:
                first_line = line.strip()
            if ":" in line:
                key, value = line.split(":", 1)
                params[key.strip().upper()] = value.strip()

    required = ["DATA_FILE", "DATA_SIZE", "BRICK_ORIGIN", "BRICK_SIZE"]
    missing = [key for key in required if key not in params]
    if missing:
        raise ValueError("Missing BOV fields %s in %s" % (missing, path))

    data_path = Path(params["DATA_FILE"])
    if not data_path.is_absolute():
        data_path = path.parent / data_path

    grid_size = tuple(int(v) for v in params["DATA_SIZE"].split())
    origin = tuple(float(v) for v in params["BRICK_ORIGIN"].split())
    brick_size = tuple(float(v) for v in params["BRICK_SIZE"].split())
    byte_offset = int(params.get("BYTE_OFFSET", "0"))
    return first_line, data_path, grid_size, origin, brick_size, byte_offset


def integrate_density_file_nearest_atom(
    path, grid_size, byte_offset, origin, brick_size, atom_coords,
):
    """Stream a BOV data file and apply ProcessDensityPointCodyMethod.

    BOV stores x fastest, then y, then z.  Only DENSITY_POINT_CHUNK_SIZE values
    and O(number_of_points_in_chunk) work arrays are resident at once.  This is
    deliberately more conservative than constructing a points-by-atoms array.
    """
    nx, ny, nz = grid_size
    if nx < 2 or ny < 2 or nz < 2:
        raise ValueError("Every BOV grid dimension must be at least 2")
    if DENSITY_POINT_CHUNK_SIZE < 1:
        raise ValueError("DENSITY_POINT_CHUNK_SIZE must be at least 1")

    dx = brick_size[0] / float(nx - 1)
    dy = brick_size[1] / float(ny - 1)
    dz = brick_size[2] / float(nz - 1)
    npoints = nx * ny * nz
    dtype = np.dtype(DENSITY_BINARY_DTYPE)
    required_bytes = byte_offset + dtype.itemsize * npoints
    actual_bytes = path.stat().st_size
    if actual_bytes < required_bytes:
        available_values = max(0, actual_bytes - byte_offset) // dtype.itemsize
        raise ValueError(
            "Expected %d density values in %s, found only %d"
            % (npoints, path, available_values)
        )

    atom_coords = np.asarray(atom_coords, dtype=np.float64)
    atom_density_sums = np.zeros(len(atom_coords), dtype=np.float64)
    cutoff2 = MAX_DENSITY_DISTANCE_A * MAX_DENSITY_DISTANCE_A
    grid_density_sum = 0.0

    with path.open("rb") as handle:
        if byte_offset:
            handle.seek(byte_offset)

        linear_start = 0
        while linear_start < npoints:
            count = min(DENSITY_POINT_CHUNK_SIZE, npoints - linear_start)
            raw = handle.read(dtype.itemsize * count)
            values_native = np.frombuffer(raw, dtype=dtype)
            if values_native.size != count:
                raise ValueError(
                    "Unexpected EOF in %s at density value %d; expected %d more, found %d"
                    % (path, linear_start, count, values_native.size)
                )
            values = values_native.astype(np.float64, copy=False)
            grid_density_sum += float(np.sum(values, dtype=np.float64))

            # Convert flat BOV indices to coordinates without a 3-D meshgrid.
            linear = np.arange(linear_start, linear_start + count, dtype=np.int64)
            ix = linear % nx
            yz = linear // nx
            iy = yz % ny
            iz = yz // ny
            px = origin[0] + ix * dx
            py = origin[1] + iy * dy
            pz = origin[2] + iz * dz

            minimum2 = np.full(count, np.inf, dtype=np.float64)
            nearest = np.full(count, -1, dtype=np.int32)
            for atom_i, atom_pos in enumerate(atom_coords):
                delta_x = px - atom_pos[0]
                delta_y = py - atom_pos[1]
                delta_z = pz - atom_pos[2]
                distance2 = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z

                # <= matches the reference loop: the last atom wins exact ties.
                closer = distance2 <= minimum2
                minimum2[closer] = distance2[closer]
                nearest[closer] = atom_i

            valid = minimum2 <= cutoff2
            if np.any(valid):
                atom_density_sums += np.bincount(
                    nearest[valid], weights=values[valid], minlength=len(atom_coords)
                )
            linear_start += count

    voxel_volume = dx * dy * dz
    per_atom_electrons = atom_density_sums * voxel_volume
    grid_total_electrons = grid_density_sum * voxel_volume
    return per_atom_electrons, grid_total_electrons, voxel_volume


def nuclear_charges_for_species(species):
    missing = sorted(set(symbol for symbol in species if symbol not in NUCLEAR_CHARGE))
    if missing:
        raise ValueError("Missing NUCLEAR_CHARGE entries for: %s" % ", ".join(missing))
    return np.asarray([NUCLEAR_CHARGE[symbol] for symbol in species], dtype=float)


class ChargeResult:
    def __init__(
        self,
        requested_time,
        density_time,
        density_index,
        bov_path,
        bov_first_line,
        xyz_frame_index,
        xyz_time,
        atom_electrons,
        atom_charges,
        assigned_total_electrons,
        grid_total_electrons,
        system_charge,
        voxel_volume,
    ):
        self.requested_time = requested_time
        self.density_time = density_time
        self.density_index = density_index
        self.bov_path = bov_path
        self.bov_first_line = bov_first_line
        self.xyz_frame_index = xyz_frame_index
        self.xyz_time = xyz_time
        self.atom_electrons = atom_electrons
        self.atom_charges = atom_charges
        self.assigned_total_electrons = assigned_total_electrons
        self.grid_total_electrons = grid_total_electrons
        self.system_charge = system_charge
        self.voxel_volume = voxel_volume


def compute_charge_result(run_dir, trajectory, requested_time_fs):
    bov_path, density_time, density_index = bov_path_for_time(run_dir, requested_time_fs)
    if not bov_path.exists():
        raise FileNotFoundError("Missing BOV file: %s" % bov_path)

    xyz_frame_i = trajectory.closest_frame_index(density_time)
    first_line, data_path, grid_size, origin, brick_size, byte_offset = read_bov_header(bov_path)
    if not data_path.exists():
        raise FileNotFoundError("Missing DAT file referenced by %s: %s" % (bov_path, data_path))
    if PRINT_DENSITY_PROGRESS:
        print(
            "    Integrating %s (%d x %d x %d, %d atoms)"
            % (data_path, grid_size[0], grid_size[1], grid_size[2], len(trajectory.species)),
            flush=True,
        )
    atom_electrons, grid_total_electrons, voxel_volume = integrate_density_file_nearest_atom(
        path=data_path,
        grid_size=grid_size,
        byte_offset=byte_offset,
        origin=origin,
        brick_size=brick_size,
        atom_coords=trajectory.coords[xyz_frame_i],
    )
    z_values = nuclear_charges_for_species(trajectory.species)
    atom_charges = z_values - atom_electrons
    assigned_total_electrons = float(np.sum(atom_electrons))
    system_charge = float(np.sum(z_values) - assigned_total_electrons)

    return ChargeResult(
        requested_time=requested_time_fs,
        density_time=density_time,
        density_index=density_index,
        bov_path=bov_path,
        bov_first_line=first_line,
        xyz_frame_index=xyz_frame_i,
        xyz_time=float(trajectory.times[xyz_frame_i]),
        atom_electrons=atom_electrons,
        atom_charges=atom_charges,
        assigned_total_electrons=assigned_total_electrons,
        grid_total_electrons=grid_total_electrons,
        system_charge=system_charge,
        voxel_volume=voxel_volume,
    )


# ============================================================================
# Event analysis and output
# ============================================================================

def load_events(path):
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "run", "dir", "h_index", "old_owner_index", "new_owner_index",
            TRANSFER_TIME_COLUMN,
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError("Event CSV is missing columns: %s" % ", ".join(missing))
        rows = []
        for event_id, row in enumerate(reader, start=1):
            copied = dict(row)
            copied["_event_id"] = event_id
            rows.append(copied)
    return rows


def float_or_none(value):
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def join_floats(values, fmt="%.9f"):
    return ";".join(fmt % float(value) for value in values)


def join_indices_one_based(indices):
    return ";".join(str(int(index) + 1) for index in indices)


def four_frame_sample_times(before_anchor_time, after_anchor_time):
    if FOUR_FRAME_COUNT < 1:
        raise ValueError("FOUR_FRAME_COUNT must be at least 1")
    if FOUR_FRAME_SPACING_FS <= 0.0:
        raise ValueError("FOUR_FRAME_SPACING_FS must be positive")

    before = [
        before_anchor_time - FOUR_FRAME_SPACING_FS * offset
        for offset in range(FOUR_FRAME_COUNT, 0, -1)
    ]
    after = [
        after_anchor_time + FOUR_FRAME_SPACING_FS * offset
        for offset in range(1, FOUR_FRAME_COUNT + 1)
    ]
    return before, after


def make_event_base_row(event):
    return {
        "event_id": event["_event_id"],
        "run": event.get("run", ""),
        "dir": event.get("dir", ""),
        "h_index": event.get("h_index", ""),
        "donor_o_index": event.get("old_owner_index", ""),
        "acceptor_o_index": event.get("new_owner_index", ""),
        "transfer_time_column": TRANSFER_TIME_COLUMN,
        "transfer_time_fs": event.get(TRANSFER_TIME_COLUMN, ""),
    }


EVENT_OUTPUT_FIELDS = [
    "event_id", "run", "dir", "h_index", "donor_o_index", "acceptor_o_index",
    "transfer_time_column", "transfer_time_fs", "status", "reason",
    "before_anchor_time_fs", "after_anchor_time_fs",
    "before_target_oh_distance_A", "after_target_oh_distance_A",
    "anchor_separation_fs",
    "donor_charge_before_e", "donor_charge_after_e",
    "acceptor_charge_before_e", "acceptor_charge_after_e",
    "qH_from_donor_e", "qH_from_acceptor_e", "qH_effective_e",
    "donor_acceptor_mismatch_e", "pair_charge_before_e", "pair_charge_after_e",
    "pair_charge_change_e", "system_charge_before_e", "system_charge_after_e",
    "system_charge_change_e", "donor_acceptor_sanity_pass",
    "system_charge_sanity_pass", "overall_sanity_pass",
    "four_frame_mean_status", "four_frame_mean_reason",
    "four_frame_before_sample_times_fs", "four_frame_after_sample_times_fs",
    "donor_charge_before_four_frame_mean_e",
    "donor_charge_after_four_frame_mean_e",
    "acceptor_charge_before_four_frame_mean_e",
    "acceptor_charge_after_four_frame_mean_e",
    "donor_charge_before_four_frame_std_e",
    "donor_charge_after_four_frame_std_e",
    "acceptor_charge_before_four_frame_std_e",
    "acceptor_charge_after_four_frame_std_e",
    "qH_from_donor_four_frame_mean_e",
    "qH_from_acceptor_four_frame_mean_e",
    "qH_effective_four_frame_mean_e",
    "qH_effective_four_frame_std_e",
    "donor_acceptor_mismatch_four_frame_mean_e",
    "pair_charge_before_four_frame_mean_e",
    "pair_charge_after_four_frame_mean_e",
    "pair_charge_change_four_frame_mean_e",
    "system_charge_before_four_frame_mean_e",
    "system_charge_after_four_frame_mean_e",
    "system_charge_change_four_frame_mean_e",
    "donor_acceptor_four_frame_sanity_pass",
    "system_charge_four_frame_sanity_pass",
    "overall_four_frame_sanity_pass",
    "qH_effective_four_frame_minus_nearest_e",
    "qH_effective_four_frame_vs_nearest_abs_difference_e",
    "donor_H_before", "donor_H_after", "acceptor_H_before", "acceptor_H_after",
    "donor_H_ejected_before_excluded_both",
    "acceptor_H_ejected_before_excluded_both",
    "donor_H_ejected_between_kept_both",
    "acceptor_H_ejected_between_kept_both",
    "target_H_ejection_type", "target_H_ejection_time_fs",
]


FRAME_OUTPUT_FIELDS = [
    "event_id", "run", "dir", "sample_set", "sample_number", "side",
    "requested_time_fs",
    "density_time_fs", "xyz_time_fs", "density_index", "bov_file",
    "target_oh_distance_A",
    "donor_charge_e", "acceptor_charge_e", "pair_charge_e", "system_charge_e",
    "assigned_total_electrons", "grid_total_electrons", "unassigned_grid_electrons",
    "donor_H_members", "acceptor_H_members",
    "donor_H_ejected_between_kept", "acceptor_H_ejected_between_kept",
]


EJECTION_OUTPUT_FIELDS = [
    "run", "dir", "h_index", "ejection_type", "ejection_time_fs",
    "ejection_frame_index_zero_based", "nearest_o_index", "nearest_oh_distance_A",
    "h_in_cap", "nearest_o_in_cap", "tracking_stop_time_fs",
    "tracking_stop_reason",
]


def reject_result(event, reason, status="REJECTED"):
    row = make_event_base_row(event)
    row["status"] = status
    row["reason"] = reason
    return row


def other_transfer_in_window(event, run_events, window_start, window_end):
    for other in run_events:
        if other["_event_id"] == event["_event_id"]:
            continue
        other_time = float_or_none(other.get(TRANSFER_TIME_COLUMN))
        if other_time is None:
            continue
        if (
            window_start - EVENT_TIME_COMPARE_TOL_FS
            <= other_time
            <= window_end + EVENT_TIME_COMPARE_TOL_FS
        ):
            return other
    return None


def fragment_charge_at_frame(charge_result, donor_o, acceptor_o, donor_h, acceptor_h):
    donor_atoms = [donor_o] + list(donor_h)
    acceptor_atoms = [acceptor_o] + list(acceptor_h)
    donor_charge = float(np.sum(charge_result.atom_charges[donor_atoms]))
    acceptor_charge = float(np.sum(charge_result.atom_charges[acceptor_atoms]))
    return {
        "donor_charge": donor_charge,
        "acceptor_charge": acceptor_charge,
        "donor_h": list(donor_h),
        "acceptor_h": list(acceptor_h),
    }


def validate_four_frame_geometry(
    trajectory, sample_times, target_h, target_o, donor_o, acceptor_o,
):
    """Validate that every averaging sample represents the intended side."""
    frame_indices = []
    target_atoms = (target_h, donor_o, acceptor_o)
    for sample_time in sample_times:
        if sample_time + EVENT_TIME_COMPARE_TOL_FS < MIN_ANALYSIS_TIME_FS:
            raise ValueError(
                "four-frame sample %.6f fs is before MIN_ANALYSIS_TIME_FS %.6f fs"
                % (sample_time, MIN_ANALYSIS_TIME_FS)
            )
        frame_i = trajectory.closest_frame_index(sample_time)
        if REQUIRE_TARGET_OH_GEOMETRY_AT_FOUR_FRAME_SAMPLES:
            if not trajectory.geometry_is_clean_for_side(frame_i, target_h, target_o):
                raise ValueError(
                    "target O-H geometry is not valid at four-frame sample %.6f fs"
                    % sample_time
                )
        if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
            in_cap_atoms = [
                atom_i for atom_i in target_atoms
                if trajectory.atom_is_in_cap(frame_i, atom_i)
            ]
            if in_cap_atoms:
                raise ValueError(
                    "four-frame sample %.6f fs has target atom(s) in CAP: %s"
                    % (sample_time, join_indices_one_based(in_cap_atoms))
                )
        frame_indices.append(frame_i)
    return frame_indices


def calculate_charge_sample(
    event, run_dir, trajectory, charge_cache, ledger, donor_o, acceptor_o,
    side, sample_set, sample_number, sample_time, target_h,
):
    """Calculate one charge sample and its detailed output row."""
    donor_h = ledger["donor_%s_h" % side]
    acceptor_h = ledger["acceptor_%s_h" % side]
    _, _, density_index = bov_path_for_time(run_dir, sample_time)
    cache_key = (str(run_dir.resolve()), density_index)
    if cache_key not in charge_cache:
        charge_cache[cache_key] = compute_charge_result(
            run_dir, trajectory, sample_time
        )
    charge_result = charge_cache[cache_key]
    fragment = fragment_charge_at_frame(
        charge_result, donor_o, acceptor_o, donor_h, acceptor_h
    )
    xyz_frame_i = charge_result.xyz_frame_index
    target_distance = float(trajectory.nearest_o_distance[xyz_frame_i, target_h])

    frame_row = {
        "event_id": event["_event_id"],
        "run": event.get("run", ""),
        "dir": event.get("dir", ""),
        "sample_set": sample_set,
        "sample_number": sample_number,
        "side": side,
        "requested_time_fs": "%.9f" % sample_time,
        "density_time_fs": "%.9f" % charge_result.density_time,
        "xyz_time_fs": "%.9f" % charge_result.xyz_time,
        "density_index": charge_result.density_index,
        "bov_file": charge_result.bov_path.name,
        "target_oh_distance_A": "%.9f" % target_distance,
        "donor_charge_e": "%.12f" % fragment["donor_charge"],
        "acceptor_charge_e": "%.12f" % fragment["acceptor_charge"],
        "pair_charge_e": "%.12f" % (
            fragment["donor_charge"] + fragment["acceptor_charge"]
        ),
        "system_charge_e": "%.12f" % charge_result.system_charge,
        "assigned_total_electrons": "%.12f" % charge_result.assigned_total_electrons,
        "grid_total_electrons": "%.12f" % charge_result.grid_total_electrons,
        "unassigned_grid_electrons": "%.12f" % (
            charge_result.grid_total_electrons
            - charge_result.assigned_total_electrons
        ),
        "donor_H_members": join_indices_one_based(fragment["donor_h"]),
        "acceptor_H_members": join_indices_one_based(fragment["acceptor_h"]),
        "donor_H_ejected_between_kept": join_indices_one_based(
            ledger["ejected_between_donor"]
        ),
        "acceptor_H_ejected_between_kept": join_indices_one_based(
            ledger["ejected_between_acceptor"]
        ),
    }
    return {
        "fragment": fragment,
        "charge_result": charge_result,
        "frame_row": frame_row,
    }


def summarize_two_sides(before_samples, after_samples):
    """Summarize either one nearest sample or equal-sized averaging sets."""
    if not before_samples or not after_samples:
        raise ValueError("both before and after samples are required")
    if len(before_samples) != len(after_samples):
        raise ValueError("before and after sample counts differ")

    donor_before_values = np.asarray([
        sample["fragment"]["donor_charge"] for sample in before_samples
    ], dtype=float)
    donor_after_values = np.asarray([
        sample["fragment"]["donor_charge"] for sample in after_samples
    ], dtype=float)
    acceptor_before_values = np.asarray([
        sample["fragment"]["acceptor_charge"] for sample in before_samples
    ], dtype=float)
    acceptor_after_values = np.asarray([
        sample["fragment"]["acceptor_charge"] for sample in after_samples
    ], dtype=float)
    system_before_values = np.asarray([
        sample["charge_result"].system_charge for sample in before_samples
    ], dtype=float)
    system_after_values = np.asarray([
        sample["charge_result"].system_charge for sample in after_samples
    ], dtype=float)

    donor_before = float(np.mean(donor_before_values))
    donor_after = float(np.mean(donor_after_values))
    acceptor_before = float(np.mean(acceptor_before_values))
    acceptor_after = float(np.mean(acceptor_after_values))
    system_before = float(np.mean(system_before_values))
    system_after = float(np.mean(system_after_values))

    qh_donor = donor_before - donor_after
    qh_acceptor = acceptor_after - acceptor_before
    qh_effective = 0.5 * (qh_donor + qh_acceptor)
    paired_qh_effective = 0.5 * (
        (donor_before_values - donor_after_values)
        + (acceptor_after_values - acceptor_before_values)
    )
    mismatch = abs(qh_donor - qh_acceptor)
    pair_before = donor_before + acceptor_before
    pair_after = donor_after + acceptor_after
    pair_change = pair_after - pair_before
    system_change = system_after - system_before
    da_pass = mismatch <= DONOR_ACCEPTOR_MISMATCH_TOL_E
    system_pass = abs(system_change) <= SYSTEM_CHARGE_CHANGE_TOL_E

    return {
        "donor_before": donor_before,
        "donor_after": donor_after,
        "acceptor_before": acceptor_before,
        "acceptor_after": acceptor_after,
        "donor_before_std": float(np.std(donor_before_values)),
        "donor_after_std": float(np.std(donor_after_values)),
        "acceptor_before_std": float(np.std(acceptor_before_values)),
        "acceptor_after_std": float(np.std(acceptor_after_values)),
        "qh_donor": qh_donor,
        "qh_acceptor": qh_acceptor,
        "qh_effective": qh_effective,
        "qh_effective_std": float(np.std(paired_qh_effective)),
        "mismatch": mismatch,
        "pair_before": pair_before,
        "pair_after": pair_after,
        "pair_change": pair_change,
        "system_before": system_before,
        "system_after": system_after,
        "system_change": system_change,
        "da_pass": da_pass,
        "system_pass": system_pass,
        "overall_pass": da_pass and system_pass,
    }


def analyze_event(event, run_events, run_dir, trajectory, charge_cache):
    transfer_time = float_or_none(event.get(TRANSFER_TIME_COLUMN))
    if transfer_time is None:
        return reject_result(event, "missing transfer time"), []

    try:
        target_h = int(event["h_index"]) - 1
        donor_o = int(event["old_owner_index"]) - 1
        acceptor_o = int(event["new_owner_index"]) - 1
    except Exception:
        return reject_result(event, "invalid atom index in event CSV"), []

    natoms = len(trajectory.species)
    if not (0 <= target_h < natoms and 0 <= donor_o < natoms and 0 <= acceptor_o < natoms):
        return reject_result(event, "event atom index outside trajectory"), []
    if trajectory.species[target_h] != "H":
        return reject_result(event, "h_index does not identify H"), []
    if trajectory.species[donor_o] != "O" or trajectory.species[acceptor_o] != "O":
        return reject_result(event, "donor/acceptor index does not identify O"), []
    if donor_o == acceptor_o:
        return reject_result(event, "donor and acceptor O are identical"), []

    before_anchor_i = trajectory.find_anchor(
        transfer_time, target_h, donor_o, direction="before"
    )
    if before_anchor_i is None:
        return reject_result(event, "no clean pre-transfer geometry anchor"), []
    after_anchor_i = trajectory.find_anchor(
        transfer_time, target_h, acceptor_o, direction="after"
    )
    if after_anchor_i is None:
        return reject_result(event, "no clean post-transfer geometry anchor"), []

    before_anchor_time = float(trajectory.times[before_anchor_i])
    after_anchor_time = float(trajectory.times[after_anchor_i])
    if after_anchor_time <= before_anchor_time + EVENT_TIME_COMPARE_TOL_FS:
        return reject_result(event, "post-transfer anchor is not later than pre-transfer anchor"), []
    before_target_distance = float(
        trajectory.nearest_o_distance[before_anchor_i, target_h]
    )
    after_target_distance = float(
        trajectory.nearest_o_distance[after_anchor_i, target_h]
    )

    if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
        target_atoms = (target_h, donor_o, acceptor_o)
        for side, frame_i in (("before", before_anchor_i), ("after", after_anchor_i)):
            in_cap_atoms = [atom_i for atom_i in target_atoms if trajectory.atom_is_in_cap(frame_i, atom_i)]
            if in_cap_atoms:
                return reject_result(
                    event,
                    "%s anchor has target atom(s) in CAP: %s"
                    % (side, join_indices_one_based(in_cap_atoms)),
                ), []

    analysis_start = before_anchor_time
    analysis_end = after_anchor_time
    if REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW:
        conflict = other_transfer_in_window(
            event, run_events, analysis_start, analysis_end
        )
        if conflict is not None:
            reason = "other H transfer event %s occurs inside analysis window" % conflict["_event_id"]
            return reject_result(event, reason), []

    ledger = trajectory.build_fixed_fragment_ledger(
        before_anchor_i, after_anchor_i, target_h, donor_o, acceptor_o
    )

    # Nearest-anchor estimate.  Existing output names continue to refer to this
    # estimate so that older downstream analysis remains compatible.
    frame_rows = []
    nearest_samples = {"before": [], "after": []}
    try:
        for side, sample_time in (
            ("before", before_anchor_time),
            ("after", after_anchor_time),
        ):
            sample = calculate_charge_sample(
                event, run_dir, trajectory, charge_cache, ledger,
                donor_o, acceptor_o, side, "nearest_anchor", 1,
                sample_time, target_h,
            )
            nearest_samples[side].append(sample)
            frame_rows.append(sample["frame_row"])
        nearest = summarize_two_sides(
            nearest_samples["before"], nearest_samples["after"]
        )
    except Exception as exc:
        return reject_result(
            event, "nearest-frame charge calculation failed: %s" % exc,
            status="ERROR",
        ), []

    # Four-frame estimate.  The same fixed fragment ledger and ejection rules
    # are used; only the density sampling times differ from the nearest result.
    four_before_times, four_after_times = four_frame_sample_times(
        before_anchor_time, after_anchor_time
    )
    four_frame_status = "ANALYZED"
    four_frame_reason = ""
    four_frame = None
    try:
        validate_four_frame_geometry(
            trajectory, four_before_times, target_h, donor_o, donor_o, acceptor_o
        )
        validate_four_frame_geometry(
            trajectory, four_after_times, target_h, acceptor_o, donor_o, acceptor_o
        )

        if REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW:
            conflict = other_transfer_in_window(
                event, run_events, min(four_before_times), max(four_after_times)
            )
            if conflict is not None:
                raise ValueError(
                    "other H transfer event %s occurs inside the four-frame window"
                    % conflict["_event_id"]
                )

        four_samples = {"before": [], "after": []}
        four_detail_rows = []
        for side, sample_times in (
            ("before", four_before_times),
            ("after", four_after_times),
        ):
            for sample_number, sample_time in enumerate(sample_times, start=1):
                sample = calculate_charge_sample(
                    event, run_dir, trajectory, charge_cache, ledger,
                    donor_o, acceptor_o, side,
                    "four_frame_average_member", sample_number,
                    sample_time, target_h,
                )
                four_samples[side].append(sample)
                four_detail_rows.append(sample["frame_row"])
        four_frame = summarize_two_sides(
            four_samples["before"], four_samples["after"]
        )
        frame_rows.extend(four_detail_rows)
    except Exception as exc:
        four_frame_status = "UNAVAILABLE"
        four_frame_reason = str(exc)

    row = make_event_base_row(event)
    row.update({
        "status": "ANALYZED",
        "reason": "",
        "before_anchor_time_fs": "%.9f" % before_anchor_time,
        "after_anchor_time_fs": "%.9f" % after_anchor_time,
        "before_target_oh_distance_A": "%.9f" % before_target_distance,
        "after_target_oh_distance_A": "%.9f" % after_target_distance,
        "anchor_separation_fs": "%.9f" % (after_anchor_time - before_anchor_time),
        "donor_charge_before_e": "%.12f" % nearest["donor_before"],
        "donor_charge_after_e": "%.12f" % nearest["donor_after"],
        "acceptor_charge_before_e": "%.12f" % nearest["acceptor_before"],
        "acceptor_charge_after_e": "%.12f" % nearest["acceptor_after"],
        "qH_from_donor_e": "%.12f" % nearest["qh_donor"],
        "qH_from_acceptor_e": "%.12f" % nearest["qh_acceptor"],
        "qH_effective_e": "%.12f" % nearest["qh_effective"],
        "donor_acceptor_mismatch_e": "%.12f" % nearest["mismatch"],
        "pair_charge_before_e": "%.12f" % nearest["pair_before"],
        "pair_charge_after_e": "%.12f" % nearest["pair_after"],
        "pair_charge_change_e": "%.12f" % nearest["pair_change"],
        "system_charge_before_e": "%.12f" % nearest["system_before"],
        "system_charge_after_e": "%.12f" % nearest["system_after"],
        "system_charge_change_e": "%.12f" % nearest["system_change"],
        "donor_acceptor_sanity_pass": int(nearest["da_pass"]),
        "system_charge_sanity_pass": int(nearest["system_pass"]),
        "overall_sanity_pass": int(nearest["overall_pass"]),
        "four_frame_mean_status": four_frame_status,
        "four_frame_mean_reason": four_frame_reason,
        "four_frame_before_sample_times_fs": join_floats(four_before_times),
        "four_frame_after_sample_times_fs": join_floats(four_after_times),
        "donor_H_before": join_indices_one_based(ledger["donor_before_h"]),
        "donor_H_after": join_indices_one_based(ledger["donor_after_h"]),
        "acceptor_H_before": join_indices_one_based(ledger["acceptor_before_h"]),
        "acceptor_H_after": join_indices_one_based(ledger["acceptor_after_h"]),
        "donor_H_ejected_before_excluded_both": join_indices_one_based(
            ledger["excluded_ejected_before_donor"]
        ),
        "acceptor_H_ejected_before_excluded_both": join_indices_one_based(
            ledger["excluded_ejected_before_acceptor"]
        ),
        "donor_H_ejected_between_kept_both": join_indices_one_based(
            ledger["ejected_between_donor"]
        ),
        "acceptor_H_ejected_between_kept_both": join_indices_one_based(
            ledger["ejected_between_acceptor"]
        ),
        "target_H_ejection_type": trajectory.ejection_type.get(target_h) or "",
        "target_H_ejection_time_fs": (
            "%.9f" % trajectory.ejection_time[target_h]
            if trajectory.ejection_time.get(target_h) is not None else ""
        ),
    })

    if four_frame is not None:
        method_difference = four_frame["qh_effective"] - nearest["qh_effective"]
        row.update({
            "donor_charge_before_four_frame_mean_e": "%.12f" % four_frame["donor_before"],
            "donor_charge_after_four_frame_mean_e": "%.12f" % four_frame["donor_after"],
            "acceptor_charge_before_four_frame_mean_e": "%.12f" % four_frame["acceptor_before"],
            "acceptor_charge_after_four_frame_mean_e": "%.12f" % four_frame["acceptor_after"],
            "donor_charge_before_four_frame_std_e": "%.12f" % four_frame["donor_before_std"],
            "donor_charge_after_four_frame_std_e": "%.12f" % four_frame["donor_after_std"],
            "acceptor_charge_before_four_frame_std_e": "%.12f" % four_frame["acceptor_before_std"],
            "acceptor_charge_after_four_frame_std_e": "%.12f" % four_frame["acceptor_after_std"],
            "qH_from_donor_four_frame_mean_e": "%.12f" % four_frame["qh_donor"],
            "qH_from_acceptor_four_frame_mean_e": "%.12f" % four_frame["qh_acceptor"],
            "qH_effective_four_frame_mean_e": "%.12f" % four_frame["qh_effective"],
            "qH_effective_four_frame_std_e": "%.12f" % four_frame["qh_effective_std"],
            "donor_acceptor_mismatch_four_frame_mean_e": "%.12f" % four_frame["mismatch"],
            "pair_charge_before_four_frame_mean_e": "%.12f" % four_frame["pair_before"],
            "pair_charge_after_four_frame_mean_e": "%.12f" % four_frame["pair_after"],
            "pair_charge_change_four_frame_mean_e": "%.12f" % four_frame["pair_change"],
            "system_charge_before_four_frame_mean_e": "%.12f" % four_frame["system_before"],
            "system_charge_after_four_frame_mean_e": "%.12f" % four_frame["system_after"],
            "system_charge_change_four_frame_mean_e": "%.12f" % four_frame["system_change"],
            "donor_acceptor_four_frame_sanity_pass": int(four_frame["da_pass"]),
            "system_charge_four_frame_sanity_pass": int(four_frame["system_pass"]),
            "overall_four_frame_sanity_pass": int(four_frame["overall_pass"]),
            "qH_effective_four_frame_minus_nearest_e": "%.12f" % method_difference,
            "qH_effective_four_frame_vs_nearest_abs_difference_e": "%.12f" % abs(method_difference),
        })
    return row, frame_rows


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ejection_rows_for_trajectory(trajectory, run, directory_name):
    rows = []
    for h_index_value in trajectory.hydrogen_indices:
        h_index = int(h_index_value)
        if trajectory.ejection_time.get(h_index) is None:
            continue
        rows.append({
            "run": run,
            "dir": directory_name,
            "h_index": h_index + 1,
            "ejection_type": trajectory.ejection_type[h_index],
            "ejection_time_fs": "%.9f" % trajectory.ejection_time[h_index],
            "ejection_frame_index_zero_based": trajectory.ejection_frame_index[h_index],
            "nearest_o_index": trajectory.ejection_nearest_o[h_index] + 1,
            "nearest_oh_distance_A": "%.9f" % trajectory.ejection_distance[h_index],
            "h_in_cap": int(trajectory.ejection_h_in_cap[h_index]),
            "nearest_o_in_cap": int(trajectory.ejection_o_in_cap[h_index]),
            "tracking_stop_time_fs": "%.9f" % trajectory.ejection_tracking_stop_time[h_index],
            "tracking_stop_reason": trajectory.ejection_tracking_stop_reason[h_index],
        })
    return rows


def write_summary(path, event_rows, total_input_events, trajectories_loaded, cache_size):
    status_counts = Counter(row.get("status", "") for row in event_rows)
    analyzed = [row for row in event_rows if row.get("status") == "ANALYZED"]
    sanity_pass = sum(int(row.get("overall_sanity_pass", 0)) for row in analyzed)

    qh_values = [float(row["qH_effective_e"]) for row in analyzed]
    four_frame_analyzed = [
        row for row in analyzed
        if row.get("four_frame_mean_status") == "ANALYZED"
        and row.get("qH_effective_four_frame_mean_e", "") != ""
    ]
    four_frame_qh_values = [
        float(row["qH_effective_four_frame_mean_e"])
        for row in four_frame_analyzed
    ]
    method_abs_differences = [
        float(row["qH_effective_four_frame_vs_nearest_abs_difference_e"])
        for row in four_frame_analyzed
    ]
    summary = {
        "input_events": total_input_events,
        "trajectories_loaded": trajectories_loaded,
        "analyzed_events": len(analyzed),
        "rejected_events": status_counts.get("REJECTED", 0),
        "error_events": status_counts.get("ERROR", 0),
        "sanity_pass_events": sanity_pass,
        "sanity_pass_fraction": (sanity_pass / len(analyzed)) if analyzed else "",
        "mean_qH_effective_e": float(np.mean(qh_values)) if qh_values else "",
        "std_qH_effective_e": float(np.std(qh_values)) if qh_values else "",
        "min_qH_effective_e": min(qh_values) if qh_values else "",
        "max_qH_effective_e": max(qh_values) if qh_values else "",
        "four_frame_analyzed_events": len(four_frame_analyzed),
        "four_frame_unavailable_events": len(analyzed) - len(four_frame_analyzed),
        "mean_qH_effective_four_frame_mean_e": (
            float(np.mean(four_frame_qh_values)) if four_frame_qh_values else ""
        ),
        "std_qH_effective_four_frame_mean_e": (
            float(np.std(four_frame_qh_values)) if four_frame_qh_values else ""
        ),
        "mean_abs_four_frame_vs_nearest_difference_e": (
            float(np.mean(method_abs_differences)) if method_abs_differences else ""
        ),
        "max_abs_four_frame_vs_nearest_difference_e": (
            max(method_abs_differences) if method_abs_differences else ""
        ),
        "unique_density_frames_integrated": cache_size,
    }
    write_csv(path, list(summary.keys()), [summary])
    return summary, status_counts


def main():
    base_dir = BASE_DIR.resolve()
    event_csv = base_dir / EVENTS_CSV_NAME
    if not event_csv.exists():
        raise FileNotFoundError("Event CSV not found: %s" % event_csv)

    events = load_events(event_csv)
    events_by_dir = defaultdict(list)
    for event in events:
        events_by_dir[event.get("dir", "")].append(event)

    event_rows = []
    frame_rows = []
    ejection_rows = []
    charge_cache = {}
    trajectories_loaded = 0

    print("=== Clean H-transfer charge analysis ===", flush=True)
    print("Event time column:", TRANSFER_TIME_COLUMN, flush=True)
    print("Input events:", len(events), flush=True)
    print("Density output interval [fs]:", DENSITY_OUTPUT_DT_FS, flush=True)

    for directory_name, run_events in events_by_dir.items():
        if not directory_name:
            for event in run_events:
                event_rows.append(reject_result(event, "missing run directory in event CSV", status="ERROR"))
            continue

        run_dir = base_dir / directory_name
        xyz_path = run_dir / XYZ_NAME
        if not xyz_path.exists():
            for event in run_events:
                event_rows.append(reject_result(event, "missing trajectory: %s" % xyz_path, status="ERROR"))
            continue

        try:
            trajectory = read_xyz_trajectory(xyz_path)
            trajectories_loaded += 1
            ejection_rows.extend(ejection_rows_for_trajectory(
                trajectory=trajectory,
                run=run_events[0].get("run", "") if run_events else "",
                directory_name=directory_name,
            ))
        except Exception as exc:
            for event in run_events:
                event_rows.append(reject_result(event, "trajectory read failed: %s" % exc, status="ERROR"))
            continue

        print("Processing %s: %d event(s)" % (directory_name, len(run_events)), flush=True)
        for event in run_events:
            result_row, result_frame_rows = analyze_event(
                event=event,
                run_events=run_events,
                run_dir=run_dir,
                trajectory=trajectory,
                charge_cache=charge_cache,
            )
            event_rows.append(result_row)
            frame_rows.extend(result_frame_rows)
            print(
                "  event %s H%s O%s->O%s: %s%s"
                % (
                    event["_event_id"], event.get("h_index", "?"),
                    event.get("old_owner_index", "?"), event.get("new_owner_index", "?"),
                    result_row.get("status", ""),
                    (" (" + result_row.get("reason", "") + ")") if result_row.get("reason") else "",
                ),
                flush=True,
            )

    out_events = base_dir / OUT_EVENTS_CSV
    out_frames = base_dir / OUT_FRAMES_CSV
    out_summary = base_dir / OUT_SUMMARY_CSV
    out_ejections = base_dir / OUT_EJECTIONS_CSV
    write_csv(out_events, EVENT_OUTPUT_FIELDS, event_rows)
    write_csv(out_frames, FRAME_OUTPUT_FIELDS, frame_rows)
    write_csv(out_ejections, EJECTION_OUTPUT_FIELDS, ejection_rows)
    summary, status_counts = write_summary(
        out_summary, event_rows, len(events), trajectories_loaded, len(charge_cache)
    )

    print("=== Overall ===", flush=True)
    for status, count in sorted(status_counts.items()):
        print("%-10s %d" % (status + ":", count), flush=True)
    print("Sanity-pass analyzed events: %s/%s" % (
        summary["sanity_pass_events"], summary["analyzed_events"]
    ), flush=True)
    print("Wrote:", out_events, flush=True)
    print("Wrote:", out_frames, flush=True)
    print("Wrote:", out_summary, flush=True)
    print("Wrote:", out_ejections, flush=True)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Estimate the effective charge carried in clean H-transfer events.

This script is designed to run after the existing H-transfer detector.  It
reads h_rearrangement_events.csv, finds the nearest pre- and post-transfer
frames at which the transferred H is covalently associated with the expected
O, and reports both (1) the charge estimate from those two nearest frames and
(2) the estimate from four-frame averages immediately outside the anchors.
The two estimates are written together so that their sampling sensitivity can
be assessed directly.

The density partition is the same nearest-nucleus ("Cody") real-space method
used in density_integrate.py: every grid point within MAX_DENSITY_DISTANCE_A of
an atom is assigned to the nearest atom.  Atomic charges are Z - N_e.

Fragment bookkeeping uses fixed atom identities.  At the pre-transfer frame,
each non-transferred H assigned to the donor or acceptor is either:
  * excluded from both frames if it had already met the H-ejection definition;
  * included in the same parent fragment at both frames otherwise, even if it
    ejects between the two frames.
The transferred H is the only atom moved from donor to acceptor.  Thus H
ejection does not create a trivial change in fragment nuclear composition.

For corrected donor/acceptor fragment charges Q_D and Q_A,

    q_H(donor)    = Q_D(before) - Q_D(after)
    q_H(acceptor) = Q_A(after)  - Q_A(before)

Agreement between these two values is the main internal sanity check.

Run from the directory that contains the tdl* folders and the event CSV:

    python clean_h_transfer_charge_analysis.py

Only NumPy and the Python standard library are required.
"""

import csv
import faulthandler
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# If a compiled numerical library or the operating system terminates Python,
# print the active Python stack instead of only the words "Segmentation fault".
faulthandler.enable(all_threads=True)

import numpy as np


# ============================================================================
# USER-TUNABLE PARAMETERS
# ============================================================================

# Input layout
BASE_DIR = Path(".")
EVENTS_CSV_NAME = "h_rearrangement_events.csv"
XYZ_NAME = "trajectory.xyz"
DENSITY_SUBDIR = ""                 # e.g. "density" if dens*.bov is elsewhere
BOV_PREFIX = "dens"
BOV_DIGITS = 5

# Which time reported by the previous script should define the transfer time.
# Usually start_time_fs is the first frame captured by the new O.  Alternatives
# in the existing CSV are stable_time_fs and confirm_time_fs.
TRANSFER_TIME_COLUMN = "start_time_fs"

# Geometry criterion for the transferred O-H only (Angstrom).  Other O-H bond
# lengths do not affect anchor selection.
TRANSFERRED_OH_BOND_CUTOFF_A = 1.5

# H-ejection definition, matched to h_ejection_by_run.py.
# Clean: first frame with nearest-O distance > 3.0 A before either relevant
# atom reaches the CAP.  CAP special case: at the first frame at which H or its
# nearest O reaches the CAP, count ejection only if nearest-O distance > 2.5 A.
CAP_BOUND_A = 11.0
H_EJECTION_CLEAN_DISTANCE_A = 3.0
H_EJECTION_CAP_DISTANCE_A = 2.5

# A target H that was already classified as ejected and was later recaptured
# is not a clean direct H-transfer event.  Keep this True for the main analysis.
REJECT_IF_TRANSFERRED_H_WAS_EJECTED = True

# Density partitioning close to the CAP is difficult to interpret.  Reject an
# anchor if the transferred H, donor O, or acceptor O is already in the CAP.
REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS = True

# Both single-frame charge samples must be at/after this time, when the
# principal ionization stage is assumed to be over.
MIN_ANALYSIS_TIME_FS = 17.5

# Search no farther than this from the transfer time for a clean anchor.  Set
# to None to search all available trajectory frames.
MAX_ANCHOR_SEARCH_FS = None

# A clean event may contain H ejection, but no second H-transfer event may
# begin from the selected pre-transfer frame through the selected
# post-transfer frame.
REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW = True

# Time matching.  Charge sample times should normally coincide with saved XYZ
# and density frames.  Increase only if output times have small rounding noise.
XYZ_TIME_MATCH_TOL_FS = 1.0e-4
DENSITY_TIME_MATCH_TOL_FS = 1.0e-6
EVENT_TIME_COMPARE_TOL_FS = 1.0e-6

# Density-file timing, copied from density_integrate(2).py.  This controls the
# mapping time -> densNNNNN.bov.
SIMULATION_TIME_STEP_FS = 0.001
TD_DENSITY_OUTPUT_FREQUENCY = 500
DENSITY_OUTPUT_DT_FS = SIMULATION_TIME_STEP_FS * TD_DENSITY_OUTPUT_FREQUENCY

# Four-frame comparison.  With the current 0.5 fs density interval, the before
# samples are anchor-2.0, -1.5, -1.0, and -0.5 fs, while the after samples are
# anchor+0.5, +1.0, +1.5, and +2.0 fs.  The anchors themselves are reported
# separately as the nearest-frame result.
FOUR_FRAME_COUNT = 4
FOUR_FRAME_SPACING_FS = DENSITY_OUTPUT_DT_FS
REQUIRE_TARGET_OH_GEOMETRY_AT_FOUR_FRAME_SAMPLES = True

# Nearest-nucleus density partition (same physical method as the uploaded
# density integration script).
MAX_DENSITY_DISTANCE_A = 10.0
DENSITY_BINARY_DTYPE = "<f4"         # little-endian float32

# Hard upper bound on the number of grid points handled at once.  The original
# vectorized implementation formed a (points x atoms x 3) temporary array,
# which could exhaust memory and be reported by some HPC systems as a bare
# segmentation fault.  The streaming implementation below loops over atoms and
# therefore uses O(points), rather than O(points * atoms), temporary memory.
DENSITY_POINT_CHUNK_SIZE = 100_000

# Print each density file immediately before it is integrated.  Keeping this on
# is useful on batch systems, whose buffered stdout can otherwise make a later
# failure look as if it happened at program startup.
PRINT_DENSITY_PROGRESS = True

# Nuclear charges needed for water.  Add entries here for other elements.
NUCLEAR_CHARGE = {"H": 1.0, "O": 8.0}

# Sanity-check tolerances (elementary charge, e)
DONOR_ACCEPTOR_MISMATCH_TOL_E = 1.25
SYSTEM_CHARGE_CHANGE_TOL_E = 1.25

# Outputs
OUT_EVENTS_CSV = "clean_h_transfer_charges.csv"
OUT_FRAMES_CSV = "clean_h_transfer_charge_frames.csv"
OUT_SUMMARY_CSV = "clean_h_transfer_charge_summary.csv"
OUT_EJECTIONS_CSV = "clean_h_transfer_detected_ejections.csv"


TIME_RE = re.compile(r"time\[fs\]\s*=\s*([0-9Ee+\-.]+)")
ITER_RE = re.compile(r"iter\s*=\s*([0-9]+)")


# ============================================================================
# XYZ trajectory and geometry analysis
# ============================================================================

def parse_xyz_time(comment, iter_index):
    match = TIME_RE.search(comment)
    if match:
        return float(match.group(1))
    if iter_index is not None:
        return iter_index * SIMULATION_TIME_STEP_FS
    raise ValueError("XYZ frame has neither time[fs] nor iter")


def read_xyz_trajectory(path):
    """Read trajectory.xyz while preserving the fixed atom ordering."""
    times = []
    iterations = []
    all_coords = []
    reference_species = None

    with path.open("r") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.strip():
                continue
            natoms = int(line.strip())
            comment = handle.readline()
            if not comment:
                raise RuntimeError("EOF after XYZ atom count in %s" % path)

            iter_match = ITER_RE.search(comment)
            iter_index = int(iter_match.group(1)) if iter_match else None
            time_fs = parse_xyz_time(comment, iter_index)

            species = []
            coords = np.empty((natoms, 3), dtype=float)
            for atom_i in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise RuntimeError("EOF inside XYZ frame in %s" % path)
                parts = atom_line.split()
                if len(parts) < 4:
                    raise RuntimeError("Bad XYZ atom line: %s" % atom_line.strip())
                species.append(parts[0].capitalize())
                coords[atom_i] = [float(parts[1]), float(parts[2]), float(parts[3])]

            if reference_species is None:
                reference_species = species
            elif species != reference_species:
                raise ValueError("Atom identities/order change between frames in %s" % path)

            times.append(time_fs)
            iterations.append(iter_index)
            all_coords.append(coords)

    if not all_coords:
        raise ValueError("No XYZ frames found in %s" % path)

    times_array = np.asarray(times, dtype=float)
    if np.any(np.diff(times_array) <= 0.0):
        raise ValueError("XYZ times are not strictly increasing in %s" % path)

    return Trajectory(
        path=path,
        times=times_array,
        iterations=iterations,
        species=reference_species,
        coords=np.stack(all_coords, axis=0),
    )


class Trajectory:
    def __init__(self, path, times, iterations, species, coords):
        self.path = path
        self.times = times
        self.iterations = iterations
        self.species = species
        self.coords = coords
        self.oxygen_indices = np.asarray(
            [i for i, symbol in enumerate(species) if symbol == "O"], dtype=int
        )
        self.hydrogen_indices = np.asarray(
            [i for i, symbol in enumerate(species) if symbol == "H"], dtype=int
        )
        if self.oxygen_indices.size == 0 or self.hydrogen_indices.size == 0:
            raise ValueError("Trajectory must contain both O and H atoms: %s" % path)

        self.nearest_o = np.full(
            (len(times), len(species)), -1, dtype=int
        )
        self.nearest_o_distance = np.full(
            (len(times), len(species)), np.inf, dtype=float
        )
        self._calculate_nearest_oxygen()
        self.ejection_time = {}
        self.ejection_frame_index = {}
        self.ejection_type = {}
        self.ejection_nearest_o = {}
        self.ejection_distance = {}
        self.ejection_h_in_cap = {}
        self.ejection_o_in_cap = {}
        self.ejection_tracking_stop_time = {}
        self.ejection_tracking_stop_reason = {}
        self._detect_h_ejections()

    def _calculate_nearest_oxygen(self):
        oxygen_coords = self.coords[:, self.oxygen_indices, :]
        for h_index in self.hydrogen_indices:
            displacement = oxygen_coords - self.coords[:, h_index, :][:, None, :]
            distances = np.sqrt(np.sum(displacement * displacement, axis=2))
            nearest_local = np.argmin(distances, axis=1)
            rows = np.arange(len(self.times))
            self.nearest_o[:, h_index] = self.oxygen_indices[nearest_local]
            self.nearest_o_distance[:, h_index] = distances[rows, nearest_local]

    def _detect_h_ejections(self):
        """Apply the exact first-crossing/CAP logic of h_ejection_by_run.py."""
        for h_index in self.hydrogen_indices:
            h_index = int(h_index)
            self.ejection_time[h_index] = None
            self.ejection_frame_index[h_index] = None
            self.ejection_type[h_index] = None
            self.ejection_nearest_o[h_index] = None
            self.ejection_distance[h_index] = None
            self.ejection_h_in_cap[h_index] = None
            self.ejection_o_in_cap[h_index] = None
            self.ejection_tracking_stop_time[h_index] = None
            self.ejection_tracking_stop_reason[h_index] = None

            for frame_i in range(len(self.times)):
                nearest_o = int(self.nearest_o[frame_i, h_index])
                nearest_distance = float(self.nearest_o_distance[frame_i, h_index])
                h_in_cap = self.atom_is_in_cap(frame_i, h_index)
                o_in_cap = self.atom_is_in_cap(frame_i, nearest_o)

                # The reference script stops tracking at the first CAP touch,
                # whether or not the 2.5 A fallback criterion is satisfied.
                if h_in_cap or o_in_cap:
                    self.ejection_tracking_stop_time[h_index] = float(self.times[frame_i])
                    if nearest_distance > H_EJECTION_CAP_DISTANCE_A:
                        self._record_ejection(
                            h_index, frame_i, "CAP", nearest_o,
                            nearest_distance, h_in_cap, o_in_cap,
                        )
                        self.ejection_tracking_stop_reason[h_index] = "CAP_EJECTION"
                    else:
                        self.ejection_tracking_stop_reason[h_index] = "CAP_NO_EJECTION"
                    break

                if nearest_distance > H_EJECTION_CLEAN_DISTANCE_A:
                    self._record_ejection(
                        h_index, frame_i, "CLEAN", nearest_o,
                        nearest_distance, h_in_cap, o_in_cap,
                    )
                    self.ejection_tracking_stop_time[h_index] = float(self.times[frame_i])
                    self.ejection_tracking_stop_reason[h_index] = "CLEAN_EJECTION"
                    break

    def _record_ejection(
        self, h_index, frame_i, event_type, nearest_o, nearest_distance,
        h_in_cap, o_in_cap,
    ):
        self.ejection_time[h_index] = float(self.times[frame_i])
        self.ejection_frame_index[h_index] = int(frame_i)
        self.ejection_type[h_index] = event_type
        self.ejection_nearest_o[h_index] = int(nearest_o)
        self.ejection_distance[h_index] = float(nearest_distance)
        self.ejection_h_in_cap[h_index] = bool(h_in_cap)
        self.ejection_o_in_cap[h_index] = bool(o_in_cap)

    def atom_is_in_cap(self, frame_i, atom_i):
        x, y, z = self.coords[frame_i, atom_i]
        return (
            abs(float(x)) >= CAP_BOUND_A
            or abs(float(y)) >= CAP_BOUND_A
            or abs(float(z)) >= CAP_BOUND_A
        )

    def closest_frame_index(self, target_time_fs, tolerance_fs=XYZ_TIME_MATCH_TOL_FS):
        insertion = int(np.searchsorted(self.times, target_time_fs))
        candidates = []
        if insertion < len(self.times):
            candidates.append(insertion)
        if insertion > 0:
            candidates.append(insertion - 1)
        if not candidates:
            raise ValueError("No XYZ frames available")
        best = min(candidates, key=lambda i: abs(self.times[i] - target_time_fs))
        mismatch = abs(self.times[best] - target_time_fs)
        if mismatch > tolerance_fs:
            raise ValueError(
                "No XYZ frame at %.9f fs; closest is %.9f fs (difference %.9g fs)"
                % (target_time_fs, self.times[best], mismatch)
            )
        return best

    def h_is_ejected_by(self, h_index, time_fs):
        eject_time = self.ejection_time.get(h_index)
        return eject_time is not None and eject_time <= time_fs + EVENT_TIME_COMPARE_TOL_FS

    def geometry_is_clean_for_side(self, frame_i, target_h, target_o):
        if self.times[frame_i] + EVENT_TIME_COMPARE_TOL_FS < MIN_ANALYSIS_TIME_FS:
            return False
        if REJECT_IF_TRANSFERRED_H_WAS_EJECTED and self.h_is_ejected_by(
            target_h, self.times[frame_i]
        ):
            return False
        if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
            if self.atom_is_in_cap(frame_i, target_h) or self.atom_is_in_cap(frame_i, target_o):
                return False
        return (
            int(self.nearest_o[frame_i, target_h]) == target_o
            and self.nearest_o_distance[frame_i, target_h]
            <= TRANSFERRED_OH_BOND_CUTOFF_A
        )

    def find_anchor(self, transfer_time_fs, target_h, target_o, direction):
        if direction == "before":
            indices = np.where(self.times <= transfer_time_fs + EVENT_TIME_COMPARE_TOL_FS)[0][::-1]
        elif direction == "after":
            indices = np.where(self.times >= transfer_time_fs - EVENT_TIME_COMPARE_TOL_FS)[0]
        else:
            raise ValueError("direction must be before or after")

        for frame_i in indices:
            if MAX_ANCHOR_SEARCH_FS is not None:
                if abs(self.times[frame_i] - transfer_time_fs) > MAX_ANCHOR_SEARCH_FS:
                    continue
            if self.geometry_is_clean_for_side(frame_i, target_h, target_o):
                return int(frame_i)
        return None

    def build_fixed_fragment_ledger(
        self, before_frame_i, after_frame_i, target_h, donor_o, acceptor_o,
    ):
        """Build before/after atom lists with only target H changing parent."""
        before_time = float(self.times[before_frame_i])
        after_time = float(self.times[after_frame_i])

        donor_fixed_h = []
        acceptor_fixed_h = []
        excluded_ejected_before_donor = []
        excluded_ejected_before_acceptor = []

        for h_index_value in self.hydrogen_indices:
            h_index = int(h_index_value)
            if h_index == target_h:
                continue
            if self.h_is_ejected_by(h_index, before_time):
                ejection_parent = self.ejection_nearest_o.get(h_index)
                if ejection_parent == donor_o:
                    excluded_ejected_before_donor.append(h_index)
                elif ejection_parent == acceptor_o:
                    excluded_ejected_before_acceptor.append(h_index)
                continue

            owner_before = int(self.nearest_o[before_frame_i, h_index])
            if owner_before == donor_o:
                donor_fixed_h.append(h_index)
            elif owner_before == acceptor_o:
                acceptor_fixed_h.append(h_index)

        ejected_between_donor = [
            h for h in donor_fixed_h
            if self.ejection_time.get(h) is not None
            and before_time < self.ejection_time[h] <= after_time + EVENT_TIME_COMPARE_TOL_FS
        ]
        ejected_between_acceptor = [
            h for h in acceptor_fixed_h
            if self.ejection_time.get(h) is not None
            and before_time < self.ejection_time[h] <= after_time + EVENT_TIME_COMPARE_TOL_FS
        ]

        return {
            "donor_before_h": sorted(donor_fixed_h + [target_h]),
            "donor_after_h": sorted(donor_fixed_h),
            "acceptor_before_h": sorted(acceptor_fixed_h),
            "acceptor_after_h": sorted(acceptor_fixed_h + [target_h]),
            "excluded_ejected_before_donor": sorted(excluded_ejected_before_donor),
            "excluded_ejected_before_acceptor": sorted(excluded_ejected_before_acceptor),
            "ejected_between_donor": sorted(ejected_between_donor),
            "ejected_between_acceptor": sorted(ejected_between_acceptor),
        }


# ============================================================================
# BOV/DAT density integration
# ============================================================================

def bov_path_for_time(run_dir, target_time_fs):
    raw_index = target_time_fs / DENSITY_OUTPUT_DT_FS
    index = int(round(raw_index))
    actual_time = index * DENSITY_OUTPUT_DT_FS
    if abs(actual_time - target_time_fs) > DENSITY_TIME_MATCH_TOL_FS:
        raise ValueError(
            "Requested %.9f fs is not a density-output time (nearest %.9f fs)"
            % (target_time_fs, actual_time)
        )
    density_dir = run_dir / DENSITY_SUBDIR if DENSITY_SUBDIR else run_dir
    filename = "%s%0*d.bov" % (BOV_PREFIX, BOV_DIGITS, index)
    return density_dir / filename, actual_time, index


def read_bov_header(path):
    params = {}
    first_line = ""
    with path.open("r") as handle:
        for line_i, line in enumerate(handle):
            if line_i == 0:
                first_line = line.strip()
            if ":" in line:
                key, value = line.split(":", 1)
                params[key.strip().upper()] = value.strip()

    required = ["DATA_FILE", "DATA_SIZE", "BRICK_ORIGIN", "BRICK_SIZE"]
    missing = [key for key in required if key not in params]
    if missing:
        raise ValueError("Missing BOV fields %s in %s" % (missing, path))

    data_path = Path(params["DATA_FILE"])
    if not data_path.is_absolute():
        data_path = path.parent / data_path

    grid_size = tuple(int(v) for v in params["DATA_SIZE"].split())
    origin = tuple(float(v) for v in params["BRICK_ORIGIN"].split())
    brick_size = tuple(float(v) for v in params["BRICK_SIZE"].split())
    byte_offset = int(params.get("BYTE_OFFSET", "0"))
    return first_line, data_path, grid_size, origin, brick_size, byte_offset


def integrate_density_file_nearest_atom(
    path, grid_size, byte_offset, origin, brick_size, atom_coords,
):
    """Stream a BOV data file and apply ProcessDensityPointCodyMethod.

    BOV stores x fastest, then y, then z.  Only DENSITY_POINT_CHUNK_SIZE values
    and O(number_of_points_in_chunk) work arrays are resident at once.  This is
    deliberately more conservative than constructing a points-by-atoms array.
    """
    nx, ny, nz = grid_size
    if nx < 2 or ny < 2 or nz < 2:
        raise ValueError("Every BOV grid dimension must be at least 2")
    if DENSITY_POINT_CHUNK_SIZE < 1:
        raise ValueError("DENSITY_POINT_CHUNK_SIZE must be at least 1")

    dx = brick_size[0] / float(nx - 1)
    dy = brick_size[1] / float(ny - 1)
    dz = brick_size[2] / float(nz - 1)
    npoints = nx * ny * nz
    dtype = np.dtype(DENSITY_BINARY_DTYPE)
    required_bytes = byte_offset + dtype.itemsize * npoints
    actual_bytes = path.stat().st_size
    if actual_bytes < required_bytes:
        available_values = max(0, actual_bytes - byte_offset) // dtype.itemsize
        raise ValueError(
            "Expected %d density values in %s, found only %d"
            % (npoints, path, available_values)
        )

    atom_coords = np.asarray(atom_coords, dtype=np.float64)
    atom_density_sums = np.zeros(len(atom_coords), dtype=np.float64)
    cutoff2 = MAX_DENSITY_DISTANCE_A * MAX_DENSITY_DISTANCE_A
    grid_density_sum = 0.0

    with path.open("rb") as handle:
        if byte_offset:
            handle.seek(byte_offset)

        linear_start = 0
        while linear_start < npoints:
            count = min(DENSITY_POINT_CHUNK_SIZE, npoints - linear_start)
            raw = handle.read(dtype.itemsize * count)
            values_native = np.frombuffer(raw, dtype=dtype)
            if values_native.size != count:
                raise ValueError(
                    "Unexpected EOF in %s at density value %d; expected %d more, found %d"
                    % (path, linear_start, count, values_native.size)
                )
            values = values_native.astype(np.float64, copy=False)
            grid_density_sum += float(np.sum(values, dtype=np.float64))

            # Convert flat BOV indices to coordinates without a 3-D meshgrid.
            linear = np.arange(linear_start, linear_start + count, dtype=np.int64)
            ix = linear % nx
            yz = linear // nx
            iy = yz % ny
            iz = yz // ny
            px = origin[0] + ix * dx
            py = origin[1] + iy * dy
            pz = origin[2] + iz * dz

            minimum2 = np.full(count, np.inf, dtype=np.float64)
            nearest = np.full(count, -1, dtype=np.int32)
            for atom_i, atom_pos in enumerate(atom_coords):
                delta_x = px - atom_pos[0]
                delta_y = py - atom_pos[1]
                delta_z = pz - atom_pos[2]
                distance2 = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z

                # <= matches the reference loop: the last atom wins exact ties.
                closer = distance2 <= minimum2
                minimum2[closer] = distance2[closer]
                nearest[closer] = atom_i

            valid = minimum2 <= cutoff2
            if np.any(valid):
                atom_density_sums += np.bincount(
                    nearest[valid], weights=values[valid], minlength=len(atom_coords)
                )
            linear_start += count

    voxel_volume = dx * dy * dz
    per_atom_electrons = atom_density_sums * voxel_volume
    grid_total_electrons = grid_density_sum * voxel_volume
    return per_atom_electrons, grid_total_electrons, voxel_volume


def nuclear_charges_for_species(species):
    missing = sorted(set(symbol for symbol in species if symbol not in NUCLEAR_CHARGE))
    if missing:
        raise ValueError("Missing NUCLEAR_CHARGE entries for: %s" % ", ".join(missing))
    return np.asarray([NUCLEAR_CHARGE[symbol] for symbol in species], dtype=float)


class ChargeResult:
    def __init__(
        self,
        requested_time,
        density_time,
        density_index,
        bov_path,
        bov_first_line,
        xyz_frame_index,
        xyz_time,
        atom_electrons,
        atom_charges,
        assigned_total_electrons,
        grid_total_electrons,
        system_charge,
        voxel_volume,
    ):
        self.requested_time = requested_time
        self.density_time = density_time
        self.density_index = density_index
        self.bov_path = bov_path
        self.bov_first_line = bov_first_line
        self.xyz_frame_index = xyz_frame_index
        self.xyz_time = xyz_time
        self.atom_electrons = atom_electrons
        self.atom_charges = atom_charges
        self.assigned_total_electrons = assigned_total_electrons
        self.grid_total_electrons = grid_total_electrons
        self.system_charge = system_charge
        self.voxel_volume = voxel_volume


def compute_charge_result(run_dir, trajectory, requested_time_fs):
    bov_path, density_time, density_index = bov_path_for_time(run_dir, requested_time_fs)
    if not bov_path.exists():
        raise FileNotFoundError("Missing BOV file: %s" % bov_path)

    xyz_frame_i = trajectory.closest_frame_index(density_time)
    first_line, data_path, grid_size, origin, brick_size, byte_offset = read_bov_header(bov_path)
    if not data_path.exists():
        raise FileNotFoundError("Missing DAT file referenced by %s: %s" % (bov_path, data_path))
    if PRINT_DENSITY_PROGRESS:
        print(
            "    Integrating %s (%d x %d x %d, %d atoms)"
            % (data_path, grid_size[0], grid_size[1], grid_size[2], len(trajectory.species)),
            flush=True,
        )
    atom_electrons, grid_total_electrons, voxel_volume = integrate_density_file_nearest_atom(
        path=data_path,
        grid_size=grid_size,
        byte_offset=byte_offset,
        origin=origin,
        brick_size=brick_size,
        atom_coords=trajectory.coords[xyz_frame_i],
    )
    z_values = nuclear_charges_for_species(trajectory.species)
    atom_charges = z_values - atom_electrons
    assigned_total_electrons = float(np.sum(atom_electrons))
    system_charge = float(np.sum(z_values) - assigned_total_electrons)

    return ChargeResult(
        requested_time=requested_time_fs,
        density_time=density_time,
        density_index=density_index,
        bov_path=bov_path,
        bov_first_line=first_line,
        xyz_frame_index=xyz_frame_i,
        xyz_time=float(trajectory.times[xyz_frame_i]),
        atom_electrons=atom_electrons,
        atom_charges=atom_charges,
        assigned_total_electrons=assigned_total_electrons,
        grid_total_electrons=grid_total_electrons,
        system_charge=system_charge,
        voxel_volume=voxel_volume,
    )


# ============================================================================
# Event analysis and output
# ============================================================================

def load_events(path):
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "run", "dir", "h_index", "old_owner_index", "new_owner_index",
            TRANSFER_TIME_COLUMN,
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError("Event CSV is missing columns: %s" % ", ".join(missing))
        rows = []
        for event_id, row in enumerate(reader, start=1):
            copied = dict(row)
            copied["_event_id"] = event_id
            rows.append(copied)
    return rows


def float_or_none(value):
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def join_floats(values, fmt="%.9f"):
    return ";".join(fmt % float(value) for value in values)


def join_indices_one_based(indices):
    return ";".join(str(int(index) + 1) for index in indices)


def four_frame_sample_times(before_anchor_time, after_anchor_time):
    if FOUR_FRAME_COUNT < 1:
        raise ValueError("FOUR_FRAME_COUNT must be at least 1")
    if FOUR_FRAME_SPACING_FS <= 0.0:
        raise ValueError("FOUR_FRAME_SPACING_FS must be positive")

    before = [
        before_anchor_time - FOUR_FRAME_SPACING_FS * offset
        for offset in range(FOUR_FRAME_COUNT, 0, -1)
    ]
    after = [
        after_anchor_time + FOUR_FRAME_SPACING_FS * offset
        for offset in range(1, FOUR_FRAME_COUNT + 1)
    ]
    return before, after


def make_event_base_row(event):
    return {
        "event_id": event["_event_id"],
        "run": event.get("run", ""),
        "dir": event.get("dir", ""),
        "h_index": event.get("h_index", ""),
        "donor_o_index": event.get("old_owner_index", ""),
        "acceptor_o_index": event.get("new_owner_index", ""),
        "transfer_time_column": TRANSFER_TIME_COLUMN,
        "transfer_time_fs": event.get(TRANSFER_TIME_COLUMN, ""),
    }


EVENT_OUTPUT_FIELDS = [
    "event_id", "run", "dir", "h_index", "donor_o_index", "acceptor_o_index",
    "transfer_time_column", "transfer_time_fs", "status", "reason",
    "before_anchor_time_fs", "after_anchor_time_fs",
    "before_target_oh_distance_A", "after_target_oh_distance_A",
    "anchor_separation_fs",
    "donor_charge_before_e", "donor_charge_after_e",
    "acceptor_charge_before_e", "acceptor_charge_after_e",
    "qH_from_donor_e", "qH_from_acceptor_e", "qH_effective_e",
    "donor_acceptor_mismatch_e", "pair_charge_before_e", "pair_charge_after_e",
    "pair_charge_change_e", "system_charge_before_e", "system_charge_after_e",
    "system_charge_change_e", "donor_acceptor_sanity_pass",
    "system_charge_sanity_pass", "overall_sanity_pass",
    "four_frame_mean_status", "four_frame_mean_reason",
    "four_frame_before_sample_times_fs", "four_frame_after_sample_times_fs",
    "donor_charge_before_four_frame_mean_e",
    "donor_charge_after_four_frame_mean_e",
    "acceptor_charge_before_four_frame_mean_e",
    "acceptor_charge_after_four_frame_mean_e",
    "donor_charge_before_four_frame_std_e",
    "donor_charge_after_four_frame_std_e",
    "acceptor_charge_before_four_frame_std_e",
    "acceptor_charge_after_four_frame_std_e",
    "qH_from_donor_four_frame_mean_e",
    "qH_from_acceptor_four_frame_mean_e",
    "qH_effective_four_frame_mean_e",
    "qH_effective_four_frame_std_e",
    "donor_acceptor_mismatch_four_frame_mean_e",
    "pair_charge_before_four_frame_mean_e",
    "pair_charge_after_four_frame_mean_e",
    "pair_charge_change_four_frame_mean_e",
    "system_charge_before_four_frame_mean_e",
    "system_charge_after_four_frame_mean_e",
    "system_charge_change_four_frame_mean_e",
    "donor_acceptor_four_frame_sanity_pass",
    "system_charge_four_frame_sanity_pass",
    "overall_four_frame_sanity_pass",
    "qH_effective_four_frame_minus_nearest_e",
    "qH_effective_four_frame_vs_nearest_abs_difference_e",
    "donor_H_before", "donor_H_after", "acceptor_H_before", "acceptor_H_after",
    "donor_H_ejected_before_excluded_both",
    "acceptor_H_ejected_before_excluded_both",
    "donor_H_ejected_between_kept_both",
    "acceptor_H_ejected_between_kept_both",
    "target_H_ejection_type", "target_H_ejection_time_fs",
]


FRAME_OUTPUT_FIELDS = [
    "event_id", "run", "dir", "sample_set", "sample_number", "side",
    "requested_time_fs",
    "density_time_fs", "xyz_time_fs", "density_index", "bov_file",
    "target_oh_distance_A",
    "donor_charge_e", "acceptor_charge_e", "pair_charge_e", "system_charge_e",
    "assigned_total_electrons", "grid_total_electrons", "unassigned_grid_electrons",
    "donor_H_members", "acceptor_H_members",
    "donor_H_ejected_between_kept", "acceptor_H_ejected_between_kept",
]


EJECTION_OUTPUT_FIELDS = [
    "run", "dir", "h_index", "ejection_type", "ejection_time_fs",
    "ejection_frame_index_zero_based", "nearest_o_index", "nearest_oh_distance_A",
    "h_in_cap", "nearest_o_in_cap", "tracking_stop_time_fs",
    "tracking_stop_reason",
]


def reject_result(event, reason, status="REJECTED"):
    row = make_event_base_row(event)
    row["status"] = status
    row["reason"] = reason
    return row


def other_transfer_in_window(event, run_events, window_start, window_end):
    for other in run_events:
        if other["_event_id"] == event["_event_id"]:
            continue
        other_time = float_or_none(other.get(TRANSFER_TIME_COLUMN))
        if other_time is None:
            continue
        if (
            window_start - EVENT_TIME_COMPARE_TOL_FS
            <= other_time
            <= window_end + EVENT_TIME_COMPARE_TOL_FS
        ):
            return other
    return None


def fragment_charge_at_frame(charge_result, donor_o, acceptor_o, donor_h, acceptor_h):
    donor_atoms = [donor_o] + list(donor_h)
    acceptor_atoms = [acceptor_o] + list(acceptor_h)
    donor_charge = float(np.sum(charge_result.atom_charges[donor_atoms]))
    acceptor_charge = float(np.sum(charge_result.atom_charges[acceptor_atoms]))
    return {
        "donor_charge": donor_charge,
        "acceptor_charge": acceptor_charge,
        "donor_h": list(donor_h),
        "acceptor_h": list(acceptor_h),
    }


def validate_four_frame_geometry(
    trajectory, sample_times, target_h, target_o, donor_o, acceptor_o,
):
    """Validate that every averaging sample represents the intended side."""
    frame_indices = []
    target_atoms = (target_h, donor_o, acceptor_o)
    for sample_time in sample_times:
        if sample_time + EVENT_TIME_COMPARE_TOL_FS < MIN_ANALYSIS_TIME_FS:
            raise ValueError(
                "four-frame sample %.6f fs is before MIN_ANALYSIS_TIME_FS %.6f fs"
                % (sample_time, MIN_ANALYSIS_TIME_FS)
            )
        frame_i = trajectory.closest_frame_index(sample_time)
        if REQUIRE_TARGET_OH_GEOMETRY_AT_FOUR_FRAME_SAMPLES:
            if not trajectory.geometry_is_clean_for_side(frame_i, target_h, target_o):
                raise ValueError(
                    "target O-H geometry is not valid at four-frame sample %.6f fs"
                    % sample_time
                )
        if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
            in_cap_atoms = [
                atom_i for atom_i in target_atoms
                if trajectory.atom_is_in_cap(frame_i, atom_i)
            ]
            if in_cap_atoms:
                raise ValueError(
                    "four-frame sample %.6f fs has target atom(s) in CAP: %s"
                    % (sample_time, join_indices_one_based(in_cap_atoms))
                )
        frame_indices.append(frame_i)
    return frame_indices


def calculate_charge_sample(
    event, run_dir, trajectory, charge_cache, ledger, donor_o, acceptor_o,
    side, sample_set, sample_number, sample_time, target_h,
):
    """Calculate one charge sample and its detailed output row."""
    donor_h = ledger["donor_%s_h" % side]
    acceptor_h = ledger["acceptor_%s_h" % side]
    _, _, density_index = bov_path_for_time(run_dir, sample_time)
    cache_key = (str(run_dir.resolve()), density_index)
    if cache_key not in charge_cache:
        charge_cache[cache_key] = compute_charge_result(
            run_dir, trajectory, sample_time
        )
    charge_result = charge_cache[cache_key]
    fragment = fragment_charge_at_frame(
        charge_result, donor_o, acceptor_o, donor_h, acceptor_h
    )
    xyz_frame_i = charge_result.xyz_frame_index
    target_distance = float(trajectory.nearest_o_distance[xyz_frame_i, target_h])

    frame_row = {
        "event_id": event["_event_id"],
        "run": event.get("run", ""),
        "dir": event.get("dir", ""),
        "sample_set": sample_set,
        "sample_number": sample_number,
        "side": side,
        "requested_time_fs": "%.9f" % sample_time,
        "density_time_fs": "%.9f" % charge_result.density_time,
        "xyz_time_fs": "%.9f" % charge_result.xyz_time,
        "density_index": charge_result.density_index,
        "bov_file": charge_result.bov_path.name,
        "target_oh_distance_A": "%.9f" % target_distance,
        "donor_charge_e": "%.12f" % fragment["donor_charge"],
        "acceptor_charge_e": "%.12f" % fragment["acceptor_charge"],
        "pair_charge_e": "%.12f" % (
            fragment["donor_charge"] + fragment["acceptor_charge"]
        ),
        "system_charge_e": "%.12f" % charge_result.system_charge,
        "assigned_total_electrons": "%.12f" % charge_result.assigned_total_electrons,
        "grid_total_electrons": "%.12f" % charge_result.grid_total_electrons,
        "unassigned_grid_electrons": "%.12f" % (
            charge_result.grid_total_electrons
            - charge_result.assigned_total_electrons
        ),
        "donor_H_members": join_indices_one_based(fragment["donor_h"]),
        "acceptor_H_members": join_indices_one_based(fragment["acceptor_h"]),
        "donor_H_ejected_between_kept": join_indices_one_based(
            ledger["ejected_between_donor"]
        ),
        "acceptor_H_ejected_between_kept": join_indices_one_based(
            ledger["ejected_between_acceptor"]
        ),
    }
    return {
        "fragment": fragment,
        "charge_result": charge_result,
        "frame_row": frame_row,
    }


def summarize_two_sides(before_samples, after_samples):
    """Summarize either one nearest sample or equal-sized averaging sets."""
    if not before_samples or not after_samples:
        raise ValueError("both before and after samples are required")
    if len(before_samples) != len(after_samples):
        raise ValueError("before and after sample counts differ")

    donor_before_values = np.asarray([
        sample["fragment"]["donor_charge"] for sample in before_samples
    ], dtype=float)
    donor_after_values = np.asarray([
        sample["fragment"]["donor_charge"] for sample in after_samples
    ], dtype=float)
    acceptor_before_values = np.asarray([
        sample["fragment"]["acceptor_charge"] for sample in before_samples
    ], dtype=float)
    acceptor_after_values = np.asarray([
        sample["fragment"]["acceptor_charge"] for sample in after_samples
    ], dtype=float)
    system_before_values = np.asarray([
        sample["charge_result"].system_charge for sample in before_samples
    ], dtype=float)
    system_after_values = np.asarray([
        sample["charge_result"].system_charge for sample in after_samples
    ], dtype=float)

    donor_before = float(np.mean(donor_before_values))
    donor_after = float(np.mean(donor_after_values))
    acceptor_before = float(np.mean(acceptor_before_values))
    acceptor_after = float(np.mean(acceptor_after_values))
    system_before = float(np.mean(system_before_values))
    system_after = float(np.mean(system_after_values))

    qh_donor = donor_before - donor_after
    qh_acceptor = acceptor_after - acceptor_before
    qh_effective = 0.5 * (qh_donor + qh_acceptor)
    paired_qh_effective = 0.5 * (
        (donor_before_values - donor_after_values)
        + (acceptor_after_values - acceptor_before_values)
    )
    mismatch = abs(qh_donor - qh_acceptor)
    pair_before = donor_before + acceptor_before
    pair_after = donor_after + acceptor_after
    pair_change = pair_after - pair_before
    system_change = system_after - system_before
    da_pass = mismatch <= DONOR_ACCEPTOR_MISMATCH_TOL_E
    system_pass = abs(system_change) <= SYSTEM_CHARGE_CHANGE_TOL_E

    return {
        "donor_before": donor_before,
        "donor_after": donor_after,
        "acceptor_before": acceptor_before,
        "acceptor_after": acceptor_after,
        "donor_before_std": float(np.std(donor_before_values)),
        "donor_after_std": float(np.std(donor_after_values)),
        "acceptor_before_std": float(np.std(acceptor_before_values)),
        "acceptor_after_std": float(np.std(acceptor_after_values)),
        "qh_donor": qh_donor,
        "qh_acceptor": qh_acceptor,
        "qh_effective": qh_effective,
        "qh_effective_std": float(np.std(paired_qh_effective)),
        "mismatch": mismatch,
        "pair_before": pair_before,
        "pair_after": pair_after,
        "pair_change": pair_change,
        "system_before": system_before,
        "system_after": system_after,
        "system_change": system_change,
        "da_pass": da_pass,
        "system_pass": system_pass,
        "overall_pass": da_pass and system_pass,
    }


def analyze_event(event, run_events, run_dir, trajectory, charge_cache):
    transfer_time = float_or_none(event.get(TRANSFER_TIME_COLUMN))
    if transfer_time is None:
        return reject_result(event, "missing transfer time"), []

    try:
        target_h = int(event["h_index"]) - 1
        donor_o = int(event["old_owner_index"]) - 1
        acceptor_o = int(event["new_owner_index"]) - 1
    except Exception:
        return reject_result(event, "invalid atom index in event CSV"), []

    natoms = len(trajectory.species)
    if not (0 <= target_h < natoms and 0 <= donor_o < natoms and 0 <= acceptor_o < natoms):
        return reject_result(event, "event atom index outside trajectory"), []
    if trajectory.species[target_h] != "H":
        return reject_result(event, "h_index does not identify H"), []
    if trajectory.species[donor_o] != "O" or trajectory.species[acceptor_o] != "O":
        return reject_result(event, "donor/acceptor index does not identify O"), []
    if donor_o == acceptor_o:
        return reject_result(event, "donor and acceptor O are identical"), []

    before_anchor_i = trajectory.find_anchor(
        transfer_time, target_h, donor_o, direction="before"
    )
    if before_anchor_i is None:
        return reject_result(event, "no clean pre-transfer geometry anchor"), []
    after_anchor_i = trajectory.find_anchor(
        transfer_time, target_h, acceptor_o, direction="after"
    )
    if after_anchor_i is None:
        return reject_result(event, "no clean post-transfer geometry anchor"), []

    before_anchor_time = float(trajectory.times[before_anchor_i])
    after_anchor_time = float(trajectory.times[after_anchor_i])
    if after_anchor_time <= before_anchor_time + EVENT_TIME_COMPARE_TOL_FS:
        return reject_result(event, "post-transfer anchor is not later than pre-transfer anchor"), []
    before_target_distance = float(
        trajectory.nearest_o_distance[before_anchor_i, target_h]
    )
    after_target_distance = float(
        trajectory.nearest_o_distance[after_anchor_i, target_h]
    )

    if REJECT_TARGET_ATOMS_IN_CAP_AT_ANCHORS:
        target_atoms = (target_h, donor_o, acceptor_o)
        for side, frame_i in (("before", before_anchor_i), ("after", after_anchor_i)):
            in_cap_atoms = [atom_i for atom_i in target_atoms if trajectory.atom_is_in_cap(frame_i, atom_i)]
            if in_cap_atoms:
                return reject_result(
                    event,
                    "%s anchor has target atom(s) in CAP: %s"
                    % (side, join_indices_one_based(in_cap_atoms)),
                ), []

    analysis_start = before_anchor_time
    analysis_end = after_anchor_time
    if REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW:
        conflict = other_transfer_in_window(
            event, run_events, analysis_start, analysis_end
        )
        if conflict is not None:
            reason = "other H transfer event %s occurs inside analysis window" % conflict["_event_id"]
            return reject_result(event, reason), []

    ledger = trajectory.build_fixed_fragment_ledger(
        before_anchor_i, after_anchor_i, target_h, donor_o, acceptor_o
    )

    # Nearest-anchor estimate.  Existing output names continue to refer to this
    # estimate so that older downstream analysis remains compatible.
    frame_rows = []
    nearest_samples = {"before": [], "after": []}
    try:
        for side, sample_time in (
            ("before", before_anchor_time),
            ("after", after_anchor_time),
        ):
            sample = calculate_charge_sample(
                event, run_dir, trajectory, charge_cache, ledger,
                donor_o, acceptor_o, side, "nearest_anchor", 1,
                sample_time, target_h,
            )
            nearest_samples[side].append(sample)
            frame_rows.append(sample["frame_row"])
        nearest = summarize_two_sides(
            nearest_samples["before"], nearest_samples["after"]
        )
    except Exception as exc:
        return reject_result(
            event, "nearest-frame charge calculation failed: %s" % exc,
            status="ERROR",
        ), []

    # Four-frame estimate.  The same fixed fragment ledger and ejection rules
    # are used; only the density sampling times differ from the nearest result.
    four_before_times, four_after_times = four_frame_sample_times(
        before_anchor_time, after_anchor_time
    )
    four_frame_status = "ANALYZED"
    four_frame_reason = ""
    four_frame = None
    try:
        validate_four_frame_geometry(
            trajectory, four_before_times, target_h, donor_o, donor_o, acceptor_o
        )
        validate_four_frame_geometry(
            trajectory, four_after_times, target_h, acceptor_o, donor_o, acceptor_o
        )

        if REJECT_OTHER_H_TRANSFER_IN_ANALYSIS_WINDOW:
            conflict = other_transfer_in_window(
                event, run_events, min(four_before_times), max(four_after_times)
            )
            if conflict is not None:
                raise ValueError(
                    "other H transfer event %s occurs inside the four-frame window"
                    % conflict["_event_id"]
                )

        four_samples = {"before": [], "after": []}
        four_detail_rows = []
        for side, sample_times in (
            ("before", four_before_times),
            ("after", four_after_times),
        ):
            for sample_number, sample_time in enumerate(sample_times, start=1):
                sample = calculate_charge_sample(
                    event, run_dir, trajectory, charge_cache, ledger,
                    donor_o, acceptor_o, side,
                    "four_frame_average_member", sample_number,
                    sample_time, target_h,
                )
                four_samples[side].append(sample)
                four_detail_rows.append(sample["frame_row"])
        four_frame = summarize_two_sides(
            four_samples["before"], four_samples["after"]
        )
        frame_rows.extend(four_detail_rows)
    except Exception as exc:
        four_frame_status = "UNAVAILABLE"
        four_frame_reason = str(exc)

    row = make_event_base_row(event)
    row.update({
        "status": "ANALYZED",
        "reason": "",
        "before_anchor_time_fs": "%.9f" % before_anchor_time,
        "after_anchor_time_fs": "%.9f" % after_anchor_time,
        "before_target_oh_distance_A": "%.9f" % before_target_distance,
        "after_target_oh_distance_A": "%.9f" % after_target_distance,
        "anchor_separation_fs": "%.9f" % (after_anchor_time - before_anchor_time),
        "donor_charge_before_e": "%.12f" % nearest["donor_before"],
        "donor_charge_after_e": "%.12f" % nearest["donor_after"],
        "acceptor_charge_before_e": "%.12f" % nearest["acceptor_before"],
        "acceptor_charge_after_e": "%.12f" % nearest["acceptor_after"],
        "qH_from_donor_e": "%.12f" % nearest["qh_donor"],
        "qH_from_acceptor_e": "%.12f" % nearest["qh_acceptor"],
        "qH_effective_e": "%.12f" % nearest["qh_effective"],
        "donor_acceptor_mismatch_e": "%.12f" % nearest["mismatch"],
        "pair_charge_before_e": "%.12f" % nearest["pair_before"],
        "pair_charge_after_e": "%.12f" % nearest["pair_after"],
        "pair_charge_change_e": "%.12f" % nearest["pair_change"],
        "system_charge_before_e": "%.12f" % nearest["system_before"],
        "system_charge_after_e": "%.12f" % nearest["system_after"],
        "system_charge_change_e": "%.12f" % nearest["system_change"],
        "donor_acceptor_sanity_pass": int(nearest["da_pass"]),
        "system_charge_sanity_pass": int(nearest["system_pass"]),
        "overall_sanity_pass": int(nearest["overall_pass"]),
        "four_frame_mean_status": four_frame_status,
        "four_frame_mean_reason": four_frame_reason,
        "four_frame_before_sample_times_fs": join_floats(four_before_times),
        "four_frame_after_sample_times_fs": join_floats(four_after_times),
        "donor_H_before": join_indices_one_based(ledger["donor_before_h"]),
        "donor_H_after": join_indices_one_based(ledger["donor_after_h"]),
        "acceptor_H_before": join_indices_one_based(ledger["acceptor_before_h"]),
        "acceptor_H_after": join_indices_one_based(ledger["acceptor_after_h"]),
        "donor_H_ejected_before_excluded_both": join_indices_one_based(
            ledger["excluded_ejected_before_donor"]
        ),
        "acceptor_H_ejected_before_excluded_both": join_indices_one_based(
            ledger["excluded_ejected_before_acceptor"]
        ),
        "donor_H_ejected_between_kept_both": join_indices_one_based(
            ledger["ejected_between_donor"]
        ),
        "acceptor_H_ejected_between_kept_both": join_indices_one_based(
            ledger["ejected_between_acceptor"]
        ),
        "target_H_ejection_type": trajectory.ejection_type.get(target_h) or "",
        "target_H_ejection_time_fs": (
            "%.9f" % trajectory.ejection_time[target_h]
            if trajectory.ejection_time.get(target_h) is not None else ""
        ),
    })

    if four_frame is not None:
        method_difference = four_frame["qh_effective"] - nearest["qh_effective"]
        row.update({
            "donor_charge_before_four_frame_mean_e": "%.12f" % four_frame["donor_before"],
            "donor_charge_after_four_frame_mean_e": "%.12f" % four_frame["donor_after"],
            "acceptor_charge_before_four_frame_mean_e": "%.12f" % four_frame["acceptor_before"],
            "acceptor_charge_after_four_frame_mean_e": "%.12f" % four_frame["acceptor_after"],
            "donor_charge_before_four_frame_std_e": "%.12f" % four_frame["donor_before_std"],
            "donor_charge_after_four_frame_std_e": "%.12f" % four_frame["donor_after_std"],
            "acceptor_charge_before_four_frame_std_e": "%.12f" % four_frame["acceptor_before_std"],
            "acceptor_charge_after_four_frame_std_e": "%.12f" % four_frame["acceptor_after_std"],
            "qH_from_donor_four_frame_mean_e": "%.12f" % four_frame["qh_donor"],
            "qH_from_acceptor_four_frame_mean_e": "%.12f" % four_frame["qh_acceptor"],
            "qH_effective_four_frame_mean_e": "%.12f" % four_frame["qh_effective"],
            "qH_effective_four_frame_std_e": "%.12f" % four_frame["qh_effective_std"],
            "donor_acceptor_mismatch_four_frame_mean_e": "%.12f" % four_frame["mismatch"],
            "pair_charge_before_four_frame_mean_e": "%.12f" % four_frame["pair_before"],
            "pair_charge_after_four_frame_mean_e": "%.12f" % four_frame["pair_after"],
            "pair_charge_change_four_frame_mean_e": "%.12f" % four_frame["pair_change"],
            "system_charge_before_four_frame_mean_e": "%.12f" % four_frame["system_before"],
            "system_charge_after_four_frame_mean_e": "%.12f" % four_frame["system_after"],
            "system_charge_change_four_frame_mean_e": "%.12f" % four_frame["system_change"],
            "donor_acceptor_four_frame_sanity_pass": int(four_frame["da_pass"]),
            "system_charge_four_frame_sanity_pass": int(four_frame["system_pass"]),
            "overall_four_frame_sanity_pass": int(four_frame["overall_pass"]),
            "qH_effective_four_frame_minus_nearest_e": "%.12f" % method_difference,
            "qH_effective_four_frame_vs_nearest_abs_difference_e": "%.12f" % abs(method_difference),
        })
    return row, frame_rows


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ejection_rows_for_trajectory(trajectory, run, directory_name):
    rows = []
    for h_index_value in trajectory.hydrogen_indices:
        h_index = int(h_index_value)
        if trajectory.ejection_time.get(h_index) is None:
            continue
        rows.append({
            "run": run,
            "dir": directory_name,
            "h_index": h_index + 1,
            "ejection_type": trajectory.ejection_type[h_index],
            "ejection_time_fs": "%.9f" % trajectory.ejection_time[h_index],
            "ejection_frame_index_zero_based": trajectory.ejection_frame_index[h_index],
            "nearest_o_index": trajectory.ejection_nearest_o[h_index] + 1,
            "nearest_oh_distance_A": "%.9f" % trajectory.ejection_distance[h_index],
            "h_in_cap": int(trajectory.ejection_h_in_cap[h_index]),
            "nearest_o_in_cap": int(trajectory.ejection_o_in_cap[h_index]),
            "tracking_stop_time_fs": "%.9f" % trajectory.ejection_tracking_stop_time[h_index],
            "tracking_stop_reason": trajectory.ejection_tracking_stop_reason[h_index],
        })
    return rows


def write_summary(path, event_rows, total_input_events, trajectories_loaded, cache_size):
    status_counts = Counter(row.get("status", "") for row in event_rows)
    analyzed = [row for row in event_rows if row.get("status") == "ANALYZED"]
    sanity_pass = sum(int(row.get("overall_sanity_pass", 0)) for row in analyzed)

    qh_values = [float(row["qH_effective_e"]) for row in analyzed]
    four_frame_analyzed = [
        row for row in analyzed
        if row.get("four_frame_mean_status") == "ANALYZED"
        and row.get("qH_effective_four_frame_mean_e", "") != ""
    ]
    four_frame_qh_values = [
        float(row["qH_effective_four_frame_mean_e"])
        for row in four_frame_analyzed
    ]
    method_abs_differences = [
        float(row["qH_effective_four_frame_vs_nearest_abs_difference_e"])
        for row in four_frame_analyzed
    ]
    summary = {
        "input_events": total_input_events,
        "trajectories_loaded": trajectories_loaded,
        "analyzed_events": len(analyzed),
        "rejected_events": status_counts.get("REJECTED", 0),
        "error_events": status_counts.get("ERROR", 0),
        "sanity_pass_events": sanity_pass,
        "sanity_pass_fraction": (sanity_pass / len(analyzed)) if analyzed else "",
        "mean_qH_effective_e": float(np.mean(qh_values)) if qh_values else "",
        "std_qH_effective_e": float(np.std(qh_values)) if qh_values else "",
        "min_qH_effective_e": min(qh_values) if qh_values else "",
        "max_qH_effective_e": max(qh_values) if qh_values else "",
        "four_frame_analyzed_events": len(four_frame_analyzed),
        "four_frame_unavailable_events": len(analyzed) - len(four_frame_analyzed),
        "mean_qH_effective_four_frame_mean_e": (
            float(np.mean(four_frame_qh_values)) if four_frame_qh_values else ""
        ),
        "std_qH_effective_four_frame_mean_e": (
            float(np.std(four_frame_qh_values)) if four_frame_qh_values else ""
        ),
        "mean_abs_four_frame_vs_nearest_difference_e": (
            float(np.mean(method_abs_differences)) if method_abs_differences else ""
        ),
        "max_abs_four_frame_vs_nearest_difference_e": (
            max(method_abs_differences) if method_abs_differences else ""
        ),
        "unique_density_frames_integrated": cache_size,
    }
    write_csv(path, list(summary.keys()), [summary])
    return summary, status_counts


def main():
    base_dir = BASE_DIR.resolve()
    event_csv = base_dir / EVENTS_CSV_NAME
    if not event_csv.exists():
        raise FileNotFoundError("Event CSV not found: %s" % event_csv)

    events = load_events(event_csv)
    events_by_dir = defaultdict(list)
    for event in events:
        events_by_dir[event.get("dir", "")].append(event)

    event_rows = []
    frame_rows = []
    ejection_rows = []
    charge_cache = {}
    trajectories_loaded = 0

    print("=== Clean H-transfer charge analysis ===", flush=True)
    print("Event time column:", TRANSFER_TIME_COLUMN, flush=True)
    print("Input events:", len(events), flush=True)
    print("Density output interval [fs]:", DENSITY_OUTPUT_DT_FS, flush=True)

    for directory_name, run_events in events_by_dir.items():
        if not directory_name:
            for event in run_events:
                event_rows.append(reject_result(event, "missing run directory in event CSV", status="ERROR"))
            continue

        run_dir = base_dir / directory_name
        xyz_path = run_dir / XYZ_NAME
        if not xyz_path.exists():
            for event in run_events:
                event_rows.append(reject_result(event, "missing trajectory: %s" % xyz_path, status="ERROR"))
            continue

        try:
            trajectory = read_xyz_trajectory(xyz_path)
            trajectories_loaded += 1
            ejection_rows.extend(ejection_rows_for_trajectory(
                trajectory=trajectory,
                run=run_events[0].get("run", "") if run_events else "",
                directory_name=directory_name,
            ))
        except Exception as exc:
            for event in run_events:
                event_rows.append(reject_result(event, "trajectory read failed: %s" % exc, status="ERROR"))
            continue

        print("Processing %s: %d event(s)" % (directory_name, len(run_events)), flush=True)
        for event in run_events:
            result_row, result_frame_rows = analyze_event(
                event=event,
                run_events=run_events,
                run_dir=run_dir,
an atom is assigned to the nearest atom.  Atomic charges are Z - N_e.
                trajectory=trajectory,
                charge_cache=charge_cache,
            )
            event_rows.append(result_row)
            frame_rows.extend(result_frame_rows)
            print(
                "  event %s H%s O%s->O%s: %s%s"
                % (
                    event["_event_id"], event.get("h_index", "?"),
                    event.get("old_owner_index", "?"), event.get("new_owner_index", "?"),
                    result_row.get("status", ""),
                    (" (" + result_row.get("reason", "") + ")") if result_row.get("reason") else "",
                ),
                flush=True,
            )
    out_events = base_dir / OUT_EVENTS_CSV
    out_frames = base_dir / OUT_FRAMES_CSV
    out_summary = base_dir / OUT_SUMMARY_CSV
    out_ejections = base_dir / OUT_EJECTIONS_CSV
    write_csv(out_events, EVENT_OUTPUT_FIELDS, event_rows)
    write_csv(out_frames, FRAME_OUTPUT_FIELDS, frame_rows)
    write_csv(out_ejections, EJECTION_OUTPUT_FIELDS, ejection_rows)
    summary, status_counts = write_summary(
        out_summary, event_rows, len(events), trajectories_loaded, len(charge_cache)
    )
    print("=== Overall ===", flush=True)
    for status, count in sorted(status_counts.items()):
        print("%-10s %d" % (status + ":", count), flush=True)
    print("Sanity-pass analyzed events: %s/%s" % (
        summary["sanity_pass_events"], summary["analyzed_events"]
    ), flush=True)
    print("Wrote:", out_events, flush=True)
    print("Wrote:", out_frames, flush=True)
    print("Wrote:", out_summary, flush=True)
    print("Wrote:", out_ejections, flush=True)
if __name__ == "__main__":
    main()
