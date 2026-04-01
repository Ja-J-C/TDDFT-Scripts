#!/usr/bin/env python3
import re
import math
import csv
from pathlib import Path
from collections import Counter

# =============================
# User-tunable parameters
# =============================
# Geometry thresholds (Angstrom)
R_ASSIGN = 1.5   # assign/re-assign an H to an O only if nearest O is within this distance
R_ENTER = 1.75    # new acceptor must stay within this distance
R_LEAVE =  2   # old donor must be farther than this distance

# Time thresholds (fs)
T_STABLE_FS = 10.0     # continuous residence on new O before it can count
T_BACKCHECK_FS = 5.0   # reject if it quickly returns to old O within this extra window

# Input / output
DT_FS_PER_ITER = 0.001
XYZ_NAME = "trajectory.xyz"
OUTCSV_RUN = "h_rearrangement_by_run.csv"
OUTCSV_EVENTS = "h_rearrangement_events.csv"
OUTCSV_SUM = "h_rearrangement_overall.csv"

TIME_RE = re.compile(r"time\[fs\]\s*=\s*([0-9Ee+\-\.]+)")
ITER_RE = re.compile(r"iter\s*=\s*([0-9]+)")


def dist(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def parse_time_fs(comment_line):
    m = TIME_RE.search(comment_line)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def parse_iter(comment_line):
    m = ITER_RE.search(comment_line)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def iter_xyz_frames(xyz_path):
    with xyz_path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                return

            line = line.strip()
            if not line:
                continue

            nat = int(line)

            comment = f.readline()
            if not comment:
                raise RuntimeError("EOF after natoms")

            time_fs = parse_time_fs(comment)
            iter_idx = parse_iter(comment)
            if time_fs is None and iter_idx is not None:
                time_fs = iter_idx * DT_FS_PER_ITER

            species = []
            coords = []
            for _ in range(nat):
                atom_line = f.readline()
                if not atom_line:
                    raise RuntimeError("EOF inside frame")

                parts = atom_line.split()
                if len(parts) < 4:
                    raise RuntimeError("Bad atom line: " + atom_line.strip())

                species.append(parts[0])
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

            yield time_fs, iter_idx, species, coords


def oxygen_distances_for_h(h_idx, coords, O_idx):
    hpos = coords[h_idx]
    pairs = []
    for oi in O_idx:
        pairs.append((oi, dist(hpos, coords[oi])))
    pairs.sort(key=lambda x: x[1])
    return pairs


def maybe_assign_owner(state, nearest_o, nearest_d, time_fs):
    if nearest_o is not None and nearest_d <= R_ASSIGN:
        state["owner"] = nearest_o
        state["owner_since_fs"] = time_fs
    else:
        state["owner"] = None
        state["owner_since_fs"] = None


def reset_candidate(state):
    state["cand_old_owner"] = None
    state["cand_new_owner"] = None
    state["cand_start_fs"] = None
    state["cand_start_iter"] = None
    state["backcheck_start_fs"] = None
    state["backcheck_start_iter"] = None
    state["in_backcheck"] = False


def start_candidate(state, old_owner, new_owner, time_fs, iter_idx):
    state["cand_old_owner"] = old_owner
    state["cand_new_owner"] = new_owner
    state["cand_start_fs"] = time_fs
    state["cand_start_iter"] = iter_idx
    state["backcheck_start_fs"] = None
    state["backcheck_start_iter"] = None
    state["in_backcheck"] = False


def elapsed_fs(t_now, t_start):
    if t_now is None or t_start is None:
        return None
    return t_now - t_start


def count_rearrangement_one_xyz(xyz_path):
    it = iter_xyz_frames(xyz_path)

    try:
        t0, it0, sp0, co0 = next(it)
    except StopIteration:
        return {
            "status": "EMPTY",
            "n_O": 0,
            "n_H": 0,
            "n_transfer_events": 0,
            "n_H_with_transfer": 0,
        }, []

    O_idx = [i for i, s in enumerate(sp0) if s.upper() == "O"]
    H_idx = [i for i, s in enumerate(sp0) if s.upper() == "H"]

    if not O_idx or not H_idx:
        return {
            "status": "NO_O_OR_H",
            "n_O": len(O_idx),
            "n_H": len(H_idx),
            "n_transfer_events": 0,
            "n_H_with_transfer": 0,
        }, []

    states = {}
    events = []
    H_with_transfer = set()

    for h in H_idx:
        dlist = oxygen_distances_for_h(h, co0, O_idx)
        nearest_o, nearest_d = dlist[0] if dlist else (None, 1.0e99)
        states[h] = {
            "owner": None,
            "owner_since_fs": None,
            "cand_old_owner": None,
            "cand_new_owner": None,
            "cand_start_fs": None,
            "cand_start_iter": None,
            "backcheck_start_fs": None,
            "backcheck_start_iter": None,
            "in_backcheck": False,
        }
        maybe_assign_owner(states[h], nearest_o, nearest_d, t0)

    def process_frame(time_fs, iter_idx, coords):
        for h in H_idx:
            state = states[h]

            dlist = oxygen_distances_for_h(h, coords, O_idx)
            if not dlist:
                continue

            hpos = coords[h]
            nearest_o, nearest_d = dlist[0]
            owner = state["owner"]

            # If no owner has been assigned yet, try to assign one and move on.
            if owner is None:
                maybe_assign_owner(state, nearest_o, nearest_d, time_fs)
                continue

            # -------- Backcheck stage --------
            if state["in_backcheck"]:
                new_owner = state["cand_new_owner"]
                old_owner = state["cand_old_owner"]
                new_d = dist(hpos, coords[new_owner])
                old_d = dist(hpos, coords[old_owner])

                # Immediate return to old owner => reject
                if nearest_o == old_owner and old_d <= R_ENTER:
                    state["owner"] = old_owner
                    state["owner_since_fs"] = time_fs
                    reset_candidate(state)
                    continue

                # Loss of stable new ownership during backcheck => reject
                if not (nearest_o == new_owner and new_d <= R_ENTER and old_d >= R_LEAVE):
                    maybe_assign_owner(state, nearest_o, nearest_d, time_fs)
                    reset_candidate(state)
                    continue

                dt_back = elapsed_fs(time_fs, state["backcheck_start_fs"])
                if dt_back is not None and dt_back >= T_BACKCHECK_FS:
                    events.append({
                        "h_index": h + 1,
                        "old_owner_index": old_owner + 1,
                        "new_owner_index": new_owner + 1,
                        "start_time_fs": state["cand_start_fs"],
                        "start_iter_idx": state["cand_start_iter"],
                        "stable_time_fs": state["backcheck_start_fs"],
                        "stable_iter_idx": state["backcheck_start_iter"],
                        "confirm_time_fs": time_fs,
                        "confirm_iter_idx": iter_idx,
                        "oh_new_A": new_d,
                        "oh_old_A": old_d,
                    })
                    H_with_transfer.add(h)
                    state["owner"] = new_owner
                    state["owner_since_fs"] = state["cand_start_fs"]
                    reset_candidate(state)
                continue

            # -------- Candidate stage (not yet stable enough) --------
            if state["cand_new_owner"] is not None:
                new_owner = state["cand_new_owner"]
                old_owner = state["cand_old_owner"]
                new_d = dist(hpos, coords[new_owner])
                old_d = dist(hpos, coords[old_owner])

                if nearest_o == new_owner and new_d <= R_ENTER and old_d >= R_LEAVE:
                    dt_cand = elapsed_fs(time_fs, state["cand_start_fs"])
                    if dt_cand is not None and dt_cand >= T_STABLE_FS:
                        if T_BACKCHECK_FS > 0.0:
                            state["in_backcheck"] = True
                            state["backcheck_start_fs"] = time_fs
                            state["backcheck_start_iter"] = iter_idx
                        else:
                            events.append({
                                "h_index": h + 1,
                                "old_owner_index": old_owner + 1,
                                "new_owner_index": new_owner + 1,
                                "start_time_fs": state["cand_start_fs"],
                                "start_iter_idx": state["cand_start_iter"],
                                "stable_time_fs": time_fs,
                                "stable_iter_idx": iter_idx,
                                "confirm_time_fs": time_fs,
                                "confirm_iter_idx": iter_idx,
                                "oh_new_A": new_d,
                                "oh_old_A": old_d,
                            })
                            H_with_transfer.add(h)
                            state["owner"] = new_owner
                            state["owner_since_fs"] = state["cand_start_fs"]
                            reset_candidate(state)
                    continue
                else:
                    # Candidate broke before becoming stable enough.
                    reset_candidate(state)
                    # fall through to possibly start a new candidate in the same frame

            # -------- No active candidate: decide whether to start one --------
            owner = state["owner"]
            if owner is None:
                maybe_assign_owner(state, nearest_o, nearest_d, time_fs)
                continue

            owner_d = dist(hpos, coords[owner])
            if nearest_o == owner:
                # Same owner remains nearest. Keep following.
                continue

            # Start candidate only if new O really captured the H and old donor lost it.
            if nearest_d <= R_ENTER and owner_d >= R_LEAVE:
                start_candidate(state, owner, nearest_o, time_fs, iter_idx)
                continue

            # Otherwise, if old owner is no longer chemically close and some other O is,
            # gently re-assign ownership so later frames do not get stuck forever.
            if owner_d > R_LEAVE and nearest_d <= R_ASSIGN:
                state["owner"] = nearest_o
                state["owner_since_fs"] = time_fs
                reset_candidate(state)

    process_frame(t0, it0, co0)
    for t, itnum, _, coords in it:
        process_frame(t, itnum, coords)

    return {
        "status": "OK",
        "n_O": len(O_idx),
        "n_H": len(H_idx),
        "n_transfer_events": len(events),
        "n_H_with_transfer": len(H_with_transfer),
    }, events


def list_tdl_dirs(base):
    out = []
    for p in base.iterdir():
        if p.is_dir() and p.name.startswith("tdl"):
            m = re.match(r"tdl0*([0-9]+)$", p.name)
            key = int(m.group(1)) if m else 10**12
            out.append((key, p))
    out.sort(key=lambda x: x[0])
    return [p for _, p in out]


def bucket(k):
    if k >= 5:
        return "5+"
    return str(k)


def main():
    base = Path(".").resolve()
    tdl_dirs = list_tdl_dirs(base)

    status_counts = Counter()
    dist_counts = Counter()

    ok_runs = 0
    runs_with_transfer = 0
    total_transfer_events = 0
    total_H_with_transfer = 0
    total_monomers = 0
    sum_events_per_monomer = 0.0
    sum_H_per_monomer = 0.0

    with open(OUTCSV_RUN, "w", newline="") as f_run:
        with open(OUTCSV_EVENTS, "w", newline="") as f_evt:
            w_run = csv.writer(f_run)
            w_evt = csv.writer(f_evt)

            w_run.writerow([
                "run",
                "dir",
                "status",
                "n_O",
                "n_H",
                "monomers",
                "transfer_events_total",
                "H_with_transfer",
                "has_transfer",
                "transfer_events_per_monomer",
                "transfer_H_per_monomer"
            ])

            w_evt.writerow([
                "run",
                "dir",
                "h_index",
                "old_owner_index",
                "new_owner_index",
                "start_time_fs",
                "start_iter_idx",
                "stable_time_fs",
                "stable_iter_idx",
                "confirm_time_fs",
                "confirm_iter_idx",
                "oh_new_A",
                "oh_old_A"
            ])

            for d in tdl_dirs:
                name = d.name
                m = re.match(r"tdl0*([0-9]+)$", name)
                run = int(m.group(1)) if m else ""

                xyz = d / XYZ_NAME
                if not xyz.exists():
                    status_counts["NO_XYZ"] += 1
                    w_run.writerow([run, name, "NO_XYZ", "", "", "", "", "", "", "", ""])
                    continue

                try:
                    res, events = count_rearrangement_one_xyz(xyz)
                except Exception:
                    status_counts["ERROR"] += 1
                    w_run.writerow([run, name, "ERROR", "", "", "", "", "", "", "", ""])
                    continue

                status_counts[res["status"]] += 1
                monomers = res["n_O"]
                ev_pm = (res["n_transfer_events"] / monomers) if monomers > 0 else ""
                h_pm = (res["n_H_with_transfer"] / monomers) if monomers > 0 else ""
                has_transfer = 1 if res["n_transfer_events"] > 0 else 0

                w_run.writerow([
                    run,
                    name,
                    res["status"],
                    res["n_O"],
                    res["n_H"],
                    monomers,
                    res["n_transfer_events"],
                    res["n_H_with_transfer"],
                    has_transfer,
                    ("%.6f" % ev_pm) if ev_pm != "" else "",
                    ("%.6f" % h_pm) if h_pm != "" else "",
                ])

                for ev in events:
                    w_evt.writerow([
                        run,
                        name,
                        ev["h_index"],
                        ev["old_owner_index"],
                        ev["new_owner_index"],
                        "" if ev["start_time_fs"] is None else ("%.6f" % ev["start_time_fs"]),
                        "" if ev["start_iter_idx"] is None else ev["start_iter_idx"],
                        "" if ev["stable_time_fs"] is None else ("%.6f" % ev["stable_time_fs"]),
                        "" if ev["stable_iter_idx"] is None else ev["stable_iter_idx"],
                        "" if ev["confirm_time_fs"] is None else ("%.6f" % ev["confirm_time_fs"]),
                        "" if ev["confirm_iter_idx"] is None else ev["confirm_iter_idx"],
                        "%.6f" % ev["oh_new_A"],
                        "%.6f" % ev["oh_old_A"],
                    ])

                if res["status"] == "OK":
                    ok_runs += 1
                    total_transfer_events += res["n_transfer_events"]
                    total_H_with_transfer += res["n_H_with_transfer"]
                    total_monomers += monomers
                    sum_events_per_monomer += (res["n_transfer_events"] / monomers) if monomers > 0 else 0.0
                    sum_H_per_monomer += (res["n_H_with_transfer"] / monomers) if monomers > 0 else 0.0
                    dist_counts[bucket(res["n_transfer_events"])] += 1
                    if res["n_transfer_events"] > 0:
                        runs_with_transfer += 1

    p_any = (runs_with_transfer / ok_runs) if ok_runs else 0.0
    avg_events_per_monomer_weighted = (total_transfer_events / total_monomers) if total_monomers else 0.0
    avg_events_per_monomer_unweighted = (sum_events_per_monomer / ok_runs) if ok_runs else 0.0
    avg_H_per_monomer_weighted = (total_H_with_transfer / total_monomers) if total_monomers else 0.0
    avg_H_per_monomer_unweighted = (sum_H_per_monomer / ok_runs) if ok_runs else 0.0

    print("=== H rearrangement overall ===")
    print("tdl dirs:", len(tdl_dirs))
    for k, v in status_counts.most_common():
        print("status %-8s : %d" % (k, v))
    print("OK runs:", ok_runs)
    print("runs with >=1 transfer:", runs_with_transfer, "P=%.3f" % p_any)
    print("total transfer events:", total_transfer_events)
    print("total unique H with >=1 transfer:", total_H_with_transfer)
    print("avg transfer-events/monomer (weighted): %.6f" % avg_events_per_monomer_weighted)
    print("avg transfer-events/monomer (per-run): %.6f" % avg_events_per_monomer_unweighted)
    print("avg transfer-H/monomer (weighted): %.6f" % avg_H_per_monomer_weighted)
    print("avg transfer-H/monomer (per-run): %.6f" % avg_H_per_monomer_unweighted)

    for key in ["0", "1", "2", "3", "4", "5+"]:
        if key in dist_counts:
            print("transfer_events=%-2s : %d" % (key, dist_counts[key]))

    with open(OUTCSV_SUM, "w", newline="") as f_sum:
        w_sum = csv.writer(f_sum)
        w_sum.writerow([
            "total_tdl_dirs",
            "ok_runs",
            "runs_with_transfer",
            "p_any_transfer",
            "total_transfer_events",
            "total_unique_H_with_transfer",
            "total_monomers",
            "avg_transfer_events_per_monomer_weighted",
            "avg_transfer_events_per_monomer_unweighted",
            "avg_transfer_H_per_monomer_weighted",
            "avg_transfer_H_per_monomer_unweighted",
            "dist0",
            "dist1",
            "dist2",
            "dist3",
            "dist4",
            "dist5plus"
        ])
        w_sum.writerow([
            len(tdl_dirs),
            ok_runs,
            runs_with_transfer,
            "%.6f" % p_any,
            total_transfer_events,
            total_H_with_transfer,
            total_monomers,
            "%.6f" % avg_events_per_monomer_weighted,
            "%.6f" % avg_events_per_monomer_unweighted,
            "%.6f" % avg_H_per_monomer_weighted,
            "%.6f" % avg_H_per_monomer_unweighted,
            dist_counts.get("0", 0),
            dist_counts.get("1", 0),
            dist_counts.get("2", 0),
            dist_counts.get("3", 0),
            dist_counts.get("4", 0),
            dist_counts.get("5+", 0)
        ])

    print("Wrote:", OUTCSV_RUN)
    print("Wrote:", OUTCSV_EVENTS)
    print("Wrote:", OUTCSV_SUM)


if __name__ == "__main__":
    main()
