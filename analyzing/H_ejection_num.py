#!/usr/bin/env python3
import re
import math
import csv
from pathlib import Path
from collections import Counter

CAP_BOUND = 11.0
D_CLEAN = 3.0
D_CAP = 2.5
FRAG_CUTOFF = 2.25
DT_FS_PER_ITER = 0.001

XYZ_NAME = "trajectory.xyz"
OUTCSV_RUN = "h_ejection_by_run.csv"
OUTCSV_EVENTS = "h_ejection_events.csv"
OUTCSV_SUM = "h_ejection_overall.csv"

TIME_RE = re.compile(r"time\[fs\]\s*=\s*([0-9Ee+\-\.]+)")
ITER_RE = re.compile(r"iter\s*=\s*([0-9]+)")


def in_cap(x, y, z, bound=CAP_BOUND):
    return (abs(x) >= bound) or (abs(y) >= bound) or (abs(z) >= bound)


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

                sym = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])

                species.append(sym)
                coords.append((x, y, z))

            yield time_fs, iter_idx, species, coords


def connected_components(coords, cutoff=FRAG_CUTOFF):
    """
    Build fragments by connecting atoms whose nuclear distance < cutoff.
    Returns list of components, each component is a list of atom indices.
    """
    n = len(coords)
    visited = [False] * n
    comps = []

    # adjacency on the fly
    for i in range(n):
        if visited[i]:
            continue

        stack = [i]
        visited[i] = True
        comp = []

        while stack:
            u = stack.pop()
            comp.append(u)

            for v in range(n):
                if visited[v]:
                    continue
                if dist(coords[u], coords[v]) < cutoff:
                    visited[v] = True
                    stack.append(v)

        comps.append(comp)

    return comps


def count_ejection_one_xyz(xyz_path):
    """
    For each H atom:
      1. Scan from the beginning of trajectory.
      2. In each frame, find its nearest O.
      3. If nearest-OH distance > 3 A and both H and nearest O are not in CAP:
         count as clean ejection and stop tracking this H.
      4. If H or nearest O reaches CAP first:
         at that first CAP-touch frame, if nearest-OH distance > 2.5 A,
         count as CAP-type ejection; otherwise do not count it.
         Stop tracking this H.

    Additional CAP contamination filter:
      - In every frame, define fragments by nuclear-distance connectivity
        (< 2.25 A).
      - If any atom in a fragment touches CAP, then all H atoms in that fragment
        are permanently marked as contaminated.
      - Contaminated H atoms are excluded from all future H-ejection counting.
    """
    it = iter_xyz_frames(xyz_path)

    try:
        t0, it0, sp0, co0 = next(it)
    except StopIteration:
        return {
            "status": "EMPTY",
            "n_O": 0,
            "n_H": 0,
            "e_tot": 0,
            "e_clean": 0,
            "e_cap": 0,
            "n_contaminated_H": 0
        }, []

    O_idx = [i for i, s in enumerate(sp0) if s.upper() == "O"]
    H_idx = [i for i, s in enumerate(sp0) if s.upper() == "H"]

    if not O_idx or not H_idx:
        return {
            "status": "NO_O_OR_H",
            "n_O": len(O_idx),
            "n_H": len(H_idx),
            "e_tot": 0,
            "e_clean": 0,
            "e_cap": 0,
            "n_contaminated_H": 0
        }, []

    active = {h: True for h in H_idx}
    contaminated = {h: False for h in H_idx}
    e_clean = 0
    e_cap = 0
    events = []

    def mark_cap_contaminated(coords):
        comps = connected_components(coords, cutoff=FRAG_CUTOFF)

        for comp in comps:
            comp_touches_cap = False
            for idx in comp:
                x, y, z = coords[idx]
                if in_cap(x, y, z):
                    comp_touches_cap = True
                    break

            if comp_touches_cap:
                for idx in comp:
                    if idx in contaminated:
                        contaminated[idx] = True
                        # once contaminated, never count this H in future
                        active[idx] = False

    def process_frame(time_fs, iter_idx, coords):
        nonlocal e_clean, e_cap

        # first apply CAP contamination memory by fragment
        mark_cap_contaminated(coords)

        O_coords = [(i, coords[i]) for i in O_idx]

        for h in H_idx:
            if not active[h]:
                continue
            if contaminated[h]:
                active[h] = False
                continue

            hpos = coords[h]

            best_d = 1.0e99
            best_o = None
            best_opos = None

            for oi, opos in O_coords:
                d = dist(hpos, opos)
                if d < best_d:
                    best_d = d
                    best_o = oi
                    best_opos = opos

            h_in = in_cap(hpos[0], hpos[1], hpos[2])
            o_in = in_cap(best_opos[0], best_opos[1], best_opos[2])

            if h_in or o_in:
                if best_d > D_CAP:
                    e_cap += 1
                    events.append({
                        "h_index": h + 1,
                        "time_fs": time_fs,
                        "iter_idx": iter_idx,
                        "event_type": "cap",
                        "nearest_o_index": best_o + 1,
                        "oh_distance_A": best_d
                    })
                active[h] = False
                continue

            if best_d > D_CLEAN:
                e_clean += 1
                events.append({
                    "h_index": h + 1,
                    "time_fs": time_fs,
                    "iter_idx": iter_idx,
                    "event_type": "clean",
                    "nearest_o_index": best_o + 1,
                    "oh_distance_A": best_d
                })
                active[h] = False
                continue

    process_frame(t0, it0, co0)

    for t, itnum, _, coords in it:
        if all(not v for v in active.values()):
            break
        process_frame(t, itnum, coords)

    e_tot = e_clean + e_cap
    n_contaminated_H = sum(1 for h in H_idx if contaminated[h])

    return {
        "status": "OK",
        "n_O": len(O_idx),
        "n_H": len(H_idx),
        "e_tot": e_tot,
        "e_clean": e_clean,
        "e_cap": e_cap,
        "n_contaminated_H": n_contaminated_H
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
    runs_with_eject = 0
    total_e = 0
    total_clean = 0
    total_cap = 0
    total_monomers = 0
    total_contaminated_H = 0
    sum_epm = 0.0

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
                "contaminated_H",
                "eject_total",
                "eject_clean",
                "eject_cap",
                "eject_per_monomer"
            ])

            w_evt.writerow([
                "run",
                "dir",
                "h_index",
                "time_fs",
                "iter_idx",
                "event_type",
                "nearest_o_index",
                "oh_distance_A"
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
                    res, events = count_ejection_one_xyz(xyz)
                except Exception:
                    status_counts["ERROR"] += 1
                    w_run.writerow([run, name, "ERROR", "", "", "", "", "", "", "", ""])
                    continue

                status_counts[res["status"]] += 1

                monomers = res["n_O"]
                epm = (res["e_tot"] / monomers) if monomers > 0 else ""

                w_run.writerow([
                    run,
                    name,
                    res["status"],
                    res["n_O"],
                    res["n_H"],
                    monomers,
                    res["n_contaminated_H"],
                    res["e_tot"],
                    res["e_clean"],
                    res["e_cap"],
                    ("%.6f" % epm) if epm != "" else ""
                ])

                for ev in events:
                    w_evt.writerow([
                        run,
                        name,
                        ev["h_index"],
                        "" if ev["time_fs"] is None else ("%.6f" % ev["time_fs"]),
                        "" if ev["iter_idx"] is None else ev["iter_idx"],
                        ev["event_type"],
                        ev["nearest_o_index"],
                        "%.6f" % ev["oh_distance_A"]
                    ])

                if res["status"] == "OK":
                    ok_runs += 1
                    total_e += res["e_tot"]
                    total_clean += res["e_clean"]
                    total_cap += res["e_cap"]
                    total_monomers += monomers
                    total_contaminated_H += res["n_contaminated_H"]
                    sum_epm += (res["e_tot"] / monomers) if monomers > 0 else 0.0
                    dist_counts[bucket(res["e_tot"])] += 1
                    if res["e_tot"] > 0:
                        runs_with_eject += 1

    p_any = (runs_with_eject / ok_runs) if ok_runs else 0.0
    avg_epm_weighted = (total_e / total_monomers) if total_monomers else 0.0
    avg_epm_unweighted = (sum_epm / ok_runs) if ok_runs else 0.0
    avg_contam_H_per_run = (total_contaminated_H / ok_runs) if ok_runs else 0.0

    print("=== H ejection overall ===")
    print("tdl dirs:", len(tdl_dirs))
    for k, v in status_counts.most_common():
        print("status %-8s : %d" % (k, v))
    print("OK runs:", ok_runs)
    print("runs with >=1 eject:", runs_with_eject, "P=%.3f" % p_any)
    print("total ejected H:", total_e, "(clean=%d, cap=%d)" % (total_clean, total_cap))
    print("total contaminated H:", total_contaminated_H)
    print("avg contaminated H per run: %.6f" % avg_contam_H_per_run)
    print("avg eject/monomer (weighted): %.6f" % avg_epm_weighted)
    print("avg eject/monomer (per-run): %.6f" % avg_epm_unweighted)

    for key in ["0", "1", "2", "3", "4", "5+"]:
        if key in dist_counts:
            print("eject=%-2s : %d" % (key, dist_counts[key]))

    with open(OUTCSV_SUM, "w", newline="") as f_sum:
        w_sum = csv.writer(f_sum)
        w_sum.writerow([
            "total_tdl_dirs",
            "ok_runs",
            "runs_with_eject",
            "p_any_eject",
            "total_ejected_H",
            "total_clean",
            "total_cap",
            "total_contaminated_H",
            "total_monomers",
            "avg_contaminated_H_per_run",
            "avg_epm_weighted",
            "avg_epm_unweighted",
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
            runs_with_eject,
            "%.6f" % p_any,
            total_e,
            total_clean,
            total_cap,
            total_contaminated_H,
            total_monomers,
            "%.6f" % avg_contam_H_per_run,
            "%.6f" % avg_epm_weighted,
            "%.6f" % avg_epm_unweighted,
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
