#!/usr/bin/env python3

# -*- coding: utf-8 -*-



"""

HCl dimer trajectory classifier (FIXED H–Cl pairing, no rmin)



Key change vs previous version:

- Instead of using rmin(H, any Cl), we define TWO fixed bonds:

    Bond A: H1–Cl1

    Bond B: H2–Cl2

  and detect their rupture times independently.



Rupture definition (your rule):

- For a bond distance d(t), rupture time t_b is the first time t0 such that:

    d(t0) >= R_CUT

  AND for the next HOLD_FS:

    d(t) >= R_CUT for ALL frames in [t0, t0+HOLD_FS]

  AND the trajectory must cover t0+HOLD_FS (avoid end-of-traj false positives).



Simultaneous vs sequential (4-body only):

- Δt = |t_b2 - t_b1|

- Δt <= DT_SIM => simultaneous; else sequential



Time:

- Prefer "time[fs]=" in header if present

- else parse "iter = ####" and do time = iter * DT_ITER_FS



Outputs:

- per run: t_break_1, bond_1, t_break_2, bond_2, classification

- CSV: break_summary.csv

- Summary stats at end



Assumptions:

- Each frame contains at least 2 H and 2 Cl atoms.

- "H1" means the first H appearing in the file; "Cl1" means the first Cl.

  If your file ordering differs, adjust PAIRING_MODE or implement custom pairing.

"""



import re

import csv

from pathlib import Path

from math import sqrt



# =========================

# ===== USER SETTINGS =====

# =========================



# Only scan these runs (edit as needed)

RUNS = [

    91, 93, 95, 97, 101, 102, 103, 106, 109, 112, 114, 115,

    117, 118, 120, 122, 125, 126, 127, 128, 134, 135, 136, 139,

    141, 143, 144, 147, 149, 154, 156, 157, 158, 160, 163, 167,

    168, 169, 179, 180, 181, 184, 185, 187, 188, 189, 193, 194,

    195, 198, 199, 201, 203, 204, 209, 215, 217, 219, 221, 222,

    225, 226, 227, 229, 231, 232, 233, 234, 235, 236, 239, 240,

    242, 244, 247, 249, 251, 253, 257, 258, 259, 261, 263, 264,

    265, 268, 269, 276, 277, 278, 282, 283, 284, 287, 289, 290,

]



# Directory containing many tdl* folders (script directory by default)

ROOT_DIR = Path(__file__).resolve().parent



TRAJ_NAME = "trajectory.xyz"

CSV_OUT = "break_summary.csv"



# time = iter * DT_ITER_FS (fs) if time[fs] missing

DT_ITER_FS = 0.001



# Rupture / classification parameters

R_CUT   = 2    # Å

HOLD_FS = 35.0   # fs  must stay above cutoff for next HOLD_FS

DT_SIM  = 6    # fs  simultaneous threshold



# Fixed pairing mode:

# "ORDERED": H1–Cl1 and H2–Cl2 (by appearance order in the xyz)

PAIRING_MODE = "ORDERED"



# =========================

# ===== IMPLEMENTATION =====

# =========================



HEADER_TIME_RE = re.compile(r"time\s*\[fs\]\s*=\s*([0-9Ee+\-\.]+)", re.IGNORECASE)

HEADER_ITER_RE = re.compile(r"iter\s*=\s*([0-9]+)", re.IGNORECASE)



def dist(a, b):

    dx = a[0] - b[0]

    dy = a[1] - b[1]

    dz = a[2] - b[2]

    return (dx*dx + dy*dy + dz*dz) ** 0.5



def mean(xs):

    return sum(xs) / len(xs) if xs else None



def std(xs):

    if not xs:

        return None

    m = mean(xs)

    return sqrt(sum((x - m) ** 2 for x in xs) / len(xs))



def minmax(xs):

    return (min(xs), max(xs)) if xs else (None, None)



def find_run_dir(root: Path, run_num: int):

    """Find folder names like tdl91, tdl0091, tdl00091 ..."""

    candidates = [

        root / f"tdl{run_num}",

        root / f"tdl{run_num:03d}",

        root / f"tdl{run_num:04d}",

        root / f"tdl{run_num:05d}",

    ]

    for c in candidates:

        if c.is_dir():

            return c



    exact = re.compile(rf"^tdl0*{run_num}$", re.IGNORECASE)

    hits = [p for p in root.iterdir() if p.is_dir() and exact.match(p.name)]

    if hits:

        hits.sort(key=lambda p: len(p.name))

        return hits[0]

    return None



def parse_trajectory(traj_path: Path, dt_iter_fs: float):

    """

    Parse blocks:

      # iter = ... (time[fs] may exist)

      H x y z

      H x y z

      Cl x y z

      Cl x y z

    Returns: (times_fs, frames)

      times_fs: list[float]

      frames: list[list[(el, (x,y,z))]]

    """

    times = []

    frames = []



    cur_time = None

    cur_atoms = []



    def flush():

        nonlocal cur_time, cur_atoms

        if cur_time is not None and cur_atoms:

            times.append(cur_time)

            frames.append(cur_atoms)

        cur_time = None

        cur_atoms = []



    with traj_path.open("r", encoding="utf-8", errors="ignore") as f:

        for line in f:

            s = line.strip()

            if not s:

                continue



            if s.startswith("#"):

                flush()



                m_time = HEADER_TIME_RE.search(s)

                if m_time:

                    try:

                        cur_time = float(m_time.group(1))

                        continue

                    except ValueError:

                        cur_time = None



                m_iter = HEADER_ITER_RE.search(s)

                if m_iter:

                    it = int(m_iter.group(1))

                    cur_time = it * dt_iter_fs

                else:

                    cur_time = None

                continue



            parts = s.split()

            if len(parts) >= 4:

                el = parts[0]

                try:

                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

                except ValueError:

                    continue

                cur_atoms.append((el, (x, y, z)))



    flush()

    return times, frames



def extract_indices(frame):

    """Return indices of all H and all Cl in the frame (by appearance order)."""

    elems = [el.strip().upper() for el, _ in frame]

    h_idx = [i for i, e in enumerate(elems) if e == "H"]

    cl_idx = [i for i, e in enumerate(elems) if e == "CL"]

    return h_idx, cl_idx



def get_fixed_pairs(h_idx, cl_idx):

    """

    Return fixed bond pairs as list of tuples:

      [(h_index, cl_index, "H1-Cl1"), (h_index, cl_index, "H2-Cl2")]

    """

    if PAIRING_MODE.upper() != "ORDERED":

        raise ValueError(f"Unknown PAIRING_MODE: {PAIRING_MODE}")



    if len(h_idx) < 2 or len(cl_idx) < 2:

        return None



    # Ordered pairing: first H with first Cl, second H with second Cl

    return [

        (h_idx[0], cl_idx[0], "H1-Cl1"),

        (h_idx[1], cl_idx[1], "H2-Cl2"),

    ]



def find_bond_rupture_time(times_fs, frames, h_index, cl_index, r_cut, hold_fs):

    """

    For a fixed bond distance d(t)=|H-Cl|:

    find first t0 such that d(t0)>=r_cut and d(t)>=r_cut for all frames until t0+hold_fs,

    and trajectory covers t0+hold_fs.

    """

    n = len(times_fs)

    if n == 0:

        return None



    # precompute bond distance

    dlist = []

    for frame in frames:

        coords = [xyz for _, xyz in frame]

        dlist.append(dist(coords[h_index], coords[cl_index]))



    for i in range(n):

        t0 = times_fs[i]

        if dlist[i] < r_cut:

            continue



        t_end = t0 + hold_fs

        if times_fs[-1] < t_end:

            continue



        ok = True

        j = i

        while j < n and times_fs[j] <= t_end:

            if dlist[j] < r_cut:

                ok = False

                break

            j += 1



        if ok:

            return t0



    return None



def classify_one_folder(tdl_dir: Path):

    """

    Returns:

      (t1, bond1, t2, bond2, classification, dt or None)

    where (t1,bond1) is the first rupture event in time order.

    """

    traj_path = tdl_dir / TRAJ_NAME

    if not traj_path.is_file():

        return None, None, None, None, "NO_TRAJ", None



    times, frames = parse_trajectory(traj_path, DT_ITER_FS)

    if not times or not frames:

        return None, None, None, None, "EMPTY_TRAJ", None



    h_idx, cl_idx = extract_indices(frames[0])

    pairs = get_fixed_pairs(h_idx, cl_idx)

    if pairs is None:

        return None, None, None, None, "BAD_ATOMS", None



    events = []

    for (hi, ci, label) in pairs:

        tb = find_bond_rupture_time(times, frames, hi, ci, R_CUT, HOLD_FS)

        if tb is not None:

            events.append((tb, label))



    events.sort(key=lambda x: x[0])



    if len(events) == 2:

        (t1, b1), (t2, b2) = events

        dt = t2 - t1

        cls = "4-body | simultaneous" if dt <= DT_SIM else "4-body | sequential"

        return t1, b1, t2, b2, cls, dt



    if len(events) == 1:

        (t1, b1) = events[0]

        return t1, b1, None, None, "3-body", None



    return None, None, None, None, "NO_BREAK", None



def fmt_time(x):

    return "" if x is None else f"{x:.3f}"



def main():

    root = ROOT_DIR.resolve()



    # ---- Summary accumulators ----

    counts = {}  # classification -> count

    dt_list_4body = []

    t1_list_4body = []

    t2_list_4body = []

    n_4body = 0

    n_3body = 0

    n_sim = 0

    n_seq = 0



    print(f"ROOT_DIR: {root}")

    print(f"Runs to scan: {len(RUNS)}")

    print(f"Fixed pairs: {PAIRING_MODE} (default: H1-Cl1, H2-Cl2)")

    print(f"Params: R_CUT={R_CUT} Å, HOLD_FS={HOLD_FS} fs, DT_SIM={DT_SIM} fs, DT_ITER_FS={DT_ITER_FS} fs\n")



    header = ["run", "folder", "t_break_1_fs", "bond_1", "t_break_2_fs", "bond_2", "classification"]

    rows = []



    print("{:<8s}  {:<16s}  {:>12s}  {:<8s}  {:>12s}  {:<8s}  {}".format(*header))

    print("-" * 108)



    for num in RUNS:

        run_dir = find_run_dir(root, num)

        if run_dir is None:

            cls = "NO_DIR"

            row = [f"tdl{num}", "", "", "", "", "", cls]

            rows.append(row)

            counts[cls] = counts.get(cls, 0) + 1

            print("{:<8s}  {:<16s}  {:>12s}  {:<8s}  {:>12s}  {:<8s}  {}".format(*row))

            continue



        t1, b1, t2, b2, cls, dt = classify_one_folder(run_dir)



        row = [

            f"tdl{num}",

            run_dir.name,

            fmt_time(t1),

            "" if b1 is None else b1,

            fmt_time(t2),

            "" if b2 is None else b2,

            cls,

        ]

        rows.append(row)

        counts[cls] = counts.get(cls, 0) + 1



        # summary buckets

        if cls.startswith("4-body"):

            n_4body += 1

            if dt is not None:

                dt_list_4body.append(dt)

            if t1 is not None:

                t1_list_4body.append(t1)

            if t2 is not None:

                t2_list_4body.append(t2)

            if "simultaneous" in cls:

                n_sim += 1

            if "sequential" in cls:

                n_seq += 1

        elif cls == "3-body":

            n_3body += 1



        print("{:<8s}  {:<16s}  {:>12s}  {:<8s}  {:>12s}  {:<8s}  {}".format(*row))



    # Write CSV

    out_path = (root / CSV_OUT).resolve()

    with out_path.open("w", newline="", encoding="utf-8") as f:

        w = csv.writer(f)

        w.writerow(header)

        w.writerows(rows)



    # ---- Print summary statistics ----

    print("\n" + "=" * 76)

    print("Summary statistics")

    print("=" * 76)



    print(f"Total runs requested: {len(RUNS)}")

    print(f"3-body: {n_3body}")

    print(f"4-body: {n_4body}   (simultaneous: {n_sim}, sequential: {n_seq})")



    print("\nCounts by classification:")

    for k in sorted(counts.keys()):

        print(f"  {k:24s} : {counts[k]}")



    if dt_list_4body:

        m = mean(dt_list_4body)

        s = std(dt_list_4body)

        mn, mx = minmax(dt_list_4body)

        print("\nΔt = t_break_2 - t_break_1 (4-body only):")

        print(f"  N={len(dt_list_4body)}  avg={m:.3f} fs  std={s:.3f} fs  range=[{mn:.3f}, {mx:.3f}] fs")

    else:
        print("\nΔt stats: no 4-body Δt values available.")
    if t1_list_4body:
        m = mean(t1_list_4body)
        s = std(t1_list_4body)
        mn, mx = minmax(t1_list_4body)
        print("\nFirst rupture time t_break_1 (4-body only):")
        print(f"  N={len(t1_list_4body)}  avg={m:.3f} fs  std={s:.3f} fs  range=[{mn:.3f}, {mx:.3f}] fs")
    if t2_list_4body:
        m = mean(t2_list_4body)
        s = std(t2_list_4body)
        mn, mx = minmax(t2_list_4body)
        print("\nSecond rupture time t_break_2 (4-body only):")
        print(f"  N={len(t2_list_4body)}  avg={m:.3f} fs  std={s:.3f} fs  range=[{mn:.3f}, {mx:.3f}] fs")
    print(f"\n✅ CSV written to: {out_path}")
if __name__ == "__main__":
    main()
