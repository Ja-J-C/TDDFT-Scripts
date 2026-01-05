#!/usr/bin/env python3

# -*- coding: utf-8 -*-



import math

import re

import csv

from pathlib import Path



# ----------------------------

# Fixed initial coordinates (Angstrom)

# Order: H1, H2, Cl1, Cl2

# ----------------------------

H1  = (-2.120318, 0.000000, 0.000000)

H2  = ( 0.682204, 0.589623, 0.000000)

Cl1 = (-1.777305, 1.231013, 0.000000)

Cl2 = ( 1.818752, 0.000000, 0.000000)



# atomic masses (amu) for COM axis (good enough for direction)

mH  = 1.00784

mCl = 35.45



VALENCE_ELECTRONS_TOTAL = 16.0

TARGET_FS = 25.0

TIME_TOL = 1e-6



re_elaser = {

    1: re.compile(r"^\s*e_laser1\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$"),

    2: re.compile(r"^\s*e_laser2\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$"),

    3: re.compile(r"^\s*e_laser3\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$"),

}



def vsub(a, b):

    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])



def vadd(a, b):

    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])



def vmul(s, a):

    return (s*a[0], s*a[1], s*a[2])



def dot(a, b):

    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]



def norm(a):

    return math.sqrt(dot(a, a))



def unit(a):

    n = norm(a)

    if n == 0.0:

        return None

    return (a[0]/n, a[1]/n, a[2]/n)



def angle_deg_axis(u, v):

    """

    angle between two *axes* (sign irrelevant): acos(|u·v|) in degrees, range [0, 90]

    """

    if u is None or v is None:

        return None

    c = abs(dot(u, v))

    # guard floating errors

    c = min(1.0, max(0.0, c))

    return math.degrees(math.acos(c))



# ----------------------------

# Define axes from initial geometry

# ----------------------------

# HCl bond axes (direction irrelevant for axis)

hcl1_axis = unit(vsub(H1, Cl1))  # H1 - Cl1

hcl2_axis = unit(vsub(H2, Cl2))  # H2 - Cl2



# COM axis between monomer COMs: COM(H1,Cl1) -> COM(H2,Cl2)

com1 = vmul(1.0/(mH+mCl), vadd(vmul(mH, H1),  vmul(mCl, Cl1)))

com2 = vmul(1.0/(mH+mCl), vadd(vmul(mH, H2),  vmul(mCl, Cl2)))

com_axis = unit(vsub(com2, com1))



def parse_laser_unitvec(control_path: Path):

    """

    Read e_laser1/2/3 from control.inp, return unit vector k_hat.

    """

    if not control_path.exists():

        return None

    vals = {1: None, 2: None, 3: None}

    with control_path.open("r", errors="ignore") as f:

        for line in f:

            for idx, rx in re_elaser.items():

                m = rx.match(line)

                if m:

                    vals[idx] = float(m.group(1))

    if any(vals[i] is None for i in (1,2,3)):

        return None

    E = (vals[1], vals[2], vals[3])

    return unit(E)



def parse_ionization_25fs(monitor_path: Path):

    """

    Find INFOLINE at t=25 fs in monitor.out:

    INFOLINE: <time> <Ne_remaining> ...

    ion25 = 16 - Ne_remaining

    """

    if not monitor_path.exists():

        return None

    ne = None

    with monitor_path.open("r", errors="ignore") as f:

        for line in f:

            if "INFOLINE:" not in line:

                continue

            parts = line.split()

            if len(parts) < 3:

                continue

            # parts[0] should be "INFOLINE:"

            try:

                t = float(parts[1])

            except ValueError:

                continue

            if abs(t - TARGET_FS) <= TIME_TOL:

                try:

                    ne = float(parts[2])

                except ValueError:

                    ne = None

                break

    if ne is None:

        return None

    return VALENCE_ELECTRONS_TOTAL - ne



def fmt(x):

    return "" if x is None else f"{x:.6f}"



def main():

    base = Path(".").resolve()

    out_csv = base / "laser_angles_ion25.csv"



    rows = []

    for d in sorted(base.glob("tdl*")):

        if not d.is_dir():

            continue



        run = d.name

        control = d / "control.inp"

        monitor = d / "monitor.out"



        k_hat = parse_laser_unitvec(control)

        ion25 = parse_ionization_25fs(monitor)



        theta_com  = angle_deg_axis(k_hat, com_axis)

        theta_hcl1 = angle_deg_axis(k_hat, hcl1_axis)

        theta_hcl2 = angle_deg_axis(k_hat, hcl2_axis)



        rows.append({

            "run": run,

            "ion25": ion25,

            "theta_com_deg": theta_com,

            "theta_hcl1_deg": theta_hcl1,

            "theta_hcl2_deg": theta_hcl2,

        })



    with out_csv.open("w", newline="") as f:

        w = csv.writer(f)
            "ion25": ion25,
        w.writerow(["run", "ion25", "theta_com_deg", "theta_hcl1_deg", "theta_hcl2_deg"])
        for r in rows:
            w.writerow([
                r["run"],
                fmt(r["ion25"]),
                fmt(r["theta_com_deg"]),
                fmt(r["theta_hcl1_deg"]),
                fmt(r["theta_hcl2_deg"]),
            ])
    print(f"Wrote: {out_csv}")
    print("First 10 lines:")
    with out_csv.open("r") as f:
        for i, line in enumerate(f):
            print(line.rstrip("\n"))
            if i >= 10:
                break
if __name__ == "__main__":
    main()
