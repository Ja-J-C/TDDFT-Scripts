#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

h2o_dimer_com_ker.py



用途：

  对 H2O dimer 的 trajectory.xyz 做末帧 KER 后处理，但只保留两个 H2O 分子的

  质心平动动能（COM translational kinetic energy）。



特点：

  1. 在脚本最上方直接指定要遍历哪些文件夹（不用命令行传 start/end）。

  2. 只用最后两帧做后向差分，得到两个 H2O 质心速度。

  3. KER = KE(COM of H2O #1) + KE(COM of H2O #2)

  4. 不计算任何角度，也不计算分子内部转动/振动能。

  5. 仅当 trajectory.xyz 的最后一帧时间约等于 TARGET_FS 时才纳入统计。



说明：

  - 默认假设 6 个原子的顺序是 O H H O H H。

  - 如果你的原子顺序不同，只需要修改 WATER1_ATOMS / WATER2_ATOMS。

  - 如果 xyz 坐标单位是 Bohr，请把 LENGTH_UNIT 改成 "bohr"。

"""



import math

import os

import re

from pathlib import Path

from collections import deque



# ======================================================================

# 用户配置区：直接在这里改

# ======================================================================

BASE_DIR = "."

TARGET_DIRS = [

"tdl3", "tdl7", "tdl13"
]

XYZ_NAME = "trajectory.xyz"

OUTPUT_TSV = "h2o_dimer_com_ker_results.tsv"



# 只纳入“跑满”的轨迹

TARGET_FS = 300.0

TIME_TOL_FS = 0.55

FRAME_DT_FS = 0.5   # 若注释行里没有 time[fs]，则用帧数和这个步长估算



# 坐标单位："A" 或 "bohr"

LENGTH_UNIT = "A"



# 是否去掉整个二聚体体系的总质心漂移

# 默认 False：直接给实验室系下两个 H2O 质心的平动动能

REMOVE_SYSTEM_COM = False



# H2O 分组（0-based index）

# 默认对应 O H H O H H

WATER1_ATOMS = [0, 1, 2]

WATER2_ATOMS = [3, 4, 5]



# 可选：检查末帧物种顺序是否符合预期；不想检查就设成 None

EXPECTED_LABELS = ["O", "H", "H", "O", "H", "H"]



# 原子质量（amu）

MASS_AMU = {

    "H": 1.00784,

    "O": 15.999,

}

# ======================================================================



E_CHARGE = 1.602176634e-19   # J per eV

AMU = 1.66053906660e-27      # kg

BOHR = 5.29177210903e-11     # m

ANG = 1.0e-10                # m



# 解析“无 E 科学计数”，如 1.23456+001 / -5.0-002

_SCIFIX = re.compile(

    r"""

    ^\s*

    (?P<sgn>[-+]?)

    (?:

      (?P<int>\d+)?

      (?:\.(?P<frac>\d+))?

     |

      \.(?P<onlyfrac>\d+)

    )

    (?P<exp>[+-]\d{3})

    \s*$

    """,

    re.VERBOSE,

)



_TIME_PAT = re.compile(r"time\s*\[\s*fs\s*\]\s*=\s*([^\s]+)", re.IGNORECASE)





def parse_num(tok):

    s = tok.strip()

    if not s:

        raise ValueError("empty token")



    s_std = s.replace("D", "E").replace("d", "E")

    try:

        return float(s_std)

    except ValueError:

        pass



    m = _SCIFIX.match(s)

    if m:

        sgn = -1.0 if m.group("sgn") == "-" else 1.0

        intp = m.group("int") or "0"

        frac = m.group("frac") or (m.group("onlyfrac") or "")

        mant = float(f"{intp}.{frac}") if frac else float(intp)

        expo = int(m.group("exp"))

        return sgn * mant * (10.0 ** expo)



    m2 = re.search(r"([+-]\d{3})\s*$", s)

    if m2:

        return float(s[:m2.start()] + "E" + s[m2.start() + 1 :])



    raise ValueError(f"cannot parse number: {tok!r}")





def parse_time_fs_from_comment(line):

    m = _TIME_PAT.search(line)

    if not m:

        return None

    try:

        return float(parse_num(m.group(1)))

    except Exception:

        return None





def read_last_k_frames_with_times(xyz_path, length_scale_m, keep_k):

    """

    返回 (frames_list, times_list, total_frames_count)

      frames_list: [ [(sym,x,y,z), ...], ... ]

      times_list : [ t_fs or None, ... ]

    """

    dq_frames = deque(maxlen=max(1, keep_k))

    dq_times = deque(maxlen=max(1, keep_k))

    total = 0



    with xyz_path.open("r", encoding="utf-8", errors="ignore") as f:

        while True:

            line = f.readline()

            if not line:

                break

            line = line.strip()

            if not line:

                continue



            try:

                nat = int(line)

            except ValueError:

                continue



            comment = f.readline()

            if not comment:

                break

            t_fs = parse_time_fs_from_comment(comment)



            atoms = []

            ok = True

            for _ in range(nat):

                atom_line = f.readline()

                if not atom_line:

                    ok = False

                    break

                parts = atom_line.split()

                if len(parts) < 4:

                    ok = False

                    break

                sym = parts[0]

                try:

                    x = parse_num(parts[1]) * length_scale_m

                    y = parse_num(parts[2]) * length_scale_m

                    z = parse_num(parts[3]) * length_scale_m

                except Exception:

                    ok = False

                    break

                atoms.append((sym, x, y, z))



            if ok and len(atoms) == nat:

                dq_frames.append(atoms)

                dq_times.append(t_fs)

                total += 1



    return list(dq_frames), list(dq_times), total





def vec_sub(a, b):

    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])





def vec_sqnorm(a):

    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2]





def molecule_mass_kg(frame, atom_indices):

    total_amu = 0.0

    for idx in atom_indices:

        sym = frame[idx][0].strip().capitalize()

        if sym not in MASS_AMU:

            raise ValueError(f"Unknown atomic symbol for mass lookup: {sym!r}")

        total_amu += MASS_AMU[sym]

    return total_amu * AMU





def molecule_com_position(frame, atom_indices):

    total_mass = 0.0

    sx = sy = sz = 0.0

    for idx in atom_indices:

        sym, x, y, z = frame[idx]

        key = sym.strip().capitalize()

        if key not in MASS_AMU:

            raise ValueError(f"Unknown atomic symbol for mass lookup: {sym!r}")

        m = MASS_AMU[key]

        total_mass += m

        sx += m * x

        sy += m * y

        sz += m * z

    if total_mass <= 0.0:

        raise ValueError("Non-positive molecular mass")

    return (sx / total_mass, sy / total_mass, sz / total_mass)





def com_velocity(frame_prev, frame_last, atom_indices, dt_s):

    r0 = molecule_com_position(frame_prev, atom_indices)

    r1 = molecule_com_position(frame_last, atom_indices)

    return ((r1[0] - r0[0]) / dt_s, (r1[1] - r0[1]) / dt_s, (r1[2] - r0[2]) / dt_s)





def ke_from_mass_and_v(mass_kg, vxyz):

    return 0.5 * mass_kg * vec_sqnorm(vxyz) / E_CHARGE





def stats(vals):

    n = len(vals)

    if n == 0:

        return (None, None, None, None, 0)

    mu = sum(vals) / n

    sd = math.sqrt(sum((x - mu) * (x - mu) for x in vals) / (n - 1)) if n >= 2 else 0.0

    return (mu, sd, min(vals), max(vals), n)





def fmt(x, nd=4):

    return f"{x:.{nd}f}" if x is not None else "NA"





def main():

    if not TARGET_DIRS:

        raise ValueError("TARGET_DIRS is empty. Please set folders at the top of the script.")



    u = LENGTH_UNIT.lower()

    if u in ("a", "angstrom", "å"):

        length_scale_m = ANG

    elif u == "bohr":

        length_scale_m = BOHR

    else:

        raise ValueError(f"Unsupported LENGTH_UNIT: {LENGTH_UNIT!r}")



    dt_s = max(FRAME_DT_FS, 1e-12) * 1e-15

    base_dir = Path(BASE_DIR).resolve()



    header = "dir\tstatus\tt_last_fs\tKER_eV\tKE_H2O1_COM_eV\tKE_H2O2_COM_eV"

    out_lines = [header]



    ker_vals = []

    ke1_vals = []

    ke2_vals = []



    n_ok = 0

    n_skip = 0



    for dirname in TARGET_DIRS:

        folder = base_dir / dirname

        xyz = folder / XYZ_NAME



        if not folder.is_dir():

            print(f"{dirname} (missing dir)")

            out_lines.append(f"{dirname}\tmissing_dir\tNA\tNA\tNA\tNA")

            continue



        if not xyz.is_file():

            print(f"{dirname} (missing {XYZ_NAME})")

            out_lines.append(f"{dirname}\tmissing_xyz\tNA\tNA\tNA\tNA")

            continue



        frames, times, total_frames = read_last_k_frames_with_times(xyz, length_scale_m, keep_k=2)

        if not frames:

            print(f"{dirname} (empty)")

            out_lines.append(f"{dirname}\tempty\tNA\tNA\tNA\tNA")

            continue



        nat = len(frames[-1])

        if nat < 6:

            print(f"{dirname} (atom count < 6, got {nat})")

            out_lines.append(f"{dirname}\ttoo_few_atoms\tNA\tNA\tNA\tNA")

            continue



        if max(WATER1_ATOMS + WATER2_ATOMS) >= nat:

            print(f"{dirname} (WATER1_ATOMS/WATER2_ATOMS exceeds atom count {nat})")

            out_lines.append(f"{dirname}\tbad_atom_indices\tNA\tNA\tNA\tNA")

            continue



        if EXPECTED_LABELS is not None:

            last_species = [atom[0].strip().capitalize() for atom in frames[-1]]

            exp_species = [s.strip().capitalize() for s in EXPECTED_LABELS]

            if len(exp_species) > nat or last_species[: len(exp_species)] != exp_species:

                print(f"{dirname} (species order mismatch with EXPECTED_LABELS) — skipped")

                out_lines.append(f"{dirname}\tbad_species\tNA\tNA\tNA\tNA")

                continue



        t_last_comment = times[-1] if times else None

        if t_last_comment is not None:

            t_last_fs = float(t_last_comment)

        else:

            t_last_fs = (max(total_frames, 1) - 1) * FRAME_DT_FS



        if abs(t_last_fs - TARGET_FS) > TIME_TOL_FS:

            print(f"{dirname}  INCOMPLETE (last t = {t_last_fs:.3f} fs, need ~{TARGET_FS} fs) — skipped")

            out_lines.append(f"{dirname}\tincomplete\t{t_last_fs:.6f}\tNA\tNA\tNA")

            n_skip += 1

            continue



        if len(frames) < 2:

            print(f"{dirname} (only one frame, cannot compute velocity)")

            out_lines.append(f"{dirname}\tone_frame\t{t_last_fs:.6f}\tNA\tNA\tNA")

            continue



        frame_prev, frame_last = frames[-2], frames[-1]



        try:

            v1 = com_velocity(frame_prev, frame_last, WATER1_ATOMS, dt_s)

            v2 = com_velocity(frame_prev, frame_last, WATER2_ATOMS, dt_s)



            m1 = molecule_mass_kg(frame_last, WATER1_ATOMS)

            m2 = molecule_mass_kg(frame_last, WATER2_ATOMS)

        except Exception as exc:

            print(f"{dirname} (failed to compute COM velocity: {exc})")

            out_lines.append(f"{dirname}\tcompute_error\t{t_last_fs:.6f}\tNA\tNA\tNA")

            continue



        if REMOVE_SYSTEM_COM:

            m_tot = m1 + m2

            v_sys = (

                (m1 * v1[0] + m2 * v2[0]) / m_tot,

                (m1 * v1[1] + m2 * v2[1]) / m_tot,

                (m1 * v1[2] + m2 * v2[2]) / m_tot,

            )

            v1 = vec_sub(v1, v_sys)

            v2 = vec_sub(v2, v_sys)



        ke1 = ke_from_mass_and_v(m1, v1)

        ke2 = ke_from_mass_and_v(m2, v2)

        ker = ke1 + ke2



        ke1_vals.append(ke1)

        ke2_vals.append(ke2)

        ker_vals.append(ker)

        n_ok += 1



        print(

            f"{dirname:<12s} t_last={t_last_fs:7.3f} fs  |  "

            f"KER={ker:10.4f} eV  |  H2O1_COM={ke1:9.4f}  H2O2_COM={ke2:9.4f}"

        )

        out_lines.append(

            f"{dirname}\tOK\t{t_last_fs:.6f}\t{ker:.10f}\t{ke1:.10f}\t{ke2:.10f}"

        )



    print("----- Summary statistics over valid (full-length) entries -----")



    def pr_stats(name, vals, nd=4):

        mu, sd, vmin, vmax, n = stats(vals)

        rng = f"[{fmt(vmin, nd)} , {fmt(vmax, nd)}]" if n > 0 else "NA"

        print(f"{name:20s} avg={fmt(mu, nd)}   std={fmt(sd, nd)}   range={rng}   N={n}")



    pr_stats("KER (eV)", ker_vals, 4)

    pr_stats("KE_H2O1_COM (eV)", ke1_vals, 4)

    pr_stats("KE_H2O2_COM (eV)", ke2_vals, 4)

    print(f"Counts: OK={n_ok}, skipped(incomplete)={n_skip}")

    def stat_tuple(vals):

        return stats(vals)[:4]

    mean_row = [

        "MEAN", "-", "-",

        fmt(stat_tuple(ker_vals)[0], 4),

        fmt(stat_tuple(ke1_vals)[0], 4),

        fmt(stat_tuple(ke2_vals)[0], 4),

    ]

    std_row = [

        "STD", "-", "-",

        fmt(stat_tuple(ker_vals)[1], 4),

        fmt(stat_tuple(ke1_vals)[1], 4),

        fmt(stat_tuple(ke2_vals)[1], 4),
    ]
    min_row = [
        "MIN", "-", "-",
        fmt(stat_tuple(ker_vals)[2], 4),
        fmt(stat_tuple(ke1_vals)[2], 4),
        fmt(stat_tuple(ke2_vals)[2], 4),
    ]
    max_row = [
        "MAX", "-", "-",
        fmt(stat_tuple(ker_vals)[3], 4),
        fmt(stat_tuple(ke1_vals)[3], 4),
        fmt(stat_tuple(ke2_vals)[3], 4),
    ]
    out_lines += ["\t".join(mean_row), "\t".join(std_row), "\t".join(min_row), "\t".join(max_row)]
    out_path = base_dir / OUTPUT_TSV
    with out_path.open("w", encoding="utf-8") as fo:
        fo.write("\n".join(out_lines) + "\n")
    print(f"Results saved to: {out_path}")
if __name__ == "__main__":
    main()
