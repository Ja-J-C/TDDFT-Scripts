#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, math
from pathlib import Path

# ===== 可调参数 =====
R_CUTOFF = 2.0          # C-H 距离阈值 (Å)
MIN_LEN_FRAMES = 1      # 最短连续帧数；如果想过滤掉很短的区间可以改成 3 或 5
NE_THRESHOLD = 9.5      # 判定“电离”用的电子数阈值
TRAJ_DT_FS = 0.5        # trajectory.xyz 每一帧间隔的时间 (fs)，从 0 fs 开始
# monitor.out 自己带有时间，不再乘 dt


# ================= 通用小工具 =================

def dist2(a, b):
    """三维点的平方距离"""
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2


def frames_from_xyz(fname):
    """
    逐帧读标准 XYZ：
        N
        comment
        El x y z
        ...
    yield: [(el, (x,y,z)), ...]
    """
    with open(fname) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            if not line.isdigit():
                # 忽略乱七八糟的行，直到读到原子数
                continue
            nat = int(line)
            _ = fh.readline()  # comment 行丢掉
            frame = []
            for _ in range(nat):
                el, *xyz = fh.readline().split()[:4]
                frame.append((el, tuple(map(float, xyz))))
            yield frame


# ============== 从 monitor.out 抽取 (time, Ne) 序列 ==============

def load_ne_sequence(monitor_path):
    """
    从 monitor.out 中读取时间和电子数序列 (t_fs, Ne)。
    INFOLINE 行示例：
    INFOLINE:  3.9999900000000002E+002  9.4880000371352313E+000  -5.97E+002 ...

    处理方式：
      - 找出该行中所有能 parse 成 float 的 token
      - 第一个 float 是时间 t_fs，第二个 float 是电子数 Ne

    返回:
        times_fs: [t0, t1, ...]  # 单位 fs
        nes:      [Ne0, Ne1, ...]
    """
    times = []
    nes = []
    with open(monitor_path) as f:
        for line in f:
            if "INFOLINE" not in line:
                continue
            parts = line.split()
            floats = []
            for tok in parts:
                try:
                    val = float(tok)
                    floats.append(val)
                except ValueError:
                    continue
            if len(floats) >= 2:
                t_fs = floats[0]
                ne = floats[1]
                times.append(t_fs)
                nes.append(ne)
    return times, nes


# ============== 找连续 r <= R_CUTOFF 的区间 ==============

def find_intervals(seq, cutoff, min_len_frames):
    """
    对给定的距离序列 seq，找出所有连续满足 r <= cutoff 的区间。

    返回: [(start_idx, end_idx), ...]
    """
    intervals = []
    start = None

    for i, r in enumerate(seq):
        if r <= cutoff:
            if start is None:
                start = i
        else:
            if start is not None:
                end = i - 1
                if end - start + 1 >= min_len_frames:
                    intervals.append((start, end))
                start = None

    # 如果一直到最后都在 cutoff 以内
    if start is not None:
        end = len(seq) - 1
        if end - start + 1 >= min_len_frames:
            intervals.append((start, end))

    return intervals


def ch_intervals_for_traj(traj_path, monitor_path=None,
                          r_cut=R_CUTOFF,
                          min_len_frames=MIN_LEN_FRAMES,
                          ne_threshold=NE_THRESHOLD,
                          traj_dt_fs=TRAJ_DT_FS):
    """
    对单个 trajectory.xyz，找出所有 C–H 对在全程模拟中
    距离连续 <= r_cut 的时间区间。

    几何时间：t_geom = frame_idx * traj_dt_fs，从 0 fs 开始。
    电离时间：直接使用 monitor.out 里的时间列（单位 fs）。

    返回:
    {
      "pair_intervals": {
          (c_idx, h_idx): [
              (start_idx, end_idx, t_start_fs, t_end_fs, duration_fs),
              ...
          ],
          ...
      },
      "t_ion_fs": 第一次 Ne < ne_threshold 的时间 (fs, 若没有则为 None),
      "n_frames": 总帧数
    }
    """

    # ---- 读几何帧 ----
    frames = list(frames_from_xyz(traj_path))
    n_frames = len(frames)
    if n_frames == 0:
        return {"pair_intervals": {}, "t_ion_fs": None, "n_frames": 0}

    # 几何帧对应的时间（fs）
    times_geom_fs = [i * traj_dt_fs for i in range(n_frames)]

    # ---- 读 Ne 序列，用于找第一次 < 阈值 ----
    t_ion_fs = None
    if monitor_path is not None and monitor_path.is_file():
        times_ne_fs, nes_all = load_ne_sequence(monitor_path)
        for t_fs, ne in zip(times_ne_fs, nes_all):
            if ne < ne_threshold:
                t_ion_fs = t_fs  # monitor.out 里已经是 fs
                break

    # ---- 找出 C / H 原子索引 ----
    first_frame = frames[0]
    elems = [el for el, _ in first_frame]
    c_idx = [i for i, el in enumerate(elems) if el.upper() == 'C']
    h_idx = [i for i, el in enumerate(elems) if el.upper() == 'H']

    pair_intervals = {}

    if not c_idx or not h_idx:
        return {
            "pair_intervals": {},
            "t_ion_fs": t_ion_fs,
            "n_frames": n_frames
        }

    # ---- 对每个 C–H pair 计算完整距离序列 ----
    for ci in c_idx:
        for hi in h_idx:
            seq = []
            for frame in frames:
                coords = [coord for _, coord in frame]
                r = math.sqrt(dist2(coords[ci], coords[hi]))
                seq.append(r)

            # 找 r <= r_cut 的连续区间
            raw_intervals = find_intervals(seq, r_cut, min_len_frames)
            if not raw_intervals:
                continue

            enriched = []
            for start_idx, end_idx in raw_intervals:
                t_start_fs = times_geom_fs[start_idx]
                t_end_fs = times_geom_fs[end_idx]
                duration_fs = t_end_fs - t_start_fs
                enriched.append(
                    (start_idx, end_idx, t_start_fs, t_end_fs, duration_fs)
                )

            pair_intervals[(ci, hi)] = enriched

    return {
        "pair_intervals": pair_intervals,
        "t_ion_fs": t_ion_fs,
        "n_frames": n_frames
    }


# ============== 遍历 td* 目录 ==============

def find_trajectories(root):
    trajs = []
    for dirpath, dirnames, filenames in os.walk(root):
        base = os.path.basename(dirpath)
        if base.startswith('td') and 'trajectory.xyz' in filenames:
            trajs.append(Path(dirpath) / 'trajectory.xyz')
    return sorted(trajs)


def find_monitor_file(td_dir):
    """
    在 td 目录下找 monitor.out
    """
    candidates = [
        "monitor.out",
        "Monitor.out",
    ]
    for name in candidates:
        p = td_dir / name
        if p.is_file():
            return p
    return None


# ============== 主程序 ==============

def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    if not root.is_dir():
        sys.exit(f"Error: {root} is not a directory")

    print(f"Scanning root directory: {root}\n")

    trajectories = find_trajectories(root)
    if not trajectories:
        sys.exit("❌ No trajectory.xyz found under any td*/ directory. Please check the path.")

    print("✅ Found the following trajectory.xyz files:")
    for t in trajectories:
        print("  -", t)
    print()

    for traj in trajectories:
        td_dir = traj.parent
        tdname = td_dir.name

        monitor_file = find_monitor_file(td_dir)

        res = ch_intervals_for_traj(
            traj_path=traj,
            monitor_path=monitor_file,
            r_cut=R_CUTOFF,
            min_len_frames=MIN_LEN_FRAMES,
            ne_threshold=NE_THRESHOLD,
            traj_dt_fs=TRAJ_DT_FS
        )

        pair_intervals = res["pair_intervals"]
        t_ion_fs = res["t_ion_fs"]
        n_frames = res["n_frames"]

        print(f"=== {tdname} ===")
        print(f"  total frames (geometry): {n_frames} "
              f"(t_end ≈ {(n_frames-1)*TRAJ_DT_FS:.2f} fs)")
        if monitor_file is not None:
            print(f"  monitor: {monitor_file.name}")
        else:
            print(f"  monitor: (none, using only geometry frames)")
        if t_ion_fs is not None:
            print(f"  first Ne < {NE_THRESHOLD}: t_ion ≈ {t_ion_fs:.3f} fs")
        else:
            print(f"  Ne never drops below {NE_THRESHOLD} (or no INFOLINE entries found)")

        if not pair_intervals:
            print(f"  -> No C–H pair has a continuous interval with r <= "
                  f"{R_CUTOFF:.2f} Å (or only very short intervals were filtered out).\n")
            continue

        print(f"  C–H distance continuous ≤ {R_CUTOFF:.2f} Å, time intervals (fs):")
        for (ci, hi), intervals in sorted(pair_intervals.items()):
            # 打印为 1-based index，更接近 XYZ 行号
            c_label = f"C({ci+1})"
            h_label = f"H({hi+1})"
            print(f"    Pair {c_label} - {h_label}:")
            for (start_idx, end_idx, t_start_fs, t_end_fs, duration_fs) in intervals:
                print(
                    f"      frames [{start_idx:5d} - {end_idx:5d}]  "
                    f"t ∈ [{t_start_fs:8.2f} , {t_end_fs:8.2f}] fs  "
                    f"(Δt = {duration_fs:8.2f} fs)"
                )
        print()

if __name__ == '__main__':
    main()
