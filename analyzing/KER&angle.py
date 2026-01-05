#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

get.py — 为 HCl dimer 库仑爆炸后处理：

  • 末帧动能 (eV)：H1/H2/Cl1/Cl2 + 总 KER

  • 发射角：θ_HH、θ_ClCl（用最后 N 帧平均的动量方向求夹角）

  • 仅当 trajectory.xyz 的“最后一帧时间 ≈ target_fs”时，才纳入统计（否则标记未跑满、并跳过）

  • 末尾输出每列的 average / range(min–max) / standard deviation（样本标准差）



用法（在 HCl 目录里）:

  python3 get.py --start 1 --end 50

常用选项：

  --remove-com               # 去质心（默认不去质心=实验室系）

  --angle-window 5           # 角度用最后 N 帧估方向（默认 3）

  --frame-dt-fs 0.5          # 帧间隔 Δt（fs），默认 0.5

  --target-fs 120            # 目标总时长（fs），默认 120

  --time-tol-fs 0.25         # “跑满”容差（fs），默认 ±0.25

  --length-unit bohr         # 坐标若是 Bohr

  --prefix td                # 子目录前缀（默认 td -> td1, td2, ...）

  --out ker_results.tsv      # 输出文件

  --labels H H Cl Cl         # 若 XYZ 顺序/大小写异常，可显式指定

"""



import argparse, os, re, sys, math



E_CHARGE = 1.602176634e-19   # J per eV

AMU      = 1.66053906660e-27 # kg

BOHR     = 5.29177210903e-11 # m

ANG      = 1.0e-10           # m



# 解析“无E科学计数”：如 1.23456+001 / -5.0-002

_scifix = re.compile(r"""

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

""", re.VERBOSE)



def parse_num(tok: str) -> float:

    s = tok.strip()

    if not s:

        raise ValueError("empty token")

    # 先试普通浮点/Fortran D 指数

    s_std = s.replace('D','E').replace('d','E')

    try:

        return float(s_std)

    except ValueError:

        pass

    # 无E指数

    m = _scifix.match(s)

    if m:

        sgn  = -1.0 if m.group('sgn')=='-' else 1.0

        intp = m.group('int') or '0'

        frac = m.group('frac') or (m.group('onlyfrac') or '')

        mant = float(f"{intp}.{frac}") if frac!='' else float(intp)

        expo = int(m.group('exp'))

        return sgn * mant * (10.0**expo)

    # 兜底：把末尾 ±NNN 改成 E±NNN

    m2 = re.search(r'([+-]\d{3})\s*$', s)

    if m2:

        return float(s[:m2.start()] + 'E' + s[m2.start()+1:])

    raise ValueError(f"cannot parse number: {tok!r}")



_time_pat = re.compile(r'time\s*\[\s*fs\s*\]\s*=\s*([^\s]+)', re.IGNORECASE)



def parse_time_fs_from_comment(line: str):

    m = _time_pat.search(line)

    if not m: return None

    try:

        return float(parse_num(m.group(1)))

    except Exception:

        return None



def read_last_k_frames_with_times(xyz_path, length_scale_m, keep_k):

    """

    读取 trajectory.xyz 的最后 keep_k 帧：

      返回 (frames_list, times_list, total_frames_count)

      * frames_list: [ [(sym,x,y,z),...], ... ]

      * times_list:  [ t_fs 或 None, ... ]（与 frames_list 同长度）

      * total_frames_count: 文件中总帧数

    """

    from collections import deque

    dq_frames = deque(maxlen=max(1, keep_k))

    dq_times  = deque(maxlen=max(1, keep_k))

    total = 0

    with open(xyz_path, 'r', encoding='utf-8', errors='ignore') as f:

        while True:

            line = f.readline()

            if not line: break

            line = line.strip()

            if not line: continue

            try:

                n = int(line)

            except ValueError:

                continue

            comment = f.readline()

            if not comment: break

            t_fs = parse_time_fs_from_comment(comment)



            atoms = []

            ok = True

            for _i in range(n):

                ln = f.readline()

                if not ln: ok=False; break

                parts = ln.split()

                if len(parts) < 4: ok=False; break

                sym = parts[0]

                try:

                    x = parse_num(parts[1])*length_scale_m

                    y = parse_num(parts[2])*length_scale_m

                    z = parse_num(parts[3])*length_scale_m

                except Exception:

                    ok=False; break

                atoms.append((sym,x,y,z))

            if ok and len(atoms)==n:

                dq_frames.append(atoms)

                dq_times.append(t_fs)

                total += 1

    return list(dq_frames), list(dq_times), total



def ke_from_v(amu_mass, vxyz):

    m = amu_mass * AMU

    vx,vy,vz = vxyz

    return 0.5 * m * (vx*vx + vy*vy + vz*vz) / E_CHARGE  # eV



def v_diff(a0, a1, i, dt_s):

    return ((a1[i][1]-a0[i][1])/dt_s, (a1[i][2]-a0[i][2])/dt_s, (a1[i][3]-a0[i][3])/dt_s)



def v_sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def v_add(a,b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def v_scale(a,s): return (a[0]*s, a[1]*s, a[2]*s)

def v_norm(a): return math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])



def angle_between(a,b):

    na, nb = v_norm(a), v_norm(b)

    if na < 1e-20 or nb < 1e-20: return None

    c = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (na*nb)

    c = 1.0 if c>1.0 else (-1.0 if c<-1.0 else c)

    return math.degrees(math.acos(c))



def stats(vals):

    """返回 (mean, std, vmin, vmax, n)。std=样本标准差（n<2 时为0）。"""

    n = len(vals)

    if n==0: return (None, None, None, None, 0)

    mu = sum(vals)/n

    sd = math.sqrt(sum((x-mu)*(x-mu) for x in vals)/(n-1)) if n>=2 else 0.0

    return (mu, sd, min(vals), max(vals), n)



def fmt(x, nd=4): return f"{x:.{nd}f}" if x is not None else "NA"



def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("--start", type=int, default=91)

    ap.add_argument("--end", type=int, default=290)

    ap.add_argument("--base", type=str, default=".")

    ap.add_argument("--prefix", type=str, default="tdl")

    ap.add_argument("--out", type=str, default="ker_results.tsv")

    ap.add_argument("--frame-dt-fs", type=float, default=0.5)

    ap.add_argument("--target-fs", type=float, default=120, help="要求的总时长 (fs)")

    ap.add_argument("--time-tol-fs", type=float, default=0.25, help="时间容差 (fs)")

    ap.add_argument("--length-unit", type=str, choices=["A","bohr","angstrom","Å"], default="A")

    ap.add_argument("--mass-H", type=float, default=1.00784)

    ap.add_argument("--mass-Cl", type=float, default=35.45)

    ap.add_argument("--remove-com", action="store_true", help="去质心（默认不去质心）")

    ap.add_argument("--labels", nargs="+", default=None, help="显式指定原子顺序，例如: --labels H H Cl Cl")

    ap.add_argument("--angle-window", type=int, default=3, help="角度估计使用最后 N 帧（默认 3）")

    args = ap.parse_args()



    # 长度单位

    u = args.length_unit.lower()

    if u in ("a","angstrom","å"): L = ANG

    elif u == "bohr":             L = BOHR

    else:

        print(f"Unsupported length unit: {args.length_unit}", file=sys.stderr); sys.exit(2)



    dt_s = max(args.frame_dt_fs, 1e-12) * 1e-15  # s



    header = "idx\tstatus\tt_last_fs\tKER_eV\tKE_H1_eV\tKE_H2_eV\tKE_Cl1_eV\tKE_Cl2_eV\ttheta_HH_deg\ttheta_ClCl_deg"

    out_lines = [header]



    L_KER, L_H1, L_H2, L_Cl1, L_Cl2, L_tHH, L_tClCl = [], [], [], [], [], [], []

    nOK, nSkip = 0, 0



    for idx in range(args.start, args.end+1):

        folder = os.path.join(args.base, f"{args.prefix}{idx}")

        xyz = os.path.join(folder, "trajectory.xyz")

        if not os.path.isdir(folder):

            print(f"{args.prefix}{idx:04d} (missing dir)")

            out_lines.append(f"{idx}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); continue

        if not os.path.isfile(xyz):

            print(f"{args.prefix}{idx:04d} (missing trajectory.xyz)")

            out_lines.append(f"{idx}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); continue



        frames, times, total_frames = read_last_k_frames_with_times(xyz, L, max(2, args.angle_window))

        if not frames:

            print(f"{args.prefix}{idx:04d} (empty)")

            out_lines.append(f"{idx}\tempty\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); continue



        nat = len(frames[-1])

        species = [a[0] for a in frames[-1]]

        if args.labels is not None:

            if len(args.labels)!=nat:

                print(f"{args.prefix}{idx:04d} (--labels length {len(args.labels)} != atom count {nat})")

                out_lines.append(f"{idx}\tbad_labels\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); continue

            species = args.labels



        h_idx  = [i for i,s in enumerate(species) if s.lower().startswith('h')][:2]

        cl_idx = [i for i,s in enumerate(species) if s.lower().startswith('cl')][:2]

        if len(h_idx)<2 or len(cl_idx)<2:

            print(f"{args.prefix}{idx:04d} (need 2 H and 2 Cl)")

            out_lines.append(f"{idx}\twrong_species\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); continue



        # ---- 检查是否“跑满 target_fs” ----

        # 末帧时间优先用注释里的 t_fs，否则用 (总帧数-1)*Δt 估算

        t_last_comment = times[-1] if len(times)>0 else None

        if t_last_comment is not None:

            t_last_fs = float(t_last_comment)

        else:

            t_last_fs = (max(total_frames,1)-1) * args.frame_dt_fs



        full = (abs(t_last_fs - args.target_fs) <= args.time_tol_fs)

        if not full:

            print(f"{args.prefix}{idx:04d}  INCOMPLETE (last t = {t_last_fs:.3f} fs, need ~{args.target_fs} fs) — skipped")

            out_lines.append(f"{idx}\tincomplete\t{t_last_fs:.6f}\tNA\tNA\tNA\tNA\tNA\tNA\tNA")

            nSkip += 1

            continue  # 不纳入统计



        # ---------- 末帧动能：用末两帧后向差分 ----------

        if len(frames) == 1:

            keH1=keH2=keCl1=keCl2=0.0; ker=0.0

            vH1_last=vH2_last=vCl1_last=vCl2_last=(0.0,0.0,0.0)

        else:

            a0, a1 = frames[-2], frames[-1]

            vH1_last  = v_diff(a0, a1, h_idx[0], dt_s)

            vH2_last  = v_diff(a0, a1, h_idx[1], dt_s)

            vCl1_last = v_diff(a0, a1, cl_idx[0], dt_s)

            vCl2_last = v_diff(a0, a1, cl_idx[1], dt_s)

            if args.remove_com:

                mH  = args.mass_H * AMU; mCl = args.mass_Cl * AMU; M = 2*mH+2*mCl

                vcm = ((mH*vH1_last[0]+mH*vH2_last[0]+mCl*vCl1_last[0]+mCl*vCl2_last[0])/M,

                       (mH*vH1_last[1]+mH*vH2_last[1]+mCl*vCl1_last[1]+mCl*vCl2_last[1])/M,

                       (mH*vH1_last[2]+mH*vH2_last[2]+mCl*vCl1_last[2]+mCl*vCl2_last[2])/M)

                vH1_last  = v_sub(vH1_last,  vcm)

                vH2_last  = v_sub(vH2_last,  vcm)

                vCl1_last = v_sub(vCl1_last, vcm)

                vCl2_last = v_sub(vCl2_last, vcm)

            keH1  = ke_from_v(args.mass_H,  vH1_last)

            keH2  = ke_from_v(args.mass_H,  vH2_last)

            keCl1 = ke_from_v(args.mass_Cl, vCl1_last)

            keCl2 = ke_from_v(args.mass_Cl, vCl2_last)

            ker   = keH1 + keH2 + keCl1 + keCl2



        # ---------- 角度：最后 N 帧步进速度平均方向 ----------

        theta_HH = theta_ClCl = None

        if len(frames) >= 2:

            steps = min(args.angle_window-1, len(frames)-1)

            steps = max(1, steps)

            start = len(frames) - 1 - steps

            vH1_acc=vH2_acc=vCl1_acc=vCl2_acc=(0.0,0.0,0.0)

            for k in range(start, len(frames)-1):

                b0, b1 = frames[k], frames[k+1]

                vh1 = v_diff(b0,b1,h_idx[0],dt_s)

                vh2 = v_diff(b0,b1,h_idx[1],dt_s)

                vc1 = v_diff(b0,b1,cl_idx[0],dt_s)

                vc2 = v_diff(b0,b1,cl_idx[1],dt_s)

                if args.remove_com:

                    mH  = args.mass_H * AMU; mCl = args.mass_Cl * AMU; M = 2*mH+2*mCl

                    vcm = ((mH*vh1[0]+mH*vh2[0]+mCl*vc1[0]+mCl*vc2[0])/M,

                           (mH*vh1[1]+mH*vh2[1]+mCl*vc1[1]+mCl*vc2[1])/M,

                           (mH*vh1[2]+mH*vh2[2]+mCl*vc1[2]+mCl*vc2[2])/M)

                    vh1 = v_sub(vh1, vcm); vh2 = v_sub(vh2, vcm)

                    vc1 = v_sub(vc1, vcm); vc2 = v_sub(vc2, vcm)

                vH1_acc  = v_add(vH1_acc,  vh1)

                vH2_acc  = v_add(vH2_acc,  vh2)

                vCl1_acc = v_add(vCl1_acc, vc1)

                vCl2_acc = v_add(vCl2_acc, vc2)

            inv = 1.0/steps

            theta_HH   = angle_between(v_scale(vH1_acc,inv),  v_scale(vH2_acc,inv))

            theta_ClCl = angle_between(v_scale(vCl1_acc,inv), v_scale(vCl2_acc,inv))



        # —— 输出、入库（纳入统计的条目）——

        L_KER.append(ker); L_H1.append(keH1); L_H2.append(keH2); L_Cl1.append(keCl1); L_Cl2.append(keCl2)

        if theta_HH   is not None: L_tHH.append(theta_HH)

        if theta_ClCl is not None: L_tClCl.append(theta_ClCl)

        nOK += 1



        thh = f"{theta_HH:7.2f}" if theta_HH is not None else "   NA  "

        tcc = f"{theta_ClCl:7.2f}" if theta_ClCl is not None else "   NA  "

        print(f"{args.prefix}{idx:04d}  t_last={t_last_fs:7.3f} fs  |  KER={ker:10.4f} eV  |  H1={keH1:9.4f}  H2={keH2:9.4f}  Cl1={keCl1:9.4f}  Cl2={keCl2:9.4f}  |  theta_HH={thh}°  theta_ClCl={tcc}°")

        out_lines.append(f"{idx}\tOK\t{t_last_fs:.6f}\t{ker:.10f}\t{keH1:.10f}\t{keH2:.10f}\t{keCl1:.10f}\t{keCl2:.10f}\t{theta_HH if theta_HH is not None else 'NA'}\t{theta_ClCl if theta_ClCl is not None else 'NA'}")



    # ------- 统计并打印 -------

    def pr_stats(name, vals, nd=4):

        mu, sd, vmin, vmax, n = stats(vals)

        rng = f"[{fmt(vmin, nd)} , {fmt(vmax, nd)}]" if n>0 else "NA"

        print(f"{name:16s}  avg={fmt(mu, nd)}   std={fmt(sd, nd)}   range={rng}   N={n}")



    print("----- Summary statistics over valid (full-length) entries -----")

    pr_stats("KER (eV)",          L_KER,  4)

    pr_stats("KE_H1 (eV)",        L_H1,   4)

    pr_stats("KE_H2 (eV)",        L_H2,   4)

    pr_stats("KE_Cl1 (eV)",       L_Cl1,  4)

    pr_stats("KE_Cl2 (eV)",       L_Cl2,  4)

    pr_stats("theta_HH (deg)",    L_tHH,  2)

    pr_stats("theta_ClCl (deg)",  L_tClCl,2)

    print(f"Counts: OK={nOK}, skipped(incomplete)={nSkip}")



    # ------- 在 TSV 末尾追加 MEAN/STD/MIN/MAX -------

    def stat_tuple(vals): return stats(vals)[:4]  # (mean,std,min,max)

    mean_row = ["MEAN","-","-",

        fmt(stat_tuple(L_KER)[0],4),  fmt(stat_tuple(L_H1)[0],4),  fmt(stat_tuple(L_H2)[0],4),

        fmt(stat_tuple(L_Cl1)[0],4),  fmt(stat_tuple(L_Cl2)[0],4),

        fmt(stat_tuple(L_tHH)[0],2),  fmt(stat_tuple(L_tClCl)[0],2)]
    std_row  = ["STD","-","-",
        fmt(stat_tuple(L_KER)[1],4),  fmt(stat_tuple(L_H1)[1],4),  fmt(stat_tuple(L_H2)[1],4),
        fmt(stat_tuple(L_Cl1)[1],4),  fmt(stat_tuple(L_Cl2)[1],4),
        fmt(stat_tuple(L_tHH)[1],2),  fmt(stat_tuple(L_tClCl)[1],2)]
    min_row  = ["MIN","-","-",
        fmt(stat_tuple(L_KER)[2],4),  fmt(stat_tuple(L_H1)[2],4),  fmt(stat_tuple(L_H2)[2],4),
        fmt(stat_tuple(L_Cl1)[2],4),  fmt(stat_tuple(L_Cl2)[2],4),
        fmt(stat_tuple(L_tHH)[2],2),  fmt(stat_tuple(L_tClCl)[2],2)]
    max_row  = ["MAX","-","-",
        fmt(stat_tuple(L_KER)[3],4),  fmt(stat_tuple(L_H1)[3],4),  fmt(stat_tuple(L_H2)[3],4),
        fmt(stat_tuple(L_Cl1)[3],4),  fmt(stat_tuple(L_Cl2)[3],4),
        fmt(stat_tuple(L_tHH)[3],2),  fmt(stat_tuple(L_tClCl)[3],2)]
    out_lines += ["\t".join(mean_row), "\t".join(std_row), "\t".join(min_row), "\t".join(max_row)]
    with open(args.out, "w", encoding="utf-8") as fo:
        fo.write("\n".join(out_lines) + "\n")
    print(f"Results saved to: {args.out}")
if __name__ == "__main__":
    main()
