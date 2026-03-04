#!/usr/bin/env bash
set -euo pipefail

# ===== 你在这里改 =====
N0=32.0          # 体系初始电子数
MONOMERS=4       # 单元数：单体=1，二聚体=2，三聚体=3，四聚体=4
PATTERN='INFOLINE:        2.5000000000000000E+001'  # 保留：匹配的 INFOLINE 时间点
# =====================

read -rp "Start run number: " START
read -rp "End run number: " END

OUTCSV="ionization_${START}_${END}.csv"
echo "run,dir,status,time_fs,remaining_e,ionized_e,ionized_per_monomer,energy" > "$OUTCSV"

for ((n=START; n<=END; n++)); do
  dir=""

  cands=(
    "tdl${n}"
    "tdl$(printf "%03d" "$n")"
    "tdl$(printf "%04d" "$n")"
    "tdl$(printf "%05d" "$n")"
  )

  for cand in "${cands[@]}"; do
    if [[ -d "$cand" ]]; then
      dir="$cand"
      break
    fi
  done

  if [[ -z "$dir" ]]; then
    echo "tdl${n}: NO_DIR"
    echo "${n},,NO_DIR,,,," >> "$OUTCSV"
    continue
  fi

  mon="${dir}/monitor.out"
  if [[ ! -f "$mon" ]]; then
    echo "${dir}: NO_MONITOR"
    echo "${n},${dir},NO_MONITOR,,,," >> "$OUTCSV"
    continue
  fi

  line="$(grep -F "$PATTERN" "$mon" | tail -n 1 || true)"
  if [[ -z "$line" ]]; then
    echo "${dir}: NO_MATCH"
    echo "${n},${dir},NO_MATCH,,,," >> "$OUTCSV"
    continue
  fi

  # 提取 time($2), remaining($3), energy($4)
  # 计算 ionized = N0 - remaining
  # 计算 per-monomer = ionized / MONOMERS
  read -r time_fs remaining_e energy_e ionized_e ion_per <<< "$(
    awk -v n0="$N0" -v nm="$MONOMERS" '{
      ion = n0 - $3
      per = ion / nm
      # time remaining energy ionized per
      printf "%s %s %s %.10f %.10f", $2, $3, $4, ion, per
    }' <<< "$line"
  )"

  printf "%s: time=%s  remaining=%s  ionized=%s  per=%s\n" "$dir" "$time_fs" "$remaining_e" "$ionized_e" "$ion_per"
  echo "${n},${dir},OK,${time_fs},${remaining_e},${ionized_e},${ion_per},${energy_e}" >> "$OUTCSV"
done

echo "CSV written: ${OUTCSV}"
