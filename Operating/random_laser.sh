#!/usr/bin/env bash

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=1

#SBATCH --time=1:59:59



# ===== 参数：目标向量长度 =====

R=9        # |S| = 9  →  Sx^2 + Sy^2 + Sz^2 = 81



# 30 个不重复随机种子

readarray -t SEEDS < <(shuf -i 1-10000 -n 50)



# 生成 |S|=R 的球面均匀随机向量

gen_vec() {

  local seed="$1" R="$2"

  awk -v s="$seed" -v R="$R" 'BEGIN{

    srand(s)

    u = 2*rand()-1

    pi = 4*atan2(1,1)

    phi = 2*pi*rand()

    t = sqrt(1 - u*u)

    sx = R * t * cos(phi)

    sy = R * t * sin(phi)

    sz = R * u

    printf("%.6f %.6f %.6f\n", sx, sy, sz)

  }'

}



# 30 组：编号 61..90（需要别的区间就改这里）

for i in $(seq 241 290); do

  dir="tdl$i"

  mkdir -p "$dir"

  cd "$dir" || exit 1



  cp ../td/control.inp .
  cp ../td/job_script.sh .
  idx=$((i-241))                  # 让 idx 从 0..29
  seed="${SEEDS[$idx]}"
  read -r sx sy sz < <(gen_vec "$seed" "$R")
  # 允许等号两边有空格并去掉 CRLF
  sed -i -e 's/\r$//' control.inp
  sed -i -E "s/^[[:space:]]*N_time_steps[[:space:]]*=.*/N_time_steps=120000/" control.inp
  sed -i -E "s/^[[:space:]]*e_laser1[[:space:]]*=.*/e_laser1=${sx}/"         control.inp
  sed -i -E "s/^[[:space:]]*e_laser2[[:space:]]*=.*/e_laser2=${sy}/"         control.inp
  sed -i -E "s/^[[:space:]]*e_laser3[[:space:]]*=.*/e_laser3=${sz}/"         control.inp
  sed -i -E "s/^[[:space:]]*ion_velocity_init_seed[[:space:]]*=.*/ion_velocity_init_seed=${seed}/" control.inp
  sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=(HCl)2_tdl${i}|" job_script.sh
  sbatch job_script.sh
  cd ..
done
