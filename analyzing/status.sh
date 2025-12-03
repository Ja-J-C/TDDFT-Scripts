#!/bin/bash



# stat.sh  ── 当前帧时间(由 iter 换算) + 碎片状态



start=1

end=1350

base_dir="."



# 每个迭代步长对应的 fs

dt_fs=0.001



#================ 工具函数 ==================



# awk 片段：基于 2.0 Å 阈值用并查集聚类，输出碎片化学式

awk_frag(){

cat <<'AWK'

BEGIN{th=2.0; th2=th*th}

function dist2(i,j){dx=x[i]-x[j];dy=y[i]-y[j];dz=z[i]-z[j]; return dx*dx+dy*dy+dz*dz}

# 并查集

function find(a){return parent[a]==a?a:parent[a]=find(parent[a])}

{

  elem[NR]=$1; x[NR]=$2; y[NR]=$3; z[NR]=$4; parent[NR]=NR

}

END{

  n=NR

  for(i=1;i<=n;i++)

    for(j=i+1;j<=n;j++)

      if(dist2(i,j)<th2){

        ri=find(i); rj=find(j);

        if(ri!=rj) parent[rj]=ri

      }



  # 汇总片段

  for(i=1;i<=n;i++){

    root=find(i)

    comp[root]=comp[root] elem[i] ","

  }



  # 组装化学式

  out=""; first=1

  for(root in comp){

    split(comp[root],arr,","); delete count

    for(k in arr){

      e=arr[k]

      if(e!=""){count[e]++}

    }

    formula=""

    for(e in count){formula=formula e count[e]}

    if(!first) out=out"+"

    out=out formula

    first=0

  }

  print out

}

AWK

}



#===========================================



for ((idx=start; idx<=end; idx++)); do

  folder="td${idx}"

  path="$base_dir/$folder"

  traj="$path/trajectory.xyz"

  control="$path/control.inp"

  monitor="$path/monitor.out"



  [[ ! -d "$path" ]] && continue



  reason=""

  complete=false

  prod=""

  time_fs="unknown"



  # -------- 读取最后一帧 & 当前碎片 / 时间 --------

  if [[ -f "$traj" ]]; then

    natoms=$(head -n 1 "$traj")

    (( natoms<1 )) && natoms=0



    if (( natoms > 0 )); then

      frame_size=$((natoms+2))



      # 抽取最后一帧（natoms 行 + 注释行）

      tail -n "$frame_size" "$traj" > __last_frame.tmp

      # 第 3..(2+natoms) 行是坐标

      tail -n "$natoms" __last_frame.tmp > __last.xyz



      # 当前碎片化学式（不管是否真正“碎裂”，都给出当前状态）

      prod=$(awk -f <(awk_frag) __last.xyz)



      # 判断是否存在 >5 Å 的原子对 → 视为“已经碎裂”，用于 status / complete

      if awk '{x[NR]=$2;y[NR]=$3;z[NR]=$4}

              END{

                for(i=1;i<=NR;i++)

                  for(j=i+1;j<=NR;j++){

                    d=(x[i]-x[j])^2+(y[i]-y[j])^2+(z[i]-z[j])^2

                    if(d>25){exit 0}

                  }

                exit 1

              }' __last.xyz; then

        complete=true

        reason="fragmented"

      fi



      # 从最后一帧的注释行中解析 iter，然后用 iter*dt_fs 算 time_fs

      last_header=$(sed -n '2p' __last_frame.tmp)

      # 匹配形如 "iter =430500" 或 "iter = 430500"

      iter_val=$(awk 'match($0,/iter[[:space:]]*=[[:space:]]*([0-9]+)/,a){print a[1]}' <<< "$last_header")



      if [[ -n "$iter_val" ]]; then

        # time_fs = iter * dt_fs

        time_fs=$(awk -v it="$iter_val" -v dt="$dt_fs" 'BEGIN{printf "%.3f", it*dt}')

      else

        time_fs="unknown"

      fi



      rm -f __last.xyz __last_frame.tmp

    fi

  else

    time_fs="0.0"

    prod=""

  fi



  # -------- 条件: Run finished --------

  if [[ "$complete" == false && -f "$monitor" ]]; then

     if tail -n1 "$monitor" | grep -q "Run finished:"; then

       complete=true

       reason="finished"

     fi
  fi
  # 温度信息（可选）
  temp=""
  if [[ -f "$control" ]]; then
    temp=$(grep -i "temperature_ions" "$control" | head -n1 | tr -d '[:space:]')
  fi
  [[ -n "$temp" ]] && temp=" [$temp]"
  # 状态字符串
  status="$reason"
  [[ -z "$status" ]] && status="running"
  # -------- 不管 complete 与否，一律输出当前 fs + 碎片 --------
  line="$folder  t=${time_fs}fs"
  [[ -n "$prod" ]] && line+="  fragments=$prod"
  line+="  status=$status$temp"
  echo "$line"
done
