#!/usr/bin/env bash

set -u



RUNS=(

  95 97 103 112 118 125 127 136 139 143 147 149 158 160 163 181 184 185 189 194 198 199 204 217 221 222 225 226 231 233 234 235 239 240 242 247 253 258 259 261 264 265 269 276 277 282 283 284 289 290
)



PATTERN='INFOLINE:        1.3000000000000000E+001'



for n in "${RUNS[@]}"; do

  # try common folder names: tdl91 / tdl0091 / tdl00091 ...

  dir=""
)
  for cand in "tdl${n}" "tdl$(printf "%03d" "$n")" "tdl$(printf "%04d" "$n")" "tdl$(printf "%05d" "$n")"; do
    if [[ -d "$cand" ]]; then
      dir="$cand"
      break
    fi
  done
  if [[ -z "$dir" ]]; then
    echo "==== tdl${n} : NO_DIR ===="
    continue
  fi
  if [[ ! -f "${dir}/monitor.out" ]]; then
    echo "==== ${dir} : NO_MONITOR ===="
    continue
  fi
  echo "==== ${dir} ===="
  grep -n "$PATTERN" "${dir}/monitor.out" || echo "(no match)"
done
