#!/usr/bin/env bash
set -euo pipefail

nevents="1000"

width_values=(5 10 15 20 30 40 50 60 80 100)
lambda_values=(10 15 20 25 30 40 50 60 75 100)
mq_values=(2.5 5 7.5 10 15 20 25 30 40 50)

scan_card="hvZ_var.cmnd"

format_value() {
  printf '%s' "$1" | tr '.' 'p'
}

double_value() {
  case "$1" in
    *.*) printf '%s' "$1" ;;
    *) printf '%s.0' "$1" ;;
  esac
}

format_seconds() {
  local seconds=$1
  local hours=$((seconds / 3600))
  local minutes=$(((seconds % 3600) / 60))
  local secs=$((seconds % 60))

  printf '%02d:%02d:%02d' "$hours" "$minutes" "$secs"
}

show_time_estimate() {
  local completed=$1
  local total=$2
  local start_time=$3
  local now elapsed estimated_total remaining

  now=$(date +%s)
  elapsed=$((now - start_time))
  estimated_total=$((elapsed * total / completed))
  remaining=$((estimated_total - elapsed))

  printf '\rCompleted %d/%d | elapsed %s | estimated total %s | remaining %s' \
    "$completed" "$total" \
    "$(format_seconds "$elapsed")" \
    "$(format_seconds "$estimated_total")" \
    "$(format_seconds "$remaining")"
}

total_runs=$((${#width_values[@]} * ${#lambda_values[@]} * ${#mq_values[@]}))
completed_runs=0
start_time=$(date +%s)

for width in "${width_values[@]}"; do
  for lambda in "${lambda_values[@]}"; do
    for mq in "${mq_values[@]}"; do
      cat > "$scan_card" <<EOF
4900023:mWidth = $(double_value "$width")
HiddenValley:Lambda = $(double_value "$lambda")
4900101:m0 = $(double_value "$mq")
EOF
      output="hvZ2_Lambda$(format_value "$lambda")_width$(format_value "$width")_mq$(format_value "$mq")"
      ./hiddenvalley -e "$nevents" -s 2 -o "$output" hvZ.cmnd "$scan_card"
      completed_runs=$((completed_runs + 1))
      show_time_estimate "$completed_runs" "$total_runs" "$start_time"
    done
  done
done

printf '\n'
