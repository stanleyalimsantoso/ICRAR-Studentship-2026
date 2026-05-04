#! /bin/bash

set -euo pipefail

python_venv="/home/stanley/icrar_env/bin/python"

remote_user="salimsan"
remote_host="ozstar.swin.edu.au"
remote_source="/fred/oz080/runvimfi.dir/tout.dat"

total_gas=0 #1 for total gas and 0 for h2gas
algo="cdktree" #cdktree or naive_search(do not choose naive search)

run=82
run_dir="/mnt/d/ICRAR Studentship Data/Run ${run}"

mkdir -p "$run_dir"

raw_dir="${run_dir}/raw"
pre_dir="${run_dir}/preprocessed"
log_dir="${run_dir}/logs"
gmc_dir="${run_dir}/gmc"
post_dir="${run_dir}/postprocessed"
output_dir="${run_dir}/output"
profiles_dir="${run_dir}/gmc_profiles"
plots_dir="${run_dir}/plots"
output_corrected_dir="${run_dir}/output_corrected"


mkdir -p "$raw_dir" "$pre_dir" "$log_dir" "$gmc_dir" "$post_dir" "$output_dir" "$profiles_dir" "$plots_dir" "$output_corrected_dir"

log_file="${log_dir}/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log_file") 2>&1

echo "Run dir: $run_dir"
echo "Log: $log_file"

echo "1) Checking for tout.dat in raw_dir"

if [ -f "${raw_dir}/tout.dat" ]; then
  echo "   Found existing ${raw_dir}/tout.dat"
  read -r -p "   tout.dat already exists. Overwrite from OzSTAR? [y/N]: " overwrite

  case "$overwrite" in
    y|Y)
      echo "   Overwriting tout.dat from OzSTAR"
      rsync -rvP --no-times --no-perms --no-owner --no-group --inplace -e "ssh" \
        "${remote_user}@${remote_host}:${remote_source}" \
        "${raw_dir}/tout.dat"
      ;;
    *)
      echo "   Keeping existing tout.dat; skipping pull"
      ;;
  esac
else
  echo "   Pulling tout.dat from OzSTAR"
  rsync -rvP --no-times --no-perms --no-owner --no-group --inplace -e "ssh" \
    "${remote_user}@${remote_host}:${remote_source}" \
    "${raw_dir}/tout.dat"
fi


max_timestamp="$(awk 'NR==1{print $3; exit}' "${raw_dir}/tout.dat")"
echo "max_timestamp: $max_timestamp"

echo "2) Classifying GMCs with AGB stars"

$python_venv "/mnt/d/ICRAR Studentship Data/parse_tout.py" \
  --input "${raw_dir}/tout.dat" \
  --preprocessed "$pre_dir" \
  --max-timestamp "$max_timestamp" \
  $( [ "$total_gas" -eq 1 ] && printf -- "--total-gas" ) \
  --run "$run" \
  --postprocessed "$post_dir" \
  --output "$output_dir" \
  --gmc-dir "$gmc_dir" \
  --plots-dir "$plots_dir"

echo "3) Building GMC profiles"

$python_venv "/mnt/d/ICRAR Studentship Data/make_gmc_profiles.py" \
  --output-dir "$output_dir" \
  --profiles-root "$profiles_dir"

best_gmc_file="$(
  awk -F': ' '/Largest AGB mass fraction file/ {print $2}' "$log_file" | tail -n 1
)"

echo "Using GMC file: $best_gmc_file"

echo "4) Creating GMC data with environment"

$python_venv "/mnt/d/ICRAR Studentship Data/reparse_tout.py" \
  --input "$best_gmc_file" \
  --output "$output_corrected_dir" \
  --run "$run"





