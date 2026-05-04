#! /bin/bash

set -euo pipefail

python_venv="/home/stanley/icrar_env/bin/python"

remote_user="salimsan"
remote_host="ozstar.swin.edu.au"
remote_source="/fred/oz080/runvimf301.dir/tout.dat"

total_gas=0 #1 for total gas and 0 for h2gas
algo="cdktree" #cdktree or naive_search(do not choose naive search)

run=81
run_dir="/mnt/d/ICRAR Studentship Data/Run ${run}"
gc_dir="${run_dir}/GC Data"
preprocessed_dir="${gc_dir}/preprocessed"
plots_dir="${gc_dir}/plots"
raw_dir="${gc_dir}/raw"


mkdir -p "$gc_dir" "$preprocessed_dir" "$plots_dir" "$raw_dir"

# echo "1) Pulling tout.dat from OzSTAR"
# rsync -rvP --no-times --no-perms --no-owner --no-group --inplace -e "ssh" \
#   "${remote_user}@${remote_host}:${remote_source}" \
#   "${raw_dir}/tout.dat"

max_timestamp="$(awk 'NR==2{print $1; exit}' "${raw_dir}/tout.dat")"
max_timestamp=10
echo "max_timestamp: $max_timestamp"

echo "2) Creating Plots"

$python_venv "/mnt/d/ICRAR Studentship Data/parse_tout_gc.py" \
  --input "${raw_dir}/tout.dat" \
  --preprocessed "$preprocessed_dir" \
  --max-timestamp "$max_timestamp" \
  --run "$run" \
  --plots-dir "$plots_dir"






