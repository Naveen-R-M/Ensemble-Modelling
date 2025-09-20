#!/usr/bin/env bash
#SBATCH --job-name=get_pdb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# --- args ---
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <folder> <USER_ID>"
  exit 1
fi
folder="$1"
USER_ID="$2"

# --- load .env (only LOG_DIR / LOG_ROOT are honored) ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
try_source_env() {
  local envfile="$1"
  [[ -f "$envfile" ]] || return 1
  set -a
  # shellcheck disable=SC1090
  . "$envfile"
  set +a
  echo "Loaded .env from: $envfile"
  return 0
}
ENV_CANDIDATES=(
  ".env"
  "$SCRIPT_DIR/.env"
  "$folder/.env"
  "$(dirname "$folder")/.env"
)
for cand in "${ENV_CANDIDATES[@]}"; do
  try_source_env "$cand" && break
done

BASE_LOG_ROOT="${LOG_DIR:-${LOG_ROOT:-"$PWD/logs"}}"
USER_LOG_ROOT="${BASE_LOG_ROOT%/}/${USER_ID}"
ENSEMBLE_LOG_DIR="${USER_LOG_ROOT}/ensemble_modelling"
mkdir -p "$ENSEMBLE_LOG_DIR"

# Mirror stdout/err into ensemble logs (besides Slurm %x-%j files)
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.err" >&2)
else
  exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.err" >&2)
fi

# --- sanity checks ---
command -v awk >/dev/null || { echo "awk not found"; exit 1; }
command -v vmd  >/dev/null || { echo "vmd not found"; exit 1; }
command -v pdb_tidy  >/dev/null || { echo "pdb_tidy (pdb-tools) not found"; exit 1; }
command -v pdb_reres >/dev/null || { echo "pdb_reres (pdb-tools) not found"; exit 1; }
command -v pdb_reatom >/dev/null || { echo "pdb_reatom (pdb-tools) not found"; exit 1; }
command -v pdb_selchain >/dev/null || { echo "pdb_selchain (pdb-tools) not found"; exit 1; }

ALIGN_TRAJECTORY_SCRIPT="$SCRIPT_DIR/align_trajectory.tcl"
COMPUTE_GLYCAN_PARAMS_SCRIPT="$SCRIPT_DIR/compute_glycan_params.py"
EXTRACT_GLYCANS_AND_CHAINS_SCRIPT="$SCRIPT_DIR/extract_glycans_and_chains.tcl"
[[ -f "$ALIGN_TRAJECTORY_SCRIPT" ]] || { echo "Missing $ALIGN_TRAJECTORY_SCRIPT"; exit 1; }
[[ -f "$COMPUTE_GLYCAN_PARAMS_SCRIPT" ]] || { echo "Missing $COMPUTE_GLYCAN_PARAMS_SCRIPT"; exit 1; }
[[ -f "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" ]] || { echo "Missing $EXTRACT_GLYCANS_AND_CHAINS_SCRIPT"; exit 1; }

# Each subfolder processed, outputs go under ENSEMBLE_LOG_DIR/<subfolder>/
for m in "$folder"/*; do
  [[ -d "$m" ]] || continue
  model_name="$(basename "$m")"
  out_dir="${ENSEMBLE_LOG_DIR}/${model_name}"
  mkdir -p "$out_dir"

  NRUNS=$(awk -F= '/^NRUNS=/{match($2,/[0-9]+/,mm); print mm[0]; exit}' "$m/input.dat")
  [[ -n "${NRUNS:-}" ]] || { echo "NRUNS not found in $m/input.dat"; exit 1; }

  read -r g1_start g1_end g2_start g2_end g3_start g3_end chain_len first_chain_len_plus_one \
    < <(python3 "$COMPUTE_GLYCAN_PARAMS_SCRIPT" "$m/align.ali" "$m/glyc.dat")

  for v in "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end"; do
    [[ "$v" =~ ^-?[0-9]+$ ]] || { echo "Bad glycan param: '$v'"; exit 1; }
  done

  out_file="${out_dir}/output_${model_name}.pdb"
  : > "$out_file"

  for i in $(seq 0 $((NRUNS - 1))); do
    filename="${m}/pred_dECALCrAS1000/siv.pdb_${i}/pm.pdb.B99990001.pdb"
    echo "Processing $filename"
    if [[ ! -f "$filename" ]]; then
      echo "Warning: missing $filename — skipping run $i."
      continue
    fi
    
    # Work inside the model's output dir so CHA.pdb, CAR*.pdb land here
    pushd "$out_dir" >/dev/null

    vmd -dispdev text -eofexit -log "${out_dir}/vmd_run_${i}.log" \
      -e "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" -args \
      "$filename" "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end"

    starting_chain_number=$(
      pdb_tidy "$filename" \
        | pdb_selchain -A \
        | awk '/^(ATOM|HETATM)/{print substr($0,23,4)+0; exit}'
    )
    if [[ -z "${starting_chain_number:-}" ]]; then
      echo "Warning: no ATOM/HETATM records for chain A in $filename — skipping run $i."
      continue
    fi

    if [[ "$chain_len" -eq 3 ]]; then
      pdb_reres -"$starting_chain_number" CHA.pdb > CHA_renum.pdb
      pdb_reres -"$starting_chain_number" CHB.pdb > CHB_renum.pdb
      pdb_reres -"$starting_chain_number" CHC.pdb > CHC_renum.pdb
      pdb_reres -1 CAR1.pdb > CAR1_renum.pdb
      pdb_reres -1 CAR2.pdb > CAR2_renum.pdb
      pdb_reres -1 CAR3.pdb > CAR3_renum.pdb
      concat_files=(CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb)
    else
      pdb_reres -"$starting_chain_number"      CHA.pdb > CHA_renum.pdb
      pdb_reres -"$first_chain_len_plus_one"   CHB.pdb > CHB_renum.pdb
      pdb_reres -"$starting_chain_number"      CHC.pdb > CHC_renum.pdb
      pdb_reres -"$first_chain_len_plus_one"   CHD.pdb > CHD_renum.pdb
      pdb_reres -"$starting_chain_number"      CHE.pdb > CHE_renum.pdb
      pdb_reres -"$first_chain_len_plus_one"   CHF.pdb > CHF_renum.pdb
      pdb_reres -1 CAR1.pdb > CAR1_renum.pdb
      pdb_reres -1 CAR2.pdb > CAR2_renum.pdb
      pdb_reres -1 CAR3.pdb > CAR3_renum.pdb
      concat_files=(CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb CHD_renum.pdb CHE_renum.pdb CHF_renum.pdb CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb)
    fi

    files=()
    for f in "${concat_files[@]}"; do
      [[ -s "$f" ]] && files+=("$f")
    done
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "Warning: no input pieces produced for ${model_name}, run $i — skipping."
      continue
    fi

    {
      awk '{
            sub(/\r$/, "");
            if ($0 ~ /^(CRYST1|REMARK|TER|END|ENDMDL)([[:space:]].*)?$/) next;
            print
          }' "${files[@]}" \
      | pdb_reatom -1
      echo 'END'
    } >> "$out_file"

    rm -f C*.pdb
  done

  echo "Aligning trajectory for ${model_name}..."
  if [[ ! -s "$out_file" ]]; then
    echo "Warning: $out_file is empty or missing — skipping alignment for ${model_name}."
  else
    vmd -dispdev text -eofexit \
        -e "$ALIGN_TRAJECTORY_SCRIPT" -args \
        "$out_file" "${out_dir}/output_aligned.pdb" "protein and backbone"
  fi
done

rm -f C*.pdb || true
