#!/usr/bin/env bash
#SBATCH --job-name=get_pdb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=short

set -euo pipefail

# --- modules & environment ---
module purge || true
module use /projects/SimBioSys/share/software/modulefiles
module load VMD/1.9.4a55
module load miniconda3/24.11.1 || true

ALLOSMOD_CONDA_ENV="${ALLOSMOD_CONDA_ENV:-/projects/SimBioSys/share/software/allosmod-env}"
echo "[get_pdb] Using ALLOSMOD_CONDA_ENV=$ALLOSMOD_CONDA_ENV"

# Initialize conda in non-interactive shells, then activate the env
if command -v conda >/dev/null 2>&1; then
  if [[ "$(type -t conda)" != "function" ]]; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]] && . "$CONDA_BASE/etc/profile.d/conda.sh"
  fi
  if ! conda activate "$ALLOSMOD_CONDA_ENV" 2>/dev/null; then
    export PATH="$ALLOSMOD_CONDA_ENV/bin:$PATH"
  fi
else
  export PATH="$ALLOSMOD_CONDA_ENV/bin:$PATH"
fi

# --- args ---
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <folder> <USER_ID>"
  exit 1
fi
folder="$1"
USER_ID="$2"

# --- .env loader (honors LOG_DIR/LOG_ROOT/OUTPUT_DIR/ALIGN_SEL) ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
try_source_env() {
  local envfile="$1"
  [[ -f "$envfile" ]] || return 1
  set -a
  # shellcheck disable=SC1090
  . "$envfile"
  set +a
  return 0
}
for cand in ".env" "$SCRIPT_DIR/.env" "$folder/.env" "$(dirname "$folder")/.env"; do
  try_source_env "$cand" && break || true
done

# --- directories (logs + outputs) ---
BASE_LOG_ROOT="${LOG_DIR:-${LOG_ROOT:-"$PWD/logs"}}"
USER_LOG_ROOT="${BASE_LOG_ROOT%/}/${USER_ID}"
ENSEMBLE_LOG_DIR="${USER_LOG_ROOT}/ensemble_modelling"
VMD_LOG_DIR="${ENSEMBLE_LOG_DIR}/VMD"
mkdir -p "$ENSEMBLE_LOG_DIR" "$VMD_LOG_DIR"

OUT_BASE="${OUTPUT_DIR:-${ENSEMBLE_LOG_DIR}}"
[[ "$OUT_BASE" = /* ]] || OUT_BASE="$PWD/${OUT_BASE}"
mkdir -p "$OUT_BASE"

SCRIPTS_BASE="${ENSEMBLE_SCRIPTS_DIR:-$SCRIPT_DIR}"

# Mirror stdout/err alongside Slurm logs
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.err" >&2)
else
  exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.err" >&2)
fi

# --- tool checks ---
for cmd in awk vmd pdb_tidy pdb_reres pdb_reatom pdb_selchain; do
  command -v "$cmd" >/dev/null || { echo "Missing required tool: $cmd"; exit 1; }
done

# --- scripts (absolute paths) ---
ALIGN_TRAJECTORY_SCRIPT="${SCRIPTS_BASE}/align_trajectory.tcl"
COMPUTE_GLYCAN_PARAMS_SCRIPT="${SCRIPTS_BASE}/compute_glycan_params.py"
EXTRACT_GLYCANS_AND_CHAINS_SCRIPT="${SCRIPTS_BASE}/extract_glycans_and_chains.tcl"

[[ -f "$ALIGN_TRAJECTORY_SCRIPT" ]] || { echo "Missing $ALIGN_TRAJECTORY_SCRIPT"; exit 1; }
[[ -f "$COMPUTE_GLYCAN_PARAMS_SCRIPT" ]] || { echo "Missing $COMPUTE_GLYCAN_PARAMS_SCRIPT"; exit 1; }
[[ -f "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" ]] || { echo "Missing $EXTRACT_GLYCANS_AND_CHAINS_SCRIPT"; exit 1; }

# Default alignment selection (override via .env ALIGN_SEL)
ALIGN_SEL="${ALIGN_SEL:-name CA}"

# Avoid reading site/user startup files in VMD
export VMDNOSTARTUP=1

# --- per-subfolder processing ---
for m in "$folder"/*; do
  [[ -d "$m" ]] || continue
  model_name="$(basename "$m")"

  # output + log dirs
  out_dir="${OUT_BASE%/}/${USER_ID}/${model_name}"
  model_vmd_log_dir="${VMD_LOG_DIR%/}/${model_name}"
  mkdir -p "$out_dir" "$model_vmd_log_dir"

  # inputs
  [[ -f "$m/input.dat" ]] || { echo "Skipping $model_name (no input.dat)"; continue; }
  NRUNS=$(awk -F= '/^NRUNS=/{match($2,/[0-9]+/,mm); print mm[0]; exit}' "$m/input.dat")
  [[ -n "${NRUNS:-}" ]] || { echo "Skipping $model_name (NRUNS not found)"; continue; }

  # glycan params
  read -r g1_start g1_end g2_start g2_end g3_start g3_end chain_len first_chain_len_plus_one \
    < <(python3 "$COMPUTE_GLYCAN_PARAMS_SCRIPT" "$m/align.ali" "$m/glyc.dat")

  for v in "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end"; do
    [[ "$v" =~ ^-?[0-9]+$ ]] || { echo "Bad glycan param '$v' for $model_name"; exit 1; }
  done

  out_file="${out_dir}/output_${model_name}.pdb"
  : > "$out_file"

  # --- per run ---
  for i in "${m%/}"/pred_dECALCrAS1000/*.pdb_*/; do   # ensure dirs
    run_id="$(basename "${i%/}")"                     # e.g., start.pdb_0
    filename="${i}pm.pdb.B99990001.pdb" 
    [[ -f "$filename" ]] || { echo "Missing $filename (run $i) — skipping"; continue; }

    (
      cd "$out_dir" || exit 1

      # Extract chains & glycans with VMD (logs to model_vmd_log_dir)
      vmd -dispdev text -eofexit \
          -e "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" -args \
          "$filename" "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end" \
          > "${model_vmd_log_dir}/vmd_run_${run_id}.log" 2>&1

      # Determine starting residue index for chain A
      starting_chain_number=$(
        pdb_tidy "$filename" \
          | pdb_selchain -A \
          | awk '/^(ATOM|HETATM)/{print substr($0,23,4)+0; exit}'
      )
      [[ -n "${starting_chain_number:-}" ]] || { rm -f C*.pdb || true; exit 0; }

      # Renumber pieces
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
      [[ ${#files[@]} -gt 0 ]] || { rm -f C*.pdb || true; exit 0; }

      # Concatenate, strip meta lines, re-atom index; NO per-run END
      awk '{
             sub(/\r$/, "")
             if ($0 ~ /^(CRYST1|REMARK|TER|END|ENDMDL)([[:space:]].*)?$/) next
             print
           }' "${files[@]}" \
        | pdb_reatom -1 >> "$out_file"

      rm -f C*.pdb || true
    )
  done

  # Single END at the end of the model file
  if [[ -s "$out_file" ]]; then
    printf 'END\n' >> "$out_file"
  fi

  # Align (selection defaults to CA)
  if [[ -s "$out_file" ]]; then
      export VMDNOSTARTUP=1
      export IN_PDB="$out_file"
      export OUT_PDB="${out_dir}/output_aligned.pdb"
      export ALIGN_SEL_ENV="${ALIGN_SEL:-name CA}"

      vmd -dispdev text -eofexit \
          -e "$ALIGN_TRAJECTORY_SCRIPT" \
          > "${model_vmd_log_dir}/align_trajectory_vmd.log" 2>&1
  fi
done

# best-effort cleanup if anything leaked in the submit dir
rm -f C*.pdb || true
