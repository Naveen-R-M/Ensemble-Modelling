#!/usr/bin/env bash
# lib/common.sh - Common functions and environment setup

# --- Global variables ---
declare -g SCRIPT_DIR
declare -g ALLOSMOD_CONDA_ENV  
declare -g BASE_LOG_ROOT USER_LOG_ROOT ENSEMBLE_LOG_DIR VMD_LOG_DIR OUT_BASE SCRIPTS_BASE

# --- Logging functions ---
log_info() {
    echo "[info] $*" >&2
}

log_warn() {
    echo "[warn] $*" >&2
}

log_error() {
    echo "[error] $*" >&2
}

# --- Environment setup ---
setup_modules() {
    log_info "Setting up modules and environment"
    
    module purge || true
    module use /projects/SimBioSys/share/software/modulefiles
    module load VMD/1.9.4a55
    module load miniconda3/24.11.1 || true

    # Avoid reading site/user startup files in VMD for all invocations
    export VMDNOSTARTUP=1

    ALLOSMOD_CONDA_ENV="${ALLOSMOD_CONDA_ENV:-/projects/SimBioSys/share/software/allosmod-env}"
    log_info "Using ALLOSMOD_CONDA_ENV=$ALLOSMOD_CONDA_ENV"

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
}

# --- Environment file loader ---
try_source_env() {
    local envfile="$1"
    [[ -f "$envfile" ]] || return 1
    set -a
    # shellcheck disable=SC1090
    . "$envfile"
    set +a
    return 0
}

load_env_files() {
    local folder="$1"
    
    for cand in ".env" "$SCRIPT_DIR/.env" "$folder/.env" "$(dirname "$folder")/.env"; do
        try_source_env "$cand" && break || true
    done
}

# --- Directory setup ---
setup_directories() {
    local user_id="$1"
    
    BASE_LOG_ROOT="${LOG_DIR:-${LOG_ROOT:-"$PWD/logs"}}"
    USER_LOG_ROOT="${BASE_LOG_ROOT%/}/${user_id}"
    ENSEMBLE_LOG_DIR="${USER_LOG_ROOT}/ensemble_modelling"
    VMD_LOG_DIR="${ENSEMBLE_LOG_DIR}/VMD"
    mkdir -p "$ENSEMBLE_LOG_DIR" "$VMD_LOG_DIR"

    OUT_BASE="${OUTPUT_DIR:-${ENSEMBLE_LOG_DIR}}"
    [[ "$OUT_BASE" = /* ]] || OUT_BASE="$PWD/${OUT_BASE}"
    mkdir -p "$OUT_BASE"

    SCRIPTS_BASE="${ENSEMBLE_SCRIPTS_DIR:-$SCRIPT_DIR}"
}

# --- Logging redirection ---
setup_logging() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_${SLURM_JOB_ID}.err" >&2)
    else
        exec > >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.out") 2> >(tee -a "${ENSEMBLE_LOG_DIR}/get_pdb_standalone.err" >&2)
    fi
}

# --- Tool validation ---
check_required_tools() {
    local missing_tools=()
    
    for cmd in awk vmd pdb_tidy pdb_reres pdb_reatom pdb_selchain; do
        if ! command -v "$cmd" >/dev/null; then
            missing_tools+=("$cmd")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        exit 1
    fi
}

# --- Script validation ---
validate_scripts() {
    local align_trajectory_script="${SCRIPTS_BASE}/align_trajectory.tcl"
    local extract_glycans_script="${SCRIPTS_BASE}/extract_glycans_and_chains.tcl"
    
    [[ -f "$align_trajectory_script" ]] || { log_error "Missing $align_trajectory_script"; exit 1; }
    [[ -f "$extract_glycans_script" ]] || { log_error "Missing $extract_glycans_script"; exit 1; }
    
    # Export for use in other modules
    export ALIGN_TRAJECTORY_SCRIPT="$align_trajectory_script"
    export EXTRACT_GLYCANS_AND_CHAINS_SCRIPT="$extract_glycans_script"
}

# --- Input validation ---
parse_input_dat() {
    local input_dat="$1"
    
    [[ -f "$input_dat" ]] || return 1
    
    local nruns
    nruns=$(awk -F= '/^NRUNS=/{match($2,/[0-9]+/,mm); print mm[0]; exit}' "$input_dat")
    [[ -n "${nruns:-}" ]] || return 1
    
    echo "$nruns"
}

# --- PDB prefix detection ---
detect_pdb_prefix() {
    local pred_dir="$1"
    
    for candidate_dir in "${pred_dir}"/pred_dECALCrAS1000/*.pdb_0/; do
        if [[ -d "$candidate_dir" ]]; then
            local candidate_name
            candidate_name=$(basename "$candidate_dir")
            echo "${candidate_name%_0}"  # Remove _0 suffix
            return 0
        fi
    done
    
    return 1
}

# --- Validation helpers ---
validate_integer() {
    local value="$1"
    local name="$2"
    
    [[ "$value" =~ ^-?[0-9]+$ ]] || { log_error "Bad integer for $name: '$value'"; exit 1; }
}

# --- File cleanup ---
cleanup_intermediate_files() {
    rm -f C*.pdb *_renum.pdb || true
}
