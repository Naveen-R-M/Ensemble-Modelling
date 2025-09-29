#!/usr/bin/env bash
#SBATCH --job-name=get_pdb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=short

# Main orchestrator script - simplified and modular

# --- Initialize ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

# Source modules
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"
# shellcheck source=lib/alignment_parser.sh  
source "$SCRIPT_DIR/lib/alignment_parser.sh"
# shellcheck source=lib/glycan_calculator.sh
source "$SCRIPT_DIR/lib/glycan_calculator.sh"
# shellcheck source=lib/pdb_processor.sh
source "$SCRIPT_DIR/lib/pdb_processor.sh"

# Load configuration
# shellcheck source=config/defaults.conf
source "$SCRIPT_DIR/config/defaults.conf"

# --- Argument validation ---
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <folder> <USER_ID>"
    exit 1
fi

folder="$1"
USER_ID="$2"

# --- Environment setup ---
log_info "Starting ensemble modeling pipeline"
log_info "Processing folder: $folder"
log_info "User ID: $USER_ID"

setup_modules
load_env_files "$folder"
setup_directories "$USER_ID"
setup_logging
check_required_tools
validate_scripts

# --- Process each model subfolder ---
for m in "$folder"/*; do
    [[ -d "$m" ]] || continue
    model_name="$(basename "$m")"
    log_info "Processing model: $model_name"

    # Setup output directories
    out_dir="${OUT_BASE%/}/${USER_ID}/${model_name}"
    model_vmd_log_dir="${VMD_LOG_DIR%/}/${model_name}"
    mkdir -p "$out_dir" "$model_vmd_log_dir"

    # Parse input.dat
    NRUNS=$(parse_input_dat "$m/$INPUT_FILE")
    if [[ -z "$NRUNS" ]]; then
        log_warn "Skipping $model_name (no valid $INPUT_FILE or NRUNS)"
        continue
    fi
    log_info "Found NRUNS=$NRUNS for model $model_name"

    # Detect PDB prefix pattern
    pdb_prefix=$(detect_pdb_prefix "$m")
    if [[ -z "$pdb_prefix" ]]; then
        log_error "Could not detect PDB prefix pattern in $m/$PDB_DIR_PATTERN/"
        continue
    fi
    log_info "Detected PDB prefix pattern: ${pdb_prefix}_*"

    # Process alignment file - declare array first
    declare -A alignment_params=()
    if ! process_alignment "$m/$ALIGNMENT_FILE" alignment_params; then
        log_error "Failed to process alignment file for $model_name"
        unset alignment_params
        continue
    fi

    # Calculate glycan parameters using Method 1 - declare array first
    sample_pdb="${m%/}/$PDB_DIR_PATTERN/${pdb_prefix}_0/$PDB_FILE_PATTERN"
    declare -A glycan_params=()
    if ! calculate_glycan_parameters "$sample_pdb" "${alignment_params[chain_len]}" glycan_params; then
        log_error "Failed to calculate glycan parameters for $model_name"
        unset alignment_params glycan_params
        continue
    fi

    # Initialize output file
    out_file="${out_dir}/$OUTPUT_FILE"
    : > "$out_file"

    frame_idx=0
    atoms_expected=""

    # Process each run with NRUNS-based iteration
    for ((run_num = 0; run_num < NRUNS; run_num++)); do
        run_dir="${m%/}/$PDB_DIR_PATTERN/${pdb_prefix}_${run_num}/"
        
        if [[ ! -d "$run_dir" ]]; then
            log_warn "Missing directory $run_dir — skipping"
            continue
        fi
        
        run_id="$(basename "${run_dir%/}")"
        filename="${run_dir}$PDB_FILE_PATTERN"
        
        if [[ ! -f "$filename" ]]; then
            log_warn "Missing $filename — skipping"
            continue
        fi

        ((++frame_idx))

        # Process single frame
        atoms_count=$(process_single_frame \
            "$filename" \
            "$run_id" \
            "$frame_idx" \
            "$out_dir" \
            "$model_vmd_log_dir" \
            alignment_params \
            glycan_params)

        # Check if processing failed
        if [[ -z "$atoms_count" ]]; then
            log_warn "Failed to process frame $frame_idx"
            continue
        fi

        # Validate atom count consistency
        if [[ -z "${atoms_expected:-}" ]]; then
            atoms_expected="$atoms_count"
            log_info "First frame atoms=$atoms_expected"
        elif [[ "$atoms_count" -ne "$atoms_expected" ]]; then
            log_warn "Skipping frame $frame_idx: atoms $atoms_count != expected $atoms_expected"
            rm -f "${out_dir}/frame_${frame_idx}.pdb" || true
            continue
        fi

        # Append frame to output
        frame_payload="${out_dir}/frame_${frame_idx}.pdb"
        append_frame_to_output "$frame_idx" "$frame_payload" "$out_file"

        # Clean up frame payload
        rm -f "$frame_payload" || true
    done

    # Align the multi-model PDB if we have output
    if [[ -s "$out_file" ]]; then
        aligned_pdb="${out_dir}/$ALIGNED_OUTPUT_FILE"
        align_log="${model_vmd_log_dir}/align_trajectory_vmd.log"
        
        if ! align_trajectory "$out_file" "$aligned_pdb" "$align_log"; then
            log_error "Failed to align trajectory for $model_name"
            exit 1
        fi
        
        log_info "Successfully processed model $model_name with aligned output"
    else
        log_warn "No valid frames processed for model $model_name"
    fi

    # Clean up any remaining intermediate files
    cleanup_intermediate_files
    
    # Clean up arrays for this model
    unset alignment_params glycan_params
done

log_info "Ensemble modeling pipeline completed successfully"
