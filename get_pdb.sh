#!/usr/bin/env bash
#SBATCH --job-name=get_pdb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=short

# ==============================================================================
# get_pdb.sh - Ensemble Modeling Pipeline for Glycoprotein Structures
# ==============================================================================
#
# PURPOSE:
#   Processes AllosMod ensemble model outputs to create multi-frame PDB files
#   with properly separated protein chains and glycan chains, then aligns them.
#
# GLYCAN PARAMETER DEFINITIONS:
#   g1_start, g1_end: Residue range defining CAR1 (1st glycan chain attachment)
#                     VMD uses this to select HETATM residues 3394-3549 for CAR1
#   g2_start, g2_end: Residue range defining CAR2 (2nd glycan chain attachment)
#                     VMD uses this to select HETATM residues 3550-3705 for CAR2
#   g3_start, g3_end: Residue range defining CAR3 (3rd glycan chain attachment)
#                     VMD uses this to select HETATM residues 3706-3861 for CAR3
#   hetatm_start, hetatm_end: Overall HETATM boundaries (e.g., 3394-3861)
#                             Stored for reference/validation, not used in extraction
#
# WORKFLOW:
#   1. Analyze alignment file to determine protein chain structure (3 or 6 chains)
#   2. Detect HETATM range from sample PDB and calculate equal glycan distribution
#   3. For each frame (start.pdb_0 to start.pdb_N):
#      - Extract protein chains (CHA, CHB, CHC[, CHD, CHE, CHF]) using VMD
#      - Extract glycan chains (CAR1, CAR2, CAR3) using calculated residue ranges
#      - Renumber residues for consistency
#      - Concatenate into single frame and append to output.pdb
#   4. Align all frames using protein backbone and write to output_aligned.pdb
#
# OUTPUT FILES:
#   output.pdb: Multi-frame PDB with all extracted frames, each ending with END
#   output_aligned.pdb: Backbone-aligned version of output.pdb
#
# ==============================================================================

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

# Input: None | Output: Sets up VMD, conda environment, exports VMDNOSTARTUP
setup_modules

# Input: folder path | Output: Sources .env files from various locations
load_env_files "$folder"

# Input: USER_ID | Output: Creates and exports LOG_DIR, OUT_BASE, SCRIPTS_BASE directory paths
setup_directories "$USER_ID"

# Input: None | Output: Redirects stdout/stderr to log files
setup_logging

# Input: None | Output: Validates required tools exist (awk, vmd, pdb_*, etc.)
check_required_tools

# Input: None | Output: Validates and exports ALIGN_TRAJECTORY_SCRIPT, EXTRACT_GLYCANS_AND_CHAINS_SCRIPT paths
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

    # Input: path to input.dat file | Output: NRUNS value (integer)
    NRUNS=$(parse_input_dat "$m/$INPUT_FILE")
    if [[ -z "$NRUNS" ]]; then
        log_warn "Skipping $model_name (no valid $INPUT_FILE or NRUNS)"
        continue
    fi
    log_info "Found NRUNS=$NRUNS for model $model_name"

    # Input: model directory path 
    # Output: PDB prefix string (e.g., "start.pdb", "bg505.pdb")
    pdb_prefix=$(detect_pdb_prefix "$m")
    if [[ -z "$pdb_prefix" ]]; then
        log_error "Could not detect PDB prefix pattern in $m/$PDB_DIR_PATTERN/"
        continue
    fi
    log_info "Detected PDB prefix pattern: ${pdb_prefix}_*"

    # Input: align.ali file path, reference to alignment_params array 
    # Output: Populates alignment_params[aa_count, slash_count, first_chain_length, chain_len, first_chain_len_plus_one]
    declare -A alignment_params=()
    if ! process_alignment "$m/$ALIGNMENT_FILE" alignment_params; then
        log_error "Failed to process alignment file for $model_name"
        unset alignment_params
        continue
    fi

    # Input: sample PDB file path, chain_len, reference to glycan_params array
    # Output: Populates glycan_params array with:
    #   g1_start, g1_end: Residue range for CAR1 (first glycan chain) - used by VMD to select and extract glycan atoms
    #   g2_start, g2_end: Residue range for CAR2 (second glycan chain) - used by VMD to select and extract glycan atoms
    #   g3_start, g3_end: Residue range for CAR3 (third glycan chain) - used by VMD to select and extract glycan atoms
    #   hetatm_start, hetatm_end: Overall HETATM residue boundaries in PDB - stored for reference/debugging (not actively used in processing)
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

        # Input: PDB filename, run_id, frame_idx, output dir, log dir, alignment_params array, glycan_params array 
        # Output: Atom count (integer) or empty on failure
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

        # Input: frame_idx, frame_payload file path, output file path 
        # Output: Appends frame to output.pdb with END terminator
        frame_payload="${out_dir}/frame_${frame_idx}.pdb"
        append_frame_to_output "$frame_idx" "$frame_payload" "$out_file"

        # Clean up frame payload
        rm -f "$frame_payload" || true
    done

    # Align the multi-model PDB if we have output
    if [[ -s "$out_file" ]]; then
        aligned_pdb="${out_dir}/$ALIGNED_OUTPUT_FILE"
        align_log="${model_vmd_log_dir}/align_trajectory_vmd.log"
        
        # Input: input PDB path, output aligned PDB path, log file path 
        # Output: Creates aligned PDB file, returns 0 on success
        if ! align_trajectory "$out_file" "$aligned_pdb" "$align_log"; then
            log_error "Failed to align trajectory for $model_name"
            exit 1
        fi
        
        log_info "Successfully processed model $model_name with aligned output"
    else
        log_warn "No valid frames processed for model $model_name"
    fi

    # Input: None | Output: Removes C*.pdb and *_renum.pdb files
    cleanup_intermediate_files
    
    # Clean up arrays for this model
    unset alignment_params glycan_params
done

log_info "Ensemble modeling pipeline completed successfully"