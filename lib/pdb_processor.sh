#!/usr/bin/env bash
# lib/pdb_processor.sh - PDB processing and VMD operations

# Ensure we can find the common.sh file
if [[ -f "${BASH_SOURCE[0]%/*}/common.sh" ]]; then
    # shellcheck source=lib/common.sh
    source "${BASH_SOURCE[0]%/*}/common.sh"
elif [[ -f "$(dirname "${BASH_SOURCE[0]}")/common.sh" ]]; then
    # shellcheck source=lib/common.sh
    source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
else
    echo "ERROR: Cannot find common.sh" >&2
    exit 1
fi

# --- VMD chain extraction ---
extract_chains_with_vmd() {
    local filename="$1"
    local g1_start="$2" g1_end="$3"
    local g2_start="$4" g2_end="$5"
    local g3_start="$6" g3_end="$7"
    local run_id="$8"
    local log_dir="$9"
    
    log_info "Extracting chains for run $run_id"
    
    if ! vmd -dispdev text -eofexit \
          -e "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" -args \
          "$filename" "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end" \
          > "${log_dir}/vmd_run_${run_id}.log" 2>&1; then
        log_error "VMD extract failed for run $run_id (see log)"
        return 1
    fi
    
    return 0
}

# --- Determine starting chain number ---
get_starting_chain_number() {
    local filename="$1"
    local chain="${2:-C}"  # Default to chain C
    
    local starting_chain_number
    starting_chain_number=$(
        pdb_tidy "$filename" \
            | pdb_selchain -"$chain" \
            | awk '/^(ATOM|HETATM)/{print substr($0,23,4)+0; exit}'
    )
    
    echo "${starting_chain_number:-}"
}

# --- Renumber PDB components ---
renumber_pdb_components() {
    local chain_len="$1"
    local starting_chain_number="$2"
    local first_chain_len_plus_one="$3"
    
    log_info "Renumbering PDB components for $chain_len-chain system"
    
    if [[ "$chain_len" -eq 3 ]]; then
        pdb_reres -"$starting_chain_number" CHA.pdb > CHA_renum.pdb
        pdb_reres -"$starting_chain_number" CHB.pdb > CHB_renum.pdb
        pdb_reres -"$starting_chain_number" CHC.pdb > CHC_renum.pdb
        pdb_reres -1 CAR1.pdb > CAR1_renum.pdb
        pdb_reres -1 CAR2.pdb > CAR2_renum.pdb
        pdb_reres -1 CAR3.pdb > CAR3_renum.pdb
        
        echo "CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb"
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
        
        echo "CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb CHD_renum.pdb CHE_renum.pdb CHF_renum.pdb CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb"
    fi
}

# --- Validate required components ---
validate_components() {
    local -a concat_files=("$@")
    
    for f in "${concat_files[@]}"; do
        if [[ ! -s "$f" ]]; then
            log_warn "Missing or empty component: $f"
            return 1
        fi
    done
    
    return 0
}

# --- Concatenate PDB files ---
concatenate_pdb_files() {
    local temp_file="$1"
    shift
    local -a concat_files=("$@")
    
    log_info "Concatenating PDB files in fixed order"
    
    {
        for file in "${concat_files[@]}"; do
            awk '{
                  sub(/\r$/, "")
                  if ($0 ~ /^(CRYST1|REMARK|TER|END|ENDMDL|MODEL)([[:space:]].*)?$/) next
                  print
                }' "$file"
        done
    } > "$temp_file"
}

# --- Reindex and count atoms ---
process_frame_payload() {
    local temp_file="$1"
    local frame_payload="$2"
    
    pdb_reatom -1 < "$temp_file" > "$frame_payload"
    
    local atoms_count
    atoms_count=$(awk '/^(ATOM|HETATM)/{n++} END{print n+0}' "$frame_payload")
    echo "$atoms_count"
}

# --- Append frame to output (MODIFIED: uses END instead of MODEL/ENDMDL) ---
append_frame_to_output() {
    local frame_idx="$1"
    local frame_payload="$2"
    local out_file="$3"
    
    {
        cat "$frame_payload"
        echo "END"
    } >> "$out_file"
}

# --- Main frame processing function ---
process_single_frame() {
    local filename="$1"
    local run_id="$2"
    local frame_idx="$3"
    local out_dir="$4"
    local log_dir="$5"
    local -n alignment_params="$6"
    local -n glycan_params="$7"
    
    log_info "Processing frame $frame_idx (run $run_id)"
    
    local temp_file="${out_dir}/output_temp.pdb"
    local frame_payload="${out_dir}/frame_${frame_idx}.pdb"
    
    # Initialize temp file
    : > "$temp_file"
    
    # Work inside out_dir (so VMD drops CHA/CHB/... here)
    pushd "$out_dir" >/dev/null
    
    # Extract chains & glycans with VMD
    if ! extract_chains_with_vmd \
        "$filename" \
        "${glycan_params[g1_start]}" "${glycan_params[g1_end]}" \
        "${glycan_params[g2_start]}" "${glycan_params[g2_end]}" \
        "${glycan_params[g3_start]}" "${glycan_params[g3_end]}" \
        "$run_id" "$log_dir"; then
        cleanup_intermediate_files
        popd >/dev/null
        return 1
    fi
    
    # Determine starting residue index
    local starting_chain_number
    starting_chain_number=$(get_starting_chain_number "$filename" "C")
    if [[ -z "$starting_chain_number" ]]; then
        log_warn "Could not determine starting chain number for $run_id"
        cleanup_intermediate_files
        popd >/dev/null
        return 1
    fi
    
    # Renumber pieces
    local concat_files_str
    concat_files_str=$(renumber_pdb_components \
        "${alignment_params[chain_len]}" \
        "$starting_chain_number" \
        "${alignment_params[first_chain_len_plus_one]}")
    
    # Convert string to array
    read -ra concat_files <<< "$concat_files_str"
    
    # Validate all required components exist
    if ! validate_components "${concat_files[@]}"; then
        log_warn "Skipping frame $frame_idx: missing component(s)"
        cleanup_intermediate_files
        popd >/dev/null
        return 1
    fi
    
    # Concatenate files
    concatenate_pdb_files "$temp_file" "${concat_files[@]}"
    
    # Process frame payload and count atoms
    local atoms_count
    atoms_count=$(process_frame_payload "$temp_file" "$frame_payload")
    
    # Cleanup intermediates
    rm -f "$temp_file" C*.pdb *_renum.pdb || true
    popd >/dev/null
    
    # Return atom count for validation
    echo "$atoms_count"
}

# --- Trajectory alignment ---
align_trajectory() {
    local in_pdb="$1"
    local out_pdb="$2"
    local align_log="$3"
    
    log_info "Aligning trajectory: $(basename "$in_pdb")"
    
    if ! vmd -dispdev text -eofexit \
          -e "$ALIGN_TRAJECTORY_SCRIPT" \
          -args "$in_pdb" "$out_pdb" \
          > "$align_log" 2>&1; then
        log_error "VMD alignment failed (see log: $align_log)"
        return 1
    fi
    
    return 0
}
