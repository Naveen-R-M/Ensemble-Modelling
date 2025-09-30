#!/usr/bin/env bash
# lib/glycan_calculator.sh - Glycan range calculation using Method 1 (Direct PDB Analysis)

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

# --- HETATM range detection ---
# Input: path to PDB file | Output: echoes "hetatm_start hetatm_end" or returns 1 if no HETATM found
get_hetatm_range() {
    local pdb_file="$1"
    
    log_info "Analyzing HETATM residues in: $(basename "$pdb_file")"
    
    # Get actual HETATM range from PDB file
    local hetatm_residues
    hetatm_residues=$(grep "^HETATM" "$pdb_file" | awk '{print substr($0,23,4)+0}' | sort -n | uniq)
    
    if [[ -z "$hetatm_residues" ]]; then
        log_error "No HETATM residues found in $pdb_file"
        return 1
    fi
    
    local hetatm_start hetatm_end
    hetatm_start=$(echo "$hetatm_residues" | head -1)
    hetatm_end=$(echo "$hetatm_residues" | tail -1)
    
    echo "$hetatm_start $hetatm_end"
}

# --- Equal glycan distribution calculation (always 3 glycan chains) ---
# Input: hetatm_start, hetatm_end (integers) | Output: echoes "g1_start g1_end g2_start g2_end g3_start g3_end"
calculate_equal_glycan_ranges() {
    local hetatm_start="$1"
    local hetatm_end="$2"
    
    local total_hetatm_residues=$((hetatm_end - hetatm_start + 1))
    log_info "HETATM residue range: $hetatm_start-$hetatm_end ($total_hetatm_residues total residues)"
    log_info "Distributing equally among 3 glycan chains (CAR1, CAR2, CAR3)"
    
    # Calculate equal distribution among 3 glycan chains
    local residues_per_chain=$((total_hetatm_residues / 3))
    local remainder=$((total_hetatm_residues % 3))
    
    # Distribute residues equally
    local g1_start="$hetatm_start"
    local g1_end=$((g1_start + residues_per_chain - 1))
    
    local g2_start=$((g1_end + 1))
    local g2_end=$((g2_start + residues_per_chain - 1))
    
    local g3_start=$((g2_end + 1))
    local g3_end="$hetatm_end"  # Use actual end to handle any remainder
    
    # Adjust for remainder by redistributing if needed
    if [[ $remainder -gt 0 ]]; then
        log_info "Distributing remainder of $remainder residues"
        if [[ $remainder -ge 1 ]]; then
            ((g1_end++))
            ((g2_start++))
            ((g2_end++))
            ((g3_start++))
        fi
        if [[ $remainder -ge 2 ]]; then
            ((g2_end++))
            ((g3_start++))
        fi
    fi
    
    # Calculate and log final distribution
    local g1_residues=$((g1_end - g1_start + 1))
    local g2_residues=$((g2_end - g2_start + 1))
    local g3_residues=$((g3_end - g3_start + 1))
    
    log_info "Equal glycan distribution:"
    log_info "  CAR1: $g1_start-$g1_end ($g1_residues residues)"
    log_info "  CAR2: $g2_start-$g2_end ($g2_residues residues)"
    log_info "  CAR3: $g3_start-$g3_end ($g3_residues residues)"
    
    # Return all ranges
    echo "$g1_start $g1_end $g2_start $g2_end $g3_start $g3_end"
}

# --- Main glycan calculation function ---
# Input: path to sample PDB, chain_len, associative array reference | Output: populates array with glycan range parameters (g1_start, g1_end, g2_start, g2_end, g3_start, g3_end, hetatm_start, hetatm_end)
calculate_glycan_parameters() {
    local sample_pdb="$1"
    local chain_len="$2"          # Chain length from alignment (for logging only)
    local -n glycan_result_ref="$3"  # Pass associative array by reference
    
    # Validate input file
    if [[ ! -f "$sample_pdb" ]]; then
        log_error "Sample PDB not found: $sample_pdb"
        return 1
    fi
    
    log_info "Protein has $chain_len chains, using 3 glycan chains (CAR1, CAR2, CAR3)"
    
    # Get HETATM range
    local hetatm_range
    hetatm_range=$(get_hetatm_range "$sample_pdb")
    if [[ -z "$hetatm_range" ]]; then
        return 1
    fi
    
    read -r hetatm_start hetatm_end <<< "$hetatm_range"
    
    # Calculate glycan ranges using Method 1 (Equal Distribution) - always 3 chains
    local glycan_ranges
    glycan_ranges=$(calculate_equal_glycan_ranges "$hetatm_start" "$hetatm_end")
    read -r g1_start g1_end g2_start g2_end g3_start g3_end <<< "$glycan_ranges"
    
    # Validate all parameters
    for param in "$g1_start" "$g1_end" "$g2_start" "$g2_end" "$g3_start" "$g3_end"; do
        validate_integer "$param" "glycan_parameter"
    done
    
    # Populate result array (always 3 glycan chains)
    glycan_result_ref[g1_start]="$g1_start"
    glycan_result_ref[g1_end]="$g1_end"
    glycan_result_ref[g2_start]="$g2_start"
    glycan_result_ref[g2_end]="$g2_end"
    glycan_result_ref[g3_start]="$g3_start"
    glycan_result_ref[g3_end]="$g3_end"
    glycan_result_ref[hetatm_start]="$hetatm_start"
    glycan_result_ref[hetatm_end]="$hetatm_end"
    
    return 0
}

# --- Alternative calculation methods (for future use) ---
# Input: PDB file path, glyc.dat path, align.ali path | Output: placeholder (not implemented)
calculate_attachment_based_ranges() {
    local pdb_file="$1"
    local glyc_file="$2"
    local ali_file="$3"
    
    # Placeholder for Method 2 implementation
    log_warn "Attachment-based glycan calculation not yet implemented"
    return 1
}

# Input: PDB file path | Output: placeholder (not implemented)
calculate_residue_type_based_ranges() {
    local pdb_file="$1"
    
    # Placeholder for Method 3 implementation  
    log_warn "Residue-type based glycan calculation not yet implemented"
    return 1
}
