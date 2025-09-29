#!/usr/bin/env bash
# lib/alignment_parser.sh - Alignment file analysis functions

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

# --- Alignment analysis (replaces compute_glycan_params.py) ---
analyze_alignment() {
    local ali_file="$1"
    
    log_info "Analyzing alignment file: $ali_file"
    
    # Initialize variables
    local aa_count=0
    local slash_count=0
    local before_first_slash_count=0
    local in_pm=0
    local seen_first_slash=0
    
    while IFS= read -r line; do
        line="${line%$'\r'}"  # Remove carriage return
        line="${line# }"      # Remove leading spaces
        
        if [[ "$line" == ">P1;pm.pdb"* ]]; then
            in_pm=1
            continue
        elif [[ "$line" == ">P1;"* && "$line" != ">P1;pm.pdb"* ]]; then
            in_pm=0
            continue
        elif [[ $in_pm -eq 1 ]]; then
            if [[ "$line" == "structureX"* || "$line" == "sequence:"* || "$line" == ">"* ]]; then
                continue
            fi
            
            # Process each character in the line
            for ((i=0; i<${#line}; i++)); do
                char="${line:$i:1}"
                if [[ "$char" == "/" ]]; then
                    if [[ $seen_first_slash -eq 0 ]]; then
                        seen_first_slash=1
                    fi
                    ((slash_count++))
                elif [[ "$char" != "*" ]]; then
                    if [[ $seen_first_slash -eq 0 ]]; then
                        ((before_first_slash_count++))
                    fi
                    ((aa_count++))
                fi
            done
        fi
    done < "$ali_file"
    
    # Return values: aa_count slash_count before_first_slash_count
    echo "$aa_count $slash_count $before_first_slash_count"
}

# --- Determine chain length from alignment ---
determine_chain_length() {
    local slash_count="$1"
    
    if [[ "$slash_count" -eq 2 ]]; then
        echo "3"
    elif [[ "$slash_count" -eq 5 ]]; then
        echo "6"
    else
        log_error "Unexpected number of slashes in alignment: $slash_count (expected 2 or 5)"
        return 1
    fi
}

# --- Main alignment processing function ---
process_alignment() {
    local ali_file="$1"
    local -n result_ref="$2"  # Pass associative array by reference
    
    # Analyze alignment
    local analysis_result
    analysis_result=$(analyze_alignment "$ali_file")
    read -r aa_count slash_count first_chain_length <<< "$analysis_result"
    
    # Determine chain length
    local chain_len
    chain_len=$(determine_chain_length "$slash_count")
    
    # Populate result array
    result_ref[aa_count]="$aa_count"
    result_ref[slash_count]="$slash_count"
    result_ref[first_chain_length]="$first_chain_length"
    result_ref[chain_len]="$chain_len"
    result_ref[first_chain_len_plus_one]="$((first_chain_length + 1))"
    
    log_info "Detected $chain_len-chain protein system"
    log_info "Alignment analysis: aa_count=$aa_count, chains=$chain_len, first_chain_len=$first_chain_length"
    
    return 0
}
