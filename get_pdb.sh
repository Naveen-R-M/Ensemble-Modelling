#!/usr/bin/env bash
set -euo pipefail

# --- usage ---
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <folder>"
  exit 1
fi
folder=$1

# --- sanity checks ---
command -v awk >/dev/null || { echo "awk not found"; exit 1; }
command -v vmd  >/dev/null || { echo "vmd not found"; exit 1; }
command -v pdb_tidy  >/dev/null || { echo "pdb_tidy (pdb-tools) not found"; exit 1; }
command -v pdb_reres >/dev/null || { echo "pdb_reres (pdb-tools) not found"; exit 1; }
command -v pdb_reatom >/dev/null || { echo "pdb_reatom (pdb-tools) not found"; exit 1; }

# helper scripts (adjust names/paths as needed)
ALIGN_TRAJECTORY_SCRIPT='align_trajectory.tcl'
COMPUTE_GLYCAN_PARAMS_SCRIPT='compute_glycan_params.py'
EXTRACT_GLYCANS_AND_CHAINS_SCRIPT='extract_glycans_and_chains.tcl'

[[ -f "$ALIGN_TRAJECTORY_SCRIPT" ]] || { echo "Missing $ALIGN_TRAJECTORY_SCRIPT"; exit 1; }
[[ -f "$COMPUTE_GLYCAN_PARAMS_SCRIPT" ]] || { echo "Missing $COMPUTE_GLYCAN_PARAMS_SCRIPT"; exit 1; }
[[ -f "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" ]] || { echo "Missing $EXTRACT_GLYCANS_AND_CHAINS_SCRIPT"; exit 1; }

# Read NRUNS from input.dat (portable, digits-only)
NRUNS=$(awk -F= '/^NRUNS=/{match($2,/[0-9]+/,m); print m[0]; exit}' "$folder/input.dat")
[[ -n "${NRUNS:-}" ]] || { echo "NRUNS not found in $folder/input.dat"; exit 1; }

for m in "$folder"/*; do
  [[ -d "$m" ]] || continue

  # Compute glycan parameters once per subfolder
  read -r g1_start g1_end g2_start g2_end g3_start g3_end chain_len first_chain_len_plus_one \
    < <(python3 "$COMPUTE_GLYCAN_PARAMS_SCRIPT" "$m/align.ali" "$m/glyc.dat")

  # One output file per subfolder (kept inside the subfolder)
  out_file="$m/output_$(basename "$m").pdb"
  : > "$out_file"

  for i in $(seq 0 $((NRUNS - 1))); do
    filename="${folder}/pred_dECALCrAS1000/siv.pdb_${i}/pm.pdb.B99990001.pdb"
    echo "Processing $filename"

    # Run VMD with computed glycan parameters
    vmd -dispdev text -e "$EXTRACT_GLYCANS_AND_CHAINS_SCRIPT" -args "$filename" \
      "{$g1_start $g1_end}" "{$g2_start $g2_end}" "{$g3_start $g3_end}"

    # Extract first residue number of chain A from the source PDB
    starting_chain_number=$(
      pdb_tidy "$filename" \
        | pdb_selchain -A \
        | awk '/^(ATOM|HETATM)/{print substr($0,23,4)+0; exit}'
    )
    if [[ -z "${starting_chain_number:-}" ]]; then
      echo "Warning: no ATOM/HETATM records for chain A in $filename — skipping run $i."
      continue
    fi

    # Renumber residue IDs in chain files as required
    if [[ "$chain_len" -eq 3 ]]; then
      pdb_reres -"$starting_chain_number" CHA.pdb > CHA_renum.pdb
      pdb_reres -"$starting_chain_number" CHB.pdb > CHB_renum.pdb
      pdb_reres -"$starting_chain_number" CHC.pdb > CHC_renum.pdb
      pdb_reres -1 CAR1.pdb > CAR1_renum.pdb
      pdb_reres -1 CAR2.pdb > CAR2_renum.pdb
      pdb_reres -1 CAR3.pdb > CAR3_renum.pdb

      concat_files=(CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb
                    CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb)
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

      concat_files=(CHA_renum.pdb CHB_renum.pdb CHC_renum.pdb
                    CHD_renum.pdb CHE_renum.pdb CHF_renum.pdb
                    CAR1_renum.pdb CAR2_renum.pdb CAR3_renum.pdb)
    fi

    # Filter out any missing/empty files before piping to awk
    files=()
    for f in "${concat_files[@]}"; do
      [[ -s "$f" ]] && files+=("$f")
    done
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "Warning: no input pieces produced for $(basename "$m"), run $i — skipping."
      continue
    fi

    # Single pipeline: filter -> reatom -> append END, then append to the run's output
    {
      awk '!/^(CRYST1|REMARK|END)$/' -- "${files[@]}" \
      | pdb_reatom -1
      echo 'END'
    } >> "$out_file"

    # Clean per-run intermediates produced by the Tcl step/renumbering
    rm -f C*.pdb
  done

  echo "Aligning trajectory for $(basename "$m")..."
  vmd -dispdev text -e "$ALIGN_TRAJECTORY_SCRIPT" -args "$out_file" "$m/output_aligned.pdb" "protein and backbone"
done

# Final cleanup (defensive)
rm -f C*.pdb || true
