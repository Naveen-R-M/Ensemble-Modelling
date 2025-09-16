#!/usr/bin/env python3
"""
compute_glycan_params.py

Reads:
  1. align.ali (alignment file)
  2. glyc.dat  (glycan specification file)

Outputs:
  g1_start g1_end g2_start g2_end g3_start g3_end chain_len first_chain_len

These values are read in getpdb bash script as:
  read g1_start g1_end g2_start g2_end g3_start g3_end chain_len first_chain_len < <(python3 compute_glycan_params.py align.ali glyc.dat)
"""

import sys

def count_amino_acids_and_chains(filepath):
    """
    Parse .ali file and return:
      - total amino acid count (aa_count)
      - number of '/' (slash_count)
      - residues before first '/' (first_chain_length)
    """
    aa_count = 0
    slash_count = 0
    before_first_slash_count = 0
    in_pm = False
    seen_first_slash = False

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">P1;pm.pdb"):
                in_pm = True
                continue
            elif line.startswith(">P1;") and not line.startswith(">P1;pm.pdb"):
                in_pm = False
                continue
            elif in_pm:
                if line.startswith("structureX") or line.startswith("sequence:") or line.startswith(">"):
                    continue

                if not seen_first_slash:
                    for c in line:
                        if c == "/":
                            seen_first_slash = True
                            break
                        elif c != "*":
                            before_first_slash_count += 1

                slash_count += line.count("/")
                aa_count += sum(1 for c in line if c not in {"/", "*"})

    return aa_count, slash_count, before_first_slash_count


def count_lines(filepath):
    """Count number of lines in glyc.dat"""
    with open(filepath, "r") as file:
        return sum(1 for _ in file)


def compute_all(ali_file, glyc_file):
    """
    Compute glycan chain positions based on .ali and glyc.dat.
    """
    pm_count, pm_slash, first_chain_length = count_amino_acids_and_chains(ali_file)
    glycan_lines = count_lines(glyc_file)

    # if glycan_lines % 3 != 0:
    #     raise ValueError("glyc.dat line count must be divisible by 3")

    # Chain length from number of slashes
    if pm_slash == 2:
        chain_length = 3
    elif pm_slash == 5:
        chain_length = 6
    else:
        raise ValueError(f"Unexpected number of slashes in alignment: {pm_slash}")

    # Each glycan section has equal number of lines
    glycan_per_chain = glycan_lines // 3

    glycan_chain1_start = pm_count + 1
    glycan_chain1_end = glycan_chain1_start + glycan_per_chain - 1

    glycan_chain2_start = glycan_chain1_end + 1
    glycan_chain2_end = glycan_chain2_start + glycan_per_chain - 1

    glycan_chain3_start = glycan_chain2_end + 1
    glycan_chain3_end = glycan_chain3_start + glycan_per_chain - 1

    return (
        glycan_chain1_start,
        glycan_chain1_end,
        glycan_chain2_start,
        glycan_chain2_end,
        glycan_chain3_start,
        glycan_chain3_end,
        chain_length,
        first_chain_length + 1  # +1 to match bash script usage
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: compute_glycan_params.py <align.ali> <glyc.dat>")
        sys.exit(1)

    ali_file = sys.argv[1]
    glyc_file = sys.argv[2]

    params = compute_all(ali_file, glyc_file)
    print(" ".join(map(str, params)))