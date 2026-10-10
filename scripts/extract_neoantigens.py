import sys
import pandas as pd

# Kyte-Doolittle Hydrophobicity Scale
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def calculate_score(peptide):
    """Calculates average Kyte-Doolittle hydrophobicity for a given peptide."""
    score = sum(HYDROPHOBICITY.get(aa, 0.0) for aa in peptide)
    return round(score / len(peptide), 3)

def generate_neoantigens(input_tsv, output_tsv, peptide_length=9):
    df = pd.read_csv(input_tsv, sep='\t')
    flank = peptide_length // 2  # 4 residues on each side for a 9-mer
    results = []

    for _, row in df.iterrows():
        pos = int(row['Protein_Position']) - 1  # 0-indexed position
        full_seq = row['Full_Protein_Seq']
        alt_aa = row['Alt_AA']
        ref_aa = row['Ref_AA']

        start = max(0, pos - flank)
        end = min(len(full_seq), pos + flank + 1)

        # Validate that the peptide length is exactly the desired 9-mer
        if (end - start) == peptide_length:
            # Construct Wild-Type (WT) and Mutant (MT) peptides
            wt_peptide = full_seq[start:end]
            mt_peptide = full_seq[start:pos] + alt_aa + full_seq[pos + 1:end]

            wt_hydro = calculate_score(wt_peptide)
            mt_hydro = calculate_score(mt_peptide)

            results.append({
                'Gene': row['Hugo_Symbol'],
                'Mutation': row['Protein_Change'],
                'WT_9mer': wt_peptide,
                'MT_Neoantigen_9mer': mt_peptide,
                'WT_Hydrophobicity': wt_hydro,
                'MT_Hydrophobicity': mt_hydro,
                'Hydrophobicity_Shift': round(mt_hydro - wt_hydro, 3)
            })

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"[✓] Extracted {len(out_df)} candidate neoantigens saved to {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_neoantigens.py <input_mutations.tsv> <output_candidates.tsv>")
        sys.exit(1)
    generate_neoantigens(sys.argv[1], sys.argv[2])
