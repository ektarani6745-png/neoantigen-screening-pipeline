# Neoantigen Screening Pipeline

A computational oncology workflow designed to identify immunogenic neoepitopes. The pipeline processes somatic variant data to extract mutant 9-mer peptides and evaluates their Kyte-Doolittle hydrophobicity to prioritize candidate neoantigens for cancer immunotherapy.

## Features
* **9-mer Peptide Extraction:** Automatically isolates 9-mer windows centered on amino acid substitutions from full protein sequences.
* **Hydrophobicity Profiling:** Calculates the Kyte-Doolittle score for both Wild-Type (WT) and Mutant (MT) peptides.
* **Shift Analysis:** Quantifies the biochemical shift in hydrophobicity (`Hydrophobicity_Shift`), a metric used to evaluate potential T-cell receptor (TCR) recognition.
* **Automated Bash Wrapper:** Includes a streamlined execution script (`run_pipeline.sh`) that processes data and outputs a formatted terminal table.

## Requirements
* **Python 3.x**
* **Libraries:** `pandas`
* **OS:** Linux / macOS

## Usage
**Option 1: Run the automated pipeline**
```bash
./run_pipeline.sh
