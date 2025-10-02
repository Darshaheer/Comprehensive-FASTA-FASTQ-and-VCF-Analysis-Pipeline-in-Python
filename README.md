# Comprehensive-FASTA-FASTQ-and-VCF-Analysis-Pipeline-in-Python
This Python pipeline analyzes genomic data from FASTA, FASTQ, and VCF files. It extracts sequences, headers, GC content, motifs, Phred scores, nucleotide frequencies, palindromes, and variant details, categorizing SNPs, insertions, deletions, and multi-allelic variants. Results and summary statistics are saved in text and JSON formats.

## Features
The pipeline integrates the following functionalities:

1. **FASTA File Analysis**
   - Extracts the header and full nucleotide sequence.
   - Calculates GC content for the sequence.
   - Identifies motif positions, such as the TATA box.

2. **FASTQ File Analysis**
   - Parses sequence reads and associated quality scores.
   - Calculates per-read Phred quality scores.
   - Computes GC content for each read.
   - Counts nucleotide frequencies (A, T, G, C).
   - Detects motif occurrences (e.g., TATA box, AAA, A{5}).
   - Identifies palindromic sequences.
   - Determines high-quality reads and maximum-quality reads.

3. **VCF File Analysis**
   - Extracts variant information including chromosome, position, reference and alternate alleles, quality scores, and filter status.
   - Computes allele frequencies and genotype information.
   - Categorizes variants into SNPs, insertions, deletions, and multi-allelic variants.
   - Identifies common, rare, and high-quality variants.
   - Extracts variants within specific genomic ranges (e.g., chr22: 15,000,000–20,000,000).

4. **Motif and Sequence Statistics**
   - Reports positions of motifs in the first set of reads.
   - Counts motifs occurrences across reads.
   - Identifies sequences ending with specific patterns (e.g., AAA).
   - Calculates GC content per read and overall averages.
   - Detects sequences of maximal length and highest GC content.
   - Reports palindromic sequences in the dataset.

5. **Summary and Output**
   - Generates detailed summary statistics across all processed files.
   - Writes comprehensive results to a text file (`Parsed Genomic File.txt`).
   - Saves parsed variant information as a JSON dictionary (`Variant Dictionary.json`).

## Input Requirements
- **FASTA File:** Contains nucleotide sequences for analysis.
- **FASTQ File:** Raw sequencing reads with Phred quality scores.
- **VCF File:** Variant Call Format file for variant extraction and analysis.

## Output Files
- `Parsed Genomic File.txt` – Text file summarizing FASTA, FASTQ, and VCF analyses.
- `Variant Dictionary.json` – JSON file containing parsed variant details.

## Usage
1. Place your input files (`FASTA`, `FASTQ`, and `VCF`) in the same directory as the pipeline script.
2. Update the filenames in the script if necessary.
3. Run the script using Python 3:
   ```bash
   python Comprehensive Parsing.py
Check the generated output files for detailed results and statistics.

## Dependencies
  Python 3.13
  Standard libraries: re, json, collections

## Applications
  Basic genomics and sequencing data analysis
  Motif detection and sequence statistics
  Variant parsing and classification
  Quality assessment of sequencing reads
  Generation of structured summary reports for downstream analysis

## Notes
  The pipeline currently processes the first 500 reads for certain sequence statistics; this can be modified according to dataset size.
  Allele frequency thresholds for common/rare variants and quality thresholds are adjustable in the script.


## License
MIT License
