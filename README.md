# alphabetMetagene
Analysis of input sequence data.

This repository contains generic standalone tools for sequence analysis.

## Installation Requirements:
* R >= 3.3.2
* gcc >= 9.2.0
* optparse
* plyr
* stringr
* reshape2
* pheatmap
* gridExtra

## Usage:

### overallAlphabetContent.R

View documentation:

    Rscript overallAlphabetContent.R -h
    
Run:

    Rscript overallAlphabetContent.R -f fasta_file1_path.fa,fasta_file2_control_path.fa,fasta_fileN_control_path.fa -m matching_first_file_to_control2.txt,matching_first_file_to_controlN.txt -t 'title' -s name_file1,name_file2,name_file3 -o output_path/file_name.pdf
    

![HepG2_ACGT_tandem_nuc_matched_controls.pdf](https://github.com/ulitskyLab/alphabetMetagene/files/7283977/HepG2_ACGT_tandem_nuc_matched_controls.pdf)

![overallAlphabetContentExample](examples/overallAlphabetContentOutputExample.pdf)
