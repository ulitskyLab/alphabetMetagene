# alphabetMetagene
Analysis of input sequence data.

This repository contains generic standalone tools for sequence analysis.

## What data is it for?

The alphabet metagene code allows the analysis of sequences of various alphabets (DNA/RNA/Protein/other alphabet) and outputs the dynamic proportions of the letter across the binned sequence.

#### Input:

* Fasta File(s) of Sequence(s) of same/varying lengths.
* Output path

##### Additional optional input includs:

* File(s) that matches between individual sequences from the first (sample) input fasta file and the other (control) sequence file(s) [Example provided]. If a single fasta file is used for the analysis- ignore this option. If multiple fasta files are used, this would allow a calculation of a paired Wilcoxon rank-sum test instead of a non-paired one (that would be used by default).
* Alphabet to consider (DNA is the default).
* Sample names
* Title
* Plot colors
* Additional options that can be seen using the `-h` option

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
    
Running example:

    Rscript overallAlphabetContent.R -f fasta_file1_path.fa,fasta_file2_control_path.fa,fasta_fileN_control_path.fa -m matching_first_file_to_control2.txt,matching_first_file_to_controlN.txt -t 'title' -s name_file1,name_file2,name_fileN -r png -o output_path/file_name.png 
    
Output example:

![overallAlphabetContentOutputExample](https://user-images.githubusercontent.com/87706940/135989900-1b70ec01-af08-4e92-9001-27896ebc5cab.png)


### alphabetMetagene.R

View documentation:

    Rscript alphabetMetagene.R -h

Running example:

    Rscript alphabetMetagene.R -f fasta_file1_path.fa,fasta_file2_control_path.fa,fasta_fileN_control_path.fa -m matching_first_file_to_control2.txt,matching_first_file_to_controlN.txt -s name_file1,name_file2,name_fileN -o output_path/file_name.pdf

Output example:

![alphabetMetageneOutputExample](https://user-images.githubusercontent.com/87706940/136018688-ee65dff3-cabe-40ff-aaed-72be686088cb.png)

