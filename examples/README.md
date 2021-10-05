## Running code using the example files:

### overall alphabet content:

    Rscript overallAlphabetContent.R -f sample.fa,control1.fa,control2.fa -m type1_matching_file_control1.txt,type1_matching_file_control2.txt -o overall_out.pdf

### alphabet metagene across the sequence

    Rscript alphabetMetagene.R -f sample.fa,control1.fa,control2.fa -m type1_matching_file_control1.txt,type1_matching_file_control2.txt -o metagene_out.pdf

### Note
* The provided matching control files are of 2 columns. Alternatively, they could be in a manipulated (extended) BED6 format, with the 4th column containing the name of the control sequence and the 7th (additional) column being the name of the matched sample sequence.
