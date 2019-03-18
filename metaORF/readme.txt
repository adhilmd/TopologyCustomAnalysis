
MetaORF was used to study the averaged enriched peak profile across all protein coding genes or specific set of genes with upstream  and downstream in the genome.

Required R packages:
plyr
dplyr
tidyr
zoo
argparse
ggplot2
gridExtra

Input file:
bed graph files

Steps:
i) Run the bdgprocess_scores_weight.sh with the directory containing the bedgraph files and required parameters:
bdg input file paths comma seperated (Mandatory) --bdgpaths=
sample name comma seperated (Mandatory) --samplename=
Strand information to be considered (Mandatory yes or no) --strand=
Feature bed file (Tab seperated) with no header. If --strand=yes, then the file should contain strand information in the fourth column (chrom,start,end,strand) (Mandatory) --bedfile=
Number of bases for avaeraging all the feature (Mandatory eg 1000) --avglength=
Number of bases to add as a flanking region (Mandatory eg 500) --flanking=
Threshold for genecount (Mandatory eg 1.5) --gthr=
Output Path (Mandatory) --outpath=

ii) Run the metaORFprocess.R with the tab seperated text files generated from previous step:
The function is to average metaORF values for a given length and create RData file

optional arguments:
  -h, --help         show this help message and exit
  -i INPSAMPLES      Input directory where all the bedgraph files are present (Mandatory)
  -sm SAMPLENAMES    comma seperated sample names
  -st ST             strand yes or no
  -ft FT             Average number of bases for feature Default = 1000
  -fb FLANKINGBASES  Flanking number of bases for feature Default = 500
  -pr PREFIX         Prefix tag for output (RData) file (Mandatory)
  -dir OUTDIR        Output directory path (Mandatory)
  
iii) Plot using the RData file generated from previous step:
This script is for plotting metaORF values

optional arguments:
  -h, --help   show this help message and exit
  -mf INMET    Comma seperated metaorf data files (RData file) (Mandatory)
  -tg TAGS     semi colan (;) seperated within the metaorf RData file and comma (,) seperated between the metaorf files (Mandatory)
  -nt NEWT     semi colan (;) seperated within the metaorf RData file and comma (,) seperated between the metaorf files (Mandatory)
  -sf STAG     suffix tag for plot (Mandatory)
  -ma MAIN     Main Title (Mandatory)
  -dir OUTDIR  Output directory path (Mandatory)
