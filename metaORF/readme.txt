
MetaORF or metaFeature to study the averaged enriched peak profile across all protein coding genes or specific set of features along with upstream and downstream profile. Python scriot is used for average profile calculation and Rscript is used for plotting.

(Note: please use python3 to run this tool)

Required python3 packages:

argparse, pybedtools, pandas, multiprocessing

Required R packages:

argparse, plyr, ggplot2, gridExtra

Input files:

(i) protein coding bed annotation file, (ii) user peak bed and bedgraph graph files

Steps:
i) Run the metaFeature.py file to generate metaFeature.tsv which contains intensity average and gene density across the average bases.

Parameters:

  -bdg BDGFILE, --bdg BDGFILE
                        bedGraph file with no header (Mandatory)
  
  -bed BEDFILE, --bed BEDFILE
                        bedfile with header (Mandatory)
  
  -ft FEATUREFILE, --feature FEATUREFILE
                        Feature bed file with header, four columns are mandatory with the column names (chrom,start,end,name) (Mandatory)
  
  -st STRAND, --strand STRAND
                        Value:[yes or no] If yes, normalize the strand orientaion for aggregation (Default = no)
  
  -stn STRNAME, --strandname STRNAME
                        Column name for strand from feature file, column should contain (+ or 1, - or -1) (Default=None)
  
  -bs BASESMERGING, --basesmerging BASESMERGING
                        Total base for feature aggragation (Mandatory)
  
  -fl FLANKING, --flanking FLANKING
                        Value:[yes or no] If yes, flanking region will be considered (Default = no)
  
  -flb FLANKINGBASES, --flankingbases FLANKINGBASES
                        Total flanking bases. (Required only if -fl is yes)
  
  -sm SAMPLENAME, --samplename SAMPLENAME
                        Sample name
  
  -sb SUBTYPE, --subtype SUBTYPE
                        Subtype column name from feature file, where subtype is used for aggregation (Default=None)
  
  -t PROCS, --multithreads PROCS
                        Comma seperated output prefix (Mandatory). Same number of items as -f argument (Default=2)
  
  -o OUTDIR, --outdir OUTDIR
                        outdir (Mandatory)


ii) Plot with script plot_metaFeature.R using the metaFeature.tsv files.

Parameters:

  -h, --help   show this help message and exit

  -mf INMET    Comma seperated text files (Mandatory)
  
  -tg TAGS     semi colan (;) seperated within the metaorf txt file and comma (,) seperated between the metaorf files
               (Example:SampleASubtype1;SampleASubtype2,SampleBSubtype1;SampleBSubtype2), If not provided all the subtypes will be used
  
  -nb NBASES   Number of bases used for merging (Mandatory)
  
  -sf STAG     suffix tag for plot (Mandatory)
  
  -ma MAIN     Main Title (Mandatory)
  
  -dir OUTDIR  Output directory path (Mandatory)
