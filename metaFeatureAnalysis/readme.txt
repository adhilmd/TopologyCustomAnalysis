
MetaFeature Analysis

MetaORF or metaFeature to study the averaged enriched peak profile across all protein coding genes or specific set of features along with upstream and downstream profile. Python script is used for average profile calculation and Rscript is used for plotting.

(Note: please use python3 to run this tool)

Required python3 packages:

argparse, pybedtools, pandas, multiprocessing

Required R packages:

argparse, plyr, ggplot2, gridExtra

Input files:

(i) protein coding bed annotation file, (ii) user peak bed and bedgraph file

Steps:
i) Run the metaFeature.py file using python3 to generate metaFeature.tsv, which contains intensity average and gene count across the averaged bases.

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
  
  -avgt {median,mean}, --avgtype {median,mean}
                        Ratio averaging type mean or median (Default = median)
  
  -sm SAMPLENAME, --samplename SAMPLENAME
                        Sample name
  
  -sb SUBTYPE, --subtype SUBTYPE
                        Subtype column name from feature file, where subtype is used for aggregation (Default=None)
  
  -t PROCS, --multithreads PROCS
                        Comma seperated output prefix (Mandatory). Same number of items as -f argument (Default=2)
  
  -o OUTDIR, --outdir OUTDIR
                        outdir (Mandatory)


ii) Plot the averaged values (metaFeature.tsv) using plot_metaFeature.R.

Parameters:

  -h, --help   show this help message and exit

  -mf INMET    Comma seperated text files (Mandatory)
  
  -tg TAGS     semi colan (;) seperated within the metaorf txt file and comma (,) seperated between the metaorf files
               (Example:SampleASubtype1;SampleASubtype2,SampleBSubtype1;SampleBSubtype2), If not provided all the subtypes will be used
 
  -rtype {l2fc,ratio}  ratio type l2fc=log2(IP/Input) or just
                       ratio=(IP/Input), if l2fc then exponential function is
                       used to convert all values to positive scale (Default =
                       l2fc)
  -sm {yes,no}         smoothing the curve using geom_smooth function
                       (Default=yes)
  -cf {yes,no}         95 percent confidence interval to be plotted
                       (Default=no)
                       
  -nb NBASES   Number of bases used for merging (Mandatory)
  
  -sf STAG     suffix tag for plot (Mandatory)
  
  -ma MAIN     Main Title (Mandatory)
  
  -dir OUTDIR  Output directory path (Mandatory)
  

Output Files:

MetaORF Gene Count plot where the y axis contains the average number of genes containing the peaks across the upstream, ORF and downstream

![MetaFeature_Count_Top2Hmo1](https://user-images.githubusercontent.com/18418058/57066718-41d3b600-6ccd-11e9-9801-8b1fefc46d33.jpeg)

MetaORF Gene Intensity plot where the y axis contains the average intensity across the upstream, ORF and downstream

![MetaFeature_Score_Top2Hmo1](https://user-images.githubusercontent.com/18418058/57066719-41d3b600-6ccd-11e9-9af5-18f3caf788c1.jpeg)

Table: metaFeature.tsv containg all the averaged values (Intesity and Gene Count)
