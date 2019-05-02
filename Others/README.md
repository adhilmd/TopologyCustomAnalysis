Peak analysis scripts
(Note: All the library and sample input files are in the Libraryfiles folder)

(i) Plotting user defined chromosome tracks using bedgraph files

Code: chromosome_trackplot.R

Required R Libraries: Gviz, GenomicRanges

Input Files: Protein coding bed annotation file, user bed and bedgraph files
(Note: Only four files are allowed per plot)

Parameters:

  -chr CHR    chromosome number to plot (Mandatory)
  
  -st ST      start site (Mandatory)
  
  -en EN      end site (Mandatory)
  
  -bdg BDG    Comma seperated bedgraph files, only four allowed (Mandatory)
  
  -bd BD      Comma seperated bed files, only four allowed and also should match the number of bedgraph files (Mandatory)
  
  -sc SC      Protein coding bed annotation file (Mandatory)
  
  -tg TG      Comma seperated tags, should match the number of bed and bedgraph files (Mandatory)
  
  -pr PR      Prefix for output file (Mandatory)
  
  -ma MA      Main Title for the plot
  
  -dir DIR    Output directory path (Mandatory)
  
 Output file:
 
 Plot containing the peak bar and peak values with gene annotation
 
 ![Top2Hmo1_genes_track](https://user-images.githubusercontent.com/18418058/57076938-0f38b600-6ceb-11e9-909d-3b0d9bad3e51.jpeg)

(ii) Correlation analysis between two bed files based on the feature overlap (feature can be Genes, ReplicationOrigin and telomeres etc)

Code: feature_peak_correlation.R

Required R Libraries: argparse, bedr, ggplot2

Input Files: Protein coding bed annotation file, two user bed files to compare the correlation (sample1VSsample2)

Parameters:

  -sc SC      bed file containing gene information with strand information, third column is strand information (Mandatory)
  
  -ud UD      Upstream and Downstream flanking bases to include
  
  -bd BD      comma seperated two bed files path for sample1VSsample2 (Mandatory)
  
  -tg TG      comma seperated tags, should be 2 tags sample1 and sample2 (Mandatory)
  
  -pr PR      Prefix for output file (Mandatory)
  
  -dir DIR    Output directory path (Mandatory)
 
 Output files:
 
 File contating the statistical measures (correlation value, pvalue) and correlation plot
 
 ![Top2Hmo1_correlation](https://user-images.githubusercontent.com/18418058/57076941-12cc3d00-6ceb-11e9-9c2c-ba60a3e208b2.jpeg)

 
 (iii) Peak annotation for given feature set
 
 Code: proteincoding_annotation.R
 
 Required R Libraries: argparse, bedr
 
 Input Files: Protein coding bed annotation file and bed file
 
 Parameters:
 
  -sc SC      Protein coding annotation bed file (Mandatory)
  
  -bd BD      bed file (Mandatory)
  
  -udb UDB    Upstream and Downstream bases for annotation (Mandatory)
  
  -pr PR      Prefix for output file (Mandatory)
  
  -dir DIR    Output directory path (Mandatory)
 
 Output files:
 
 Output annotated peak file and pie chart showing the percentage of annotated peaks

Top2

![Top2_updown500_annotation](https://user-images.githubusercontent.com/18418058/57076967-1c55a500-6ceb-11e9-99f9-a6963e372c3d.jpeg)

Hmo1
![Hmo1_updown500_annotation](https://user-images.githubusercontent.com/18418058/57076976-1f509580-6ceb-11e9-9be1-3a9bb5a09fb5.jpeg)
