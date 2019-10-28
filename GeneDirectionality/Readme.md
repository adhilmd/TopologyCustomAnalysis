Gene Pairs Directionality Analysis

This tool annotate the intergenic spaces of the whole genome based on the upstream and downstream gene directionality into convergent, divergent and codirectional± groups. Then calculate the base percentage of given bed file for the annotated intergenic space groups. 

Code: genedirectionality.R 

Required R libraries: argparse, bedr, ggplot2

Input Files: Protein coding annotation file and bed files  

Arguments:

-sc SC bed file containing gene information with strand information, third column is strand information, if not provided then -dc should be provided

-dc DC Directionality annotated file, if not provided then -sc should be provided

-bc BC comma seperated bed sample file path (Mandatory)

-sn SN sample name comma seperated for tag (Mandatory)

-pr PR  Prefix tag for output file (Mandatory)

-dir DIR Output directory path (Mandatory)
  
Output files:

 Intergenic spaces of the whole genome based on the gene directionality in yeast (saccer2011)
 
 ![Genearrangement](https://user-images.githubusercontent.com/18418058/57011196-dc19f800-6c00-11e9-8d12-952a97d13fbb.jpeg)
  
 Hmo1 protein peak base percentage in the annotated groups convergent, divergent and codirectional± (GSE16258)
 
 ![Hmo1_directionalityIntergenicspaces](https://user-images.githubusercontent.com/18418058/67693674-6241ac00-f9a2-11e9-908a-c67818b1f9fe.jpeg)
 
 Top2 protein peak base percentage in the annotated groups convergent, divergent and codirectional± (GSE16258)
 
  ![Top2_directionalityIntergenicspaces](https://user-images.githubusercontent.com/18418058/67693771-92894a80-f9a2-11e9-8e43-5b96d7240306.jpeg)

