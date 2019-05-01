Gene Pairs Directionality Analysis

This tool annotate the intergenic spaces of the whole genome based on the upstream and downstream gene directionality into convergent, divergent and codirectional± groups. Then calculate the base percentage of given bed file for the annotated intergenic space groups. 

Code: genedirectionality.R 

Required R libraries:

argparse
bedr
ggplot2

optional arguments:
  -h, --help  show this help message and exit
  -sc SC      bed file containing gene information with strand information, third column is strand information, if not provided then -dc should be provided
  -dc DC      Directionality annotated file, if not provided then -sc should be provided
  -bc BC      comma seperated bed sample file path
  -sn SN      sample name comma seperated for tag
  -pr PR      Prefix tag for output file (Mandatory)
  -dir DIR    Output directory path (Mandatory)
  
 Intergenic spaces of the whole genome based on the gene directionality in yeast (saccer2011)
 
 ![Genearrangement](https://user-images.githubusercontent.com/18418058/57011196-dc19f800-6c00-11e9-8d12-952a97d13fbb.jpeg)
  
 Hmo1 protein peak base percentage in the annotated groups convergent, divergent and codirectional± (GSE16258)
 
 ![Hmo1_directionalityIntergenicspaces](https://user-images.githubusercontent.com/18418058/57011331-74b07800-6c01-11e9-9d1e-441b3eaf65c9.jpeg)
 
 Top2 protein peak base percentage in the annotated groups convergent, divergent and codirectional± (GSE16258)
 
 ![Top2_directionalityIntergenicspaces](https://user-images.githubusercontent.com/18418058/57011421-d7097880-6c01-11e9-956a-f144f1f67aa9.jpeg)

  

