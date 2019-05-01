This tool annotate the intergenic spaces of the whole genome based on the upstream and downstream gene directionality into convergent, divergent and codirectional groups. Then calculate the base percentage of given bed file into those groups. 

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
  
  
  

