Replication-Transcription Collision Analysis

This tool annotate the genes for collision (headon or codirectional) type based on the gene strand orientation and the nearby replication origin (present in either upstream or downstream of the gene). The center point is calculated for the replication origin and the nearby genes with varying distance (250b, 500b, 1kb, 2kb, 5kb, 10kb) from the origin center are used for the annotation

R Code: reptrans_collision.R

Required Libraries: argparse, bedr, ggplot2

Input Parameters:
  -sc SC      bed file containing gene information with strand information, third column is strand information (Mandatory)
  
  -rep REP    Replication Origin bed file, if not provided then -arep file should be given
  
  -arep AREP  Annotated Replication Origin bed file, if not provided then -rep should be provided
  
  -bc BC      comma seperated bed sample file path (Mandatory)
  
  -sn SN      sample name comma seperated for collisiontype (Mandatory)
  
  -pr PR      Prefix collisiontype for output file (Mandatory)
  
  -dir DIR    Output directory path (Mandatory)
