Replication-Transcription Collision Analysis

This tool annotate the genes for collision (headon or codirectional) type based on the gene strand orientation and the nearby replication origin (present in either upstream or downstream of the gene). The center point is calculated for the replication origin and the nearby genes with varying distance (250b, 500b, 1kb, 2kb, 5kb, 10kb) from the origin center are used for the annotation

R Code: 

reptrans_collision.R

Required R Libraries: 

argparse, bedr, ggplot2

Input Files: 

Genes bed annotation file, Replication Origin bed annotation file, User input bed files for collision analysis. 
(Note: Example bed files are present in the libraryfiles folder)

Input Parameters:

  -sc SC      bed file containing gene information with strand information, third column is strand information (Mandatory)
  
  -rep REP    Replication Origin bed file, if not provided then -arep file should be given
  
  -arep AREP  Annotated Replication Origin bed file, if not provided then -rep should be provided
  
  -bc BC      comma seperated bed sample file path (Mandatory)
  
  -sn SN      sample name comma seperated for collisiontype (Mandatory)
  
  -pr PR      Prefix collisiontype for output file (Mandatory)
  
  -dir DIR    Output directory path (Mandatory)

Output:

  Replication Transcription Collision Gene Percentage for different distances (250b, 500b, 1kb, 2kb, 5kb, 10kb)
  
  ![ReplicationTranscriptionCollision](https://user-images.githubusercontent.com/18418058/57063280-0633ee80-6cc3-11e9-916b-3374114c8483.jpeg)
  
  Top2 protein peaks base percentage in the annotated collision genes (Data from GEO: GSE16258)
  
  ![Top2_RepTrans_Collision_Percentage](https://user-images.githubusercontent.com/18418058/57063346-35e2f680-6cc3-11e9-9701-b4e0b867864e.jpeg)
 
  Hmo1 protein peaks base percentage in the annotated collision genes (Data from GEO: GSE16258)
  
  ![Hmo1_RepTrans_Collision_Percentage](https://user-images.githubusercontent.com/18418058/57063391-5d39c380-6cc3-11e9-8dc2-2f3fa005a8ae.jpeg)
   
   Text files: 
   
   (i) allgene_reptranscollision.txt, ReplicationTranscriptionCollision_aggregate.txt contains the gene annotation with the collision type
   
   (ii) RepTrans_Collision_allgenes.txt, RepTrans_Collision_aggregate.txt contains the peak (User bed file) base percentage overlap in each collision type
