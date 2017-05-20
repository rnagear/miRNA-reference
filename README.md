# **miRNA-reference**

Files in this repository contain support tables and scripts used in the generation
of the reference table which serves as the basis of the metaMIR algorithm. 

**metaMIR** is a framework to predict in human interactions between 
microRNAs (miRNA) and clusters of genes. The user provides a set of genes to be targeted, 
and optionally genes not to be targeted. The analysis is performed to identify 
miRNAs that may simultaneously interact with a number of genes.

As the basis of these predictions, various machine learning algorithms were
used to integrate predictions from multiple other previously established
miRNA:target gene predictions (updated 2013-2016) which were available for download:
  * Diana micro-CDS v4.0 (Sep 2013 release) [here](http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=microT_CDS/index)
  * miRMap (version 201301) [here](http:77cegg.unige.ch/mirmap)
  * miRTarget 3.0 (prediction algorithm of miRDB v5.0) [here](http://mirdb.org/download.html)
  * MIRZA-G (note web access may be discontinued as of 30 Jun 2017) [here](http://www.clipz.unibas.ch/index.php?r=tools/mirza/Submission/index)
  * PACCMIT-CDS [here](http://paccmit.epfl.ch/)
  * TargetScan 7.0 [here](http://www.targetscan.org/vert_70/)

The miRNA:target data these algorithms were processed to align gene and miRNA names and identifiers,
and were trimmed to remove deprecated Entrez gene IDs (Ensembl Genes 83, GRCh38.p5). The predictions 
were used to classify miRNA:target gene interactions as positive (predicted to occur) or negative. 
The scores reported by the source algorithms are used as training features for supervised learning 
via a selection of classification algorithms.

This repository contains the script used to process the input data into a single R list object.
Note that execution of this script will require the source databases to be separately obtained from
the sources indicated above. The required directory structure for the location of the source
files is indicates in the `fread` statements in the script where the source files are read. Further
information is provided in comments in the script.

Also contained are flat text files derived from sets of validated miRNA:target gene interactions
consisting of:
  * the complete collection of validated interactions including the source of the experimental data
  * the training set used for machine learning
  * a test set of interactions that were unseen during training, used for subsequent performance evaluation
  
