# check for installed packages. Install missing as needed
pkgs = c("data.table","ggplot2","RColorBrewer","corrplot","stringr")
if(length(new.pkgs <- setdiff(pkgs, rownames(installed.packages())))>0) install.packages(new.pkgs)
rm(pkgs,new.pkgs)
require(data.table)
# require(ggplot2)
# require(corrplot)
# require(gridExtra)
require(stringr)

# run scripts for ID (Ensembl, Entrez, miRNA alias) conversion
# note that the supporting reference tables must be in a subfolder called
# mapping under the current folder
source("conversion_toolsU.R")



################################################################################
# load TargetScan 7.0 data. Two files are provided of Context++ scores. A
# smaller files of conserved sites and a much larger file of non-conserved.
# Process both separately.
#
# redux with data.table fread. Do not need all columns and some need to be renamed
# to remove white spaces that will make referencing challenging later
# source file contains 12 columns:
# Classes 'data.table' and 'data.frame':  38309978 obs. of  12 variables:
# $ Gene ID                            : chr  "ENSG00000004059.6" "ENSG0...
# $ Gene Symbol                        : chr  "ARF5" "ARF5" "ARF5" "ARF5" ...
# $ Transcript ID                      : chr  "ENST00000000233.5" "ENST000...
# $ Gene Tax ID                        : int  9606 9606 10090 13616 9544 9...
# $ miRNA                              : chr  "hsa-miR-548az-5p" "hsa-miR-...
# $ Site Type                          : int  2 2 2 2 2 2 3 3 1 2 ...
# $ UTR_start                          : int  95 95 435 425 388 392 272 27...
# $ UTR end                            : int  101 101 441 431 394 398 279 ...
# $ context++ score                    : chr  "-0.0800" "-0.0800" "-0.2020...
# $ context++ score percentile         : chr  "91" "91" "96" "98" ...
# $ weighted context++ score           : chr  "-0.0800" "-0.0800" "-0.0010...
# $ weighted context++ score percentile: chr  "94" "94" "57" "56" ...
# - attr(*, ".internal.selfref")=<externalptr>
#
# retain only columns 1 (GeneID), 2 (GeneSym), 4 (Taxon), 5 (miRID),
# 6 (SiteType), 9 (CS), 11 (Score)
print("Load TargetScan, Non-Conserved data")
tscolnames <- c("ENSG","GeneSym","Taxon","miRNA","SiteType","CS","Score")
tscan <- fread("hsa/TargetScan7_hsa/Nonconserved_Site_Context_Scores.txt",
               select = c(1,2,4,5,6,9,11), col.names = tscolnames,
               na.strings="NULL")
tscan <- tscan[Taxon==9606 & !(is.na(CS) & is.na(Score))]
tscan <- tscan[, Taxon := NULL]
tscan[, ENSG := substring(ENSG,1,15)]

# convert Ensembl IDs
setkey(tscan,ENSG)
setkey(lookEnsg,Ensembl)
convE_tscan <- lookEnsg[tscan,nomatch=0]

# update miRNA IDs, ensuring miRBase v20
setkey(convE_tscan,miRNA)
setkey(lookmiRN,Alias)
convEM_tscan <- lookmiRN[convE_tscan,nomatch=0]

# remove duplicates, grouping by GeneName and miRNA taking best (min)
# weighted context++ score and retaining site type
targetscan_NC <- copy(convEM_tscan[,.(Score=sum(Score)),by = .(GeneName,miRID)])

#repeat for the much smaller conserved file]
print("Load TargetScan, Conserved data")
tscan <- fread("hsa/TargetScan7_hsa/Conserved_Site_Context_Scores.txt",
               select = c(1,2,4,5,6,9,11), col.names = tscolnames,
               na.strings="NULL")
tscan <- tscan[Taxon==9606 & !(is.na(CS) & is.na(Score))]
tscan <- tscan[, Taxon := NULL]
tscan[, ENSG := substring(ENSG,1,15)]

# convert Ensembl IDs
setkey(tscan,ENSG)
setkey(lookEnsg,Ensembl)
convE_tscan <- lookEnsg[tscan,nomatch=0]

# update miRNA IDs, ensuring miRBase v20
setkey(convE_tscan,miRNA)
setkey(lookmiRN,Alias)
convEM_tscan <- lookmiRN[convE_tscan,nomatch=0]

# remove duplicates, grouping by GeneName and miRNA taking best (min)
# weighted context++ score and retaining site type
targetscan_C <- copy(convEM_tscan[,.(Score=sum(Score)),by = .(GeneName,miRID)])

rm(convE_tscan,convEM_tscan,tscan,tscolnames)
gc(verbose=F)

################################################################################
# load PACCMIT
print("Load PACCMIT data")
temp <- fread("hsa/paccmit/predictions_human.txt",
              col.names=c("ENST","cons","Score","rawmiRs"))

setkey(temp,ENST)
setkey(lookEnst,Ensembl)
temp <- lookEnst[temp,nomatch=0]
rawmiRs <- strsplit(temp[["rawmiRs"]],split=",")
maxLen <- max(sapply(rawmiRs,length))
miRs <- t(sapply(rawmiRs,function(x){c(x,rep(NA,maxLen-length(x)))}))
temp2 <- do.call(cbind,list(temp[,.(GeneName,cons,Score)],miRs))
temp <- melt(temp2,id=c("GeneName","cons","Score"),value.name="miRNA",na.rm=T)
paccmit <- temp[,variable := NULL]
setkey(paccmit,miRNA)
setkey(lookmiRN,Alias)
paccmit <- lookmiRN[paccmit,nomatch=0]

setkey(paccmit,GeneName,miRID)
paccmit <- paccmit[,.(Score = min(Score)),by=.(GeneName,miRID)]
rm(temp, temp2, maxLen, rawmiRs, miRs)
gc(verbose=F)

################################################################################
# load MIRZA
# here, data provided in two separate tab delimited files, with and without
# seed consideration: mirza-g_all_mirnas_per_gene_scores.tab.tar.gz
# and seed-mirza-g_all_mirnas_per_gene_scores.tab.tar.gz
# unpack and begin with (slightly) smaller seed file
print("Load MIRZA-G data")
temp <- fread("hsa/mirza/seed-mirza-g_all_mirnas_per_gene_scores.tab",
              col.names=c("GeneID","miRNA","Score","ScoreNC"))

setkey(lookEntr,Entrez)
setkey(temp,GeneID)
convE_mirza <- lookEntr[temp,nomatch=0]
setkey(convE_mirza,miRNA)
setkey(lookmiRN,Alias)
convEM_mirza <- lookmiRN[convE_mirza,nomatch=0]

setkey(convEM_mirza,GeneName,miRID)
mirza_S <- convEM_mirza[,.(Score = max(Score)),by=.(GeneName, miRID)]
rm(temp,convE_mirza,convEM_mirza)
temp <- fread("hsa/mirza/mirza-g_all_mirnas_per_gene_scores.tab",
              col.names=c("GeneID","miRNA","Score","ScoreNC"))

setkey(lookEntr,Entrez)
setkey(temp,GeneID)
convE_mirza <- lookEntr[temp,nomatch=0]
setkey(convE_mirza,miRNA)
setkey(lookmiRN,Alias)
convEM_mirza <- lookmiRN[convE_mirza,nomatch=0]

setkey(convEM_mirza,GeneName,miRID)
mirza_NS <- convEM_mirza[,.(Score = max(Score)),by=.(GeneName, miRID)]

rm(convE_mirza,convEM_mirza,temp)
gc(verbose=F)

################################################################################
# load MirTarget3 (source data behind miRDB v5.0) released Aug 2014
# downloaded flat file,
mirtarget <- fread("hsa/mirtarget3/miRDB_v5.0_prediction_result.txt",
                   col.names=c("miRNA","RefID","Score"))
mirtarget <- mirtarget[grep("hsa",miRNA)]
mirtarget <- mirtarget[grep("NM",RefID)]
setkey(mirtarget,miRNA)
setkey(lookmiRN,Alias)
mirtarget <- lookmiRN[mirtarget,nomatch=0]
# lookRef2 <- unique(lookRef[,c(3:4),with=F])
setkey(mirtarget,RefID)
setkey(lookRef,Refseq)
mirtarget <- lookRef[mirtarget,nomatch=0,allow.cartesian=T]
mirtarget <- mirtarget[,.(GeneName,miRID,Score)]
mirtarget[,"Score":=max(Score),by=.(GeneName,miRID)]
setkey(mirtarget,GeneName,miRID)
mirtarget <- unique(mirtarget)

################################################################################
# load miRMap
# downloaded flat file, mirmap201301e_homsap_targets_1to1.csv.xz
# from http:77cegg.unige.ch/mirmap
# file structure (27 columns):
# Classes 'data.table' and 'data.frame':  22932919 obs. of  27 variables:
# $ mirna_id            : int  11021 11023 11030 11031 11039 ...
# $ mature_name         : chr  "hsa-let-7a-3p" "hsa-let-7b-3p...
# $ transcript_id       : int  24705 24705 24705 24705 24705 ...
# $ transcript_stable_id: chr  "ENST00000373020" "ENST0000037 ...
# $ transcript_canonical: chr  "t" "t" "t" "t" ...
# $ transcript_chr      : chr  "X" "X" "X" "X" ...
# $ transcript_strand   : int  -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
# $ gene_stable_id      : chr  "ENSG00000000003" "ENSG0000000...
# $ gene_name           : chr  "TSPAN6" "TSPAN6" "TSPAN6" "TS...
# $ name                : chr  "TSPAN6" "TSPAN6" "TSPAN6" "TS...
# $ site_nb             : int  2 2 2 2 2 1 1 1 1 1 ...
# $ seed6_nb            : int  1 1 1 2 1 1 1 1 1 1 ...
# $ seed7_nb            : int  1 1 1 0 1 0 0 0 0 0 ...
# $ tgs_au              : num  0.682 0.682 0.682 0.672 0.732 ...
# $ tgs_position        : int  150 150 150 150 260 386 232 39...
# $ tgs_pairing3p       : num  2 2.5 2.5 2 3 1 3.5 1 2 3.5 ...
# $ dg_duplex           : num  -10.7 -7.6 -12.3 -10.7 -15.6 -...
# $ dg_binding          : num  -11.05 -8.08 -12.21 -10.37 -16...
# $ dg_duplex_seed      : num  -3.8 -3.8 -3.8 -2.9 -7.5 -2.5 ...
# $ dg_binding_seed     : num  -4.38 -4.38 -4.38 -3.34 -7.87 ...
# $ dg_open             : num  17.8 17.9 17.9 17.9 13.3 ...
# $ dg_total            : num  8.145 10.519 6.658 8.258 0.699...
# $ prob_exact          : num  0.171 0.171 0.171 0.171 0.145 ...
# $ prob_binomial       : num  0.158 0.158 0.158 0.377 0.208 ...
# $ cons_bls            : num  1.024 1.024 1.024 1.024 0.763 ...
# $ selec_phylop        : num  0.132 0.132 0.132 0.425 0.333 ...
# $ mirmap_score        : num  0.0774 0.071 0.0593 -0.0322 -0...
# - attr(*, ".internal.selfref")=<externalptr>
# retain cols mature_name (2), gene_stable_id (8), gene_name (9),
# mirmap_score (27)
print("Load miRMap data")
mmcolnames <- c("miRNA","ENSG","MNam","Score")
temp <- fread("hsa/miRMap/mirmap201301e_homsap_targets_1to1.csv",
              select = c(2,8,9,27), col.names = mmcolnames)
setkey(lookEnsg,Ensembl)
setkey(temp,ENSG)
convE_mirmap <- lookEnsg[temp,nomatch=0]
setkey(convE_mirmap,miRNA)
setkey(lookmiRN,Alias)
convEM_mirmap <- lookmiRN[convE_mirmap]

setkey(convEM_mirmap,GeneName,miRID)
mirmap <- convEM_mirmap[,.(Score = min(Score)),by=.(GeneName, miRID)]
mirmap <- mirmap[Score<0.05]
rm(convEM_mirmap,convE_mirmap,temp,mmcolnames)
gc(verbose=F)

################################################################################
# load Diana microT-CDS
# source data from 2013 (Diana microT-CDS, 2.1Gb file) apparently also form the
# base for the most recent version of the Taverna v5.0 plug-in
# with no fastreadline equivalent in data.table, use base R to read data. Process
# data to retain only rows with ENS string (to exclude rows of site location
# data) and resave as plain text. Can subsequently use fread to process

# if working from the original Diana file, uncomment this block to preprocess
# temp=readLines("hsa/Diana_microtv4_data/diana_microT_CDS_data.csv")
# temp=temp[grep("ENS",temp)]
# temp=append("ENST,ENSG_GenName,miRNA_ver,Score",temp,after=1)
# write.table(temp,
#             file="hsa/Diana_microtv4_data/processed_Diana_microT_CDS_data.csv",
#             quote=F,row.names=F,col.names=F)
print("Load Diana microT-CDS v4 data")
temp <- fread("hsa/Diana_microtv4_data/processed_Diana_microT_CDS_data.csv")
keepers <- grep("hsa",temp[["miRNA_ver"]])
temp_trim <- data.table(miRNA_ver = temp[["miRNA_ver"]][keepers])
for (i in colnames(temp)[c(2,4)]){
  temp_trim[,(i) := temp[[i]][keepers]]
}

temp_trim[,`:=`(ENSG = substring(ENSG_GenName,1,15),
                DNam = substring(ENSG_GenName,17,(nchar(ENSG_GenName)-1)))]
temp_trim[,`:=`(miRNA = substring(miRNA_ver,1,(nchar(miRNA_ver)-4)),
                mRver = substring(miRNA_ver,(nchar(miRNA_ver)-2),(nchar(miRNA_ver)-1)))]
temp <- temp_trim[,.(ENSG,DNam,miRNA,Score)]

setkey(temp,ENSG)
setkey(lookEnsg,Ensembl)
convE_diana <- lookEnsg[temp,nomatch=0]
setkey(convE_diana,miRNA)
setkey(lookmiRN,Alias)
convEM_diana <- lookmiRN[convE_diana,nomatch=0]

setkey(convEM_diana,GeneName,miRID)
diana <- convEM_diana[,.(Score = max(Score)),by=.(GeneName,miRID)]
rm(temp,temp_trim,convE_diana,convEM_diana,i,keepers)
gc(verbose=F)

################################################################################
# load prediction set from ComiRNet released Jun 2015
# downloaded flat file,
print("Load ComiRNet data")
comirnet <- fread("hsa/ComiRNet/ourDataset.txt",
              col.names = c("miRNA","Rawname","Score"))
# authors converted all entries in flat file to lowercase, including
# miRNA names. Substitute -miR- for -mir-
comirnet[,"miRNA":=str_replace(miRNA,pattern = "mir",replacement = "miR")]

setkey(comirnet,miRNA)
setkey(lookmiRN,Alias)
comirnet <- lookmiRN[comirnet,nomatch=0]
# the gene names in the ComiRNet dataset are all lowercase. Convert
# to uppercase, then check for presence in reference lookup tables
comirnet[,"Rawname":=toupper(Rawname)]
# comirnet names match better to Primary column in lookup table
# take only comirnet rows with a "Rawname" that appears in lookUnip's
# "Primary" column
setkey(lookUnip,Primary)
setkey(comirnet,Rawname)
comirnet <- unique(lookUnip[,.(Primary)])[comirnet,nomatch=0]
comirnet <- comirnet[,.(GeneName=Primary,miRID,Score)]
comirnet[,"Score":=max(Score),by=.(GeneName,miRID)]
setkey(comirnet,GeneName,miRID)
comirnet <- unique(comirnet)


################################################################################
# collect matrices in a single list object for ease of manipulation
print("Generate list object of input databases")

targetscan_U <- do.call(rbind,list(targetscan_C,targetscan_NC))
setkey(targetscan_U,GeneName,miRID)
targetscan_U[,Score := sum(Score),by=.(GeneName,miRID)]
targetscan_U <- unique(targetscan_U)

mirza_U <- do.call(rbind,list(mirza_S,mirza_NS))
setkey(mirza_U,GeneName,miRID)
mirza_U[,Score:=max(Score),by=.(GeneName,miRID)]
mirza_U <- unique(mirza_U)

dbase <- list(targetscan_C=copy(targetscan_C),targetscan_NC=copy(targetscan_NC),
              targetscan_U=copy(targetscan_U),mirza_S=copy(mirza_S),
              mirza_NS=copy(mirza_NS),mirza_U=copy(mirza_U),
              mirtarget=copy(mirtarget),mirmap=copy(mirmap),
              paccmit=copy(paccmit),diana=copy(diana))
# clean individual objects
rm(diana,mirmap,mirtarget,paccmit,targetscan_C,targetscan_NC,targetscan_U,
   mirza_NS,mirza_S,mirza_U)
for (i in 1:length(dbase)){
  dbase[[i]][,"dbase":=names(dbase)[i]]
  print(paste(names(dbase)[i],paste0(as.character(dim(dbase[[i]])))))
}

setkey(dbase[[1]],GeneName,miRID)
uni_set <- copy(dbase[[1]])
int_set <- copy(dbase[[1]])
int_set <- int_set[,.(GeneName,miRID)]
# create a long version of the data to collect prediction scores
for (i in 2:length(dbase)){
  dbase[[i]] <- dbase[[i]][!is.na(miRID)]
  setkey(dbase[[i]],GeneName,miRID)
  int_set <- do.call(merge,list(int_set,dbase[[i]][,.(GeneName,miRID)],by=c("GeneName","miRID")))
  uni_set <- do.call(rbind,list(uni_set,dbase[[i]]))
  print(paste(names(dbase)[i],as.character(dim(dbase[[i]])),sep=":"))
}
# write the dbase list object to file if necessary
if(!("dbase_object.RData" %in% list.files())){
  save(dbase,file="dbase_object.RData")
}
# remove dbase object
rm(dbase)
gc()

master_tab <- dcast.data.table(uni_set,GeneName + miRID ~ dbase,value.var="Score")
# write the master table object to file, if necessary
if(!("Master_raw_algo_table.RData" %in% list.files())){
  save(master_tab,file="Master_raw_algo_table.RData")
}

# variable and function cleanup
rm(list=ls(pattern="look"))
rm(int_set,uni_set,i)
gc()

##############################################################################
# To use the previous algorithm scores as input for machine learning models,
# need to deal with NAs, where a given algorithm (eg diana) has no prediction
# score for a given miRNA/GeneName pair

# load the master table object, if not already in memory
if (!exists("master_tab")){
  load("Master_raw_algo_table.RData")
}
master_tab.N <- copy(master_tab)
rm(master_tab)

# in final analysis, will separately consider targetscan C and NC, but
# due to significant overlap (Spearman correlation of 0.81) the union
# set of mirza will be used. Remove the targetscan_U and mirza_S, mirza_NS
master_tab.N[,c("targetscan_U","mirza_NS","mirza_S"):=NULL]
algo_bounds <- data.frame(matrix(nrow=(length(names(master_tab.N))-2),ncol=3))
colnames(algo_bounds) <- c("Algorithm","Min","Max")
for (i in 1:length(names(master_tab.N)[-c(1:2)])){
  algo_bounds$Algorithm[i] <- names(master_tab.N)[2+i]
  algo_bounds$Max[i] <- max(master_tab.N[,names(master_tab.N)[i+2],with=F],na.rm=T)
  algo_bounds$Min[i] <- min(master_tab.N[,names(master_tab.N)[i+2],with=F],na.rm=T)
}

# reorient mirmap,paccmit, and targetscan so that best score is the highest
# to keep positive, mirmap -> -mirmap, paccmit -> 1-paccmit, targetscan -> -targetscan
master_tab.N[,c("mirmap","paccmit","targetscan_C","targetscan_NC"):=
               list((-mirmap+0.05),(1.0-paccmit),(-targetscan_C),(-targetscan_NC))]
# floor the scores for each algorithm. Determine the spread and subtract
# one half of the spread from the minimum value. The result will replace
# the NA values in the array
for (colidx in names(master_tab.N)[-c(1:2)]){
  raw_dat <- master_tab.N[,colidx,with=F]
  max_dat <- max(raw_dat,na.rm=T)
  min_dat <- min(raw_dat,na.rm=T)
  repl <- min_dat-(max_dat-min_dat)/2
  raw_dat[is.na(raw_dat)] <- repl
  master_tab.N[,(colidx):=raw_dat]
}
rm(raw_dat,max_dat,min_dat,repl,i,colidx)
# resulting table is the "floored" data, where NAs have been replaced with a
# value 50% of the range lower than the "lowest" score (ie least likely to
# indicate miRNA targeting of a gene)
master_tab.Fl <- copy(master_tab.N)
# write the object to disk if not already saved
if(!("Range_floored_master_tab.RData" %in% list.files())){
  save(master_tab.Fl,file="Range_floored_master_tab.RData")
}

rm(master_tab.N)

