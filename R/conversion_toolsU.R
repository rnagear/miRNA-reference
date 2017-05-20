# check for installed packages. Install missing as needed
pkgs = c("data.table")
if(length(new.pkgs <- setdiff(pkgs, rownames(installed.packages())))>0) install.packages(new.pkgs)
rm(pkgs,new.pkgs)
require(data.table)

# from Ensembl biomart, select Genes 83, Homo Sapiens Genes GRCh38.p5
# under attributes, check Ensembl Gene ID, Associated Gene Name,
# filter by protein coding genes
lookEnsg <- fread("hsa/mapping/ensembl_GID_GeneName_mart_export_20160103.txt",
                  col.names=c("Ensembl","GeneName"))
setkey(lookEnsg,NULL)
lookEnsg <- unique(lookEnsg)

# from Ensembl biomart, select Genes 83, Homo Sapiens Genes GRCh38.p5
# under attributes, check Ensembl Transcript ID, Associated Gene Name,
# filter by protein coding transcripts
lookEnst <- fread("hsa/mapping/ensembl_TID_GeneName_mart_export_20160107.txt",
                  col.names=c("Ensembl","GeneName"))
setkey(lookEnst,NULL)
lookEnst <- unique(lookEnst)

# for Entrez IDs, visit http://www.genenames.org/cgi-bin/download
# check options for approved gene symbol, Entrez ID, and filter for those
# entries with Approved status
lookEntr <- fread("hsa/mapping/HGNC_AppSym_EntrezID_20160106.txt",
                  col.names=c("GeneName","Entrez"))
setkey(lookEntr,NULL)
lookEntr <- unique(lookEntr)
lookEntr <- lookEntr[!is.na(Entrez)]

# create lookup table for refseq IDs
lookRef <- fread("hsa/mapping/ENSG_ENST_Refseq_GN_mart_export_20161019.txt",
              col.names=c("ENSG","ENST","Refseq","GeneName"))
lookRef <- lookRef[Refseq!=""]
lookRef <- lookRef[grep("NM_",Refseq)]
setkey(lookRef,Refseq,GeneName)
lookRef <- unique(lookRef)
setkey(lookRef,Refseq)

# load UniProt data to check Primary Gene Names against these lists
# query on uniprot.org, taxonomy:9606, selecting "Human" as species
# (to exclude other data such as from neanderthal which is also caught
# a taxon:9606 query). Ensure columns for revision status,
# Gene Names (primary), and Gene Names (Synonym) are selected
lookUnip <- fread("hsa/mapping/uniprot-organism_Homo+sapiens+Human.tab",
                  col.names=c("UniProtID","RevStat","Aliases", "Primary"))
# trim entries with no Primary Gene Name
lookUnip <- lookUnip[Primary!=""]

# trim entries with slashes, spaces, semicolons or <> in the Primary Name
rejects <- grep(";",lookUnip[["Primary"]])
rejects <- union(rejects, grep(" ",lookUnip[["Primary"]]))
rejects <- union(rejects, grep("\\/",lookUnip[["Primary"]]))
rejects <- union(rejects, grep("<",lookUnip[["Primary"]]))
keepers <- setdiff(1:nrow(lookUnip),rejects)
Unip.sub <- data.table(Primary = lookUnip[["Primary"]][keepers])
colkeep <- colnames(lookUnip)[c(3,1,2)]
for (i in colkeep){
  Unip.sub[, (i) := lookUnip[[i]][keepers]]
}
lookUnip <- Unip.sub
rm(rejects,keepers,colkeep,i,Unip.sub)

# remove duplicates with keys on Primary and Aliases
setkey(lookUnip,Primary,Aliases)
lookUnip <- unique(lookUnip)

###############################################################################
# Generate conversion table for miRNA IDs. Visit database downloads at mirbase.org
# ftp site, ftp://mirbase.org/pub/mirbase/CURRENT/ and download the alias.txt.gz file. Rename
temp <- fread("hsa/mapping/aliases.txt",col.names=c("Acc","Aliases"))
temp <- temp[intersect(grep("hsa",Aliases),grep("MIMAT",Acc))]

temp2 <- readLines("hsa/mapping/mature.fa")
temp2 <- temp2[grep("hsa",temp2)]
temp2 <- gsub(">","",temp2)
temp2 <- data.table(t(sapply(temp2,function(x){unlist(strsplit(x,split=" "))})))
temp2[,c("V3","V4","V5") := NULL]
colnames(temp2) <- c("miRID","Acc")
lookmiRAc <- temp2

setkey(temp,Acc)
setkey(temp2,Acc)
temp <- temp[temp2]

temp[miRID=="hsa-miR-299-3p", Aliases := "hsa-miR-299;hsa-miR-299-3p"]
temp[miRID=="hsa-miR-34b-5p",  Aliases := "hsa-miR-34b*;hsa-miR-34b-5p"]
temp[miRID=="hsa-miR-378a-5p", Aliases := "hsa-miR-378*;hsa-miR-378a-5p"]
temp[miRID=="hsa-miR-340-3p",  Aliases := "hsa-miR-340*;hsa-miR-340-3p"]
temp[miRID=="hsa-miR-425-3p",  Aliases := "hsa-miR-425*;hsa-miR-425-3p"]
temp[miRID=="hsa-miR-488-5p",  Aliases := "hsa-miR-488*;hsa-miR-488-5p"]
temp[miRID=="hsa-miR-493-5p",  Aliases := "hsa-miR-493-5p;hsa-miR-493*"]
temp[miRID=="hsa-miR-516b-3p", Aliases := "hsa-miR-516b*;hsa-miR-516-3p;hsa-miR-516b-3p"]
temp[miRID=="hsa-miR-499b-5p", Aliases := "hsa-miR-499b-5p"]
temp[miRID=="hsa-miR-500a-3p", Aliases := "hsa-miR-500*;hsa-miR-500a*;hsa-miR-500a-3p"]
temp[miRID=="hsa-miR-589-3p",  Aliases := "hsa-miR-589*;hsa-miR-589-3p"]
temp[miRID=="hsa-miR-550a-3p", Aliases := "hsa-miR-550*;hsa-miR-550a*;hsa-miR-550a-3p"]
temp[miRID=="hsa-miR-593-5p",  Aliases := "hsa-miR-593*;hsa-miR-593-5p"]
temp[miRID=="hsa-miR-616-5p",  Aliases := "hsa-miR-616*;hsa-miR-616-5p"]
temp[miRID=="hsa-miR-624-5p",  Aliases := "hsa-miR-624*;hsa-miR-624-5p"]
temp[miRID=="hsa-miR-629-3p",  Aliases := "hsa-miR-629*;hsa-miR-629-3p"]
temp[miRID=="hsa-miR-499b-3p", Aliases := "hsa-miR-499b-3p"]

# tmp_miRs <- temp[["Aliases"]]
tmp_miRs <- strsplit(temp[["Aliases"]],split=";")
maxLen <- max(sapply(tmp_miRs,length))
tmp_miRs <- t(sapply(tmp_miRs,function(x){c(x,rep(NA,maxLen-length(x)))}))
temp <- do.call(cbind,list(temp[,.(Acc,miRID)],tmp_miRs))
temp <- melt(temp,id.vars=c("Acc","miRID"),value.name="Alias",na.rm=T)
temp[,c("variable","Acc") := NULL]
setkey(temp,NULL)
lookmiRN <- unique(temp)
rm(temp,temp2,maxLen,tmp_miRs)
