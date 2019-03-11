getgCSI <-
  function (verbose=FALSE,
            nthread=1){
options(stringsAsFactors = FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[()]|[)]|[()]|[/]"

gCSI_V1_2_raw <- read.csv("/pfs/getgcsi_2018/gCSI_GRvalues_v1.2.tsv", sep="\t")
#gCSI_v1_1_metric <- read.csv("../data/gCSI_GRmetrics_v1.tsv", sep="\t")
gCSI_v1_2_metric <- read.csv("/pfs/getgcsi_2018/gCSI_GRmetrics_v1.2.tsv", sep="\t")

cell_all <- read.csv("/pfs/getgcsi_2018/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
drug_all <- read.csv("/pfs/getgcsi_2018/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

curationCell <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationTissue <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "gCSI.cellid")]
curationTissue <- curationTissue[ , c("unique.tissueid", "gCSI.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


####drug matching####

matchToIDTableDRUG <- function(ids,tbl, column) {
  sapply(ids,function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating drug ids")
    }
    return(tbl[myx, "unique.drugid"])
  })
}

matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}



druginfo <- unique(gCSI_V1_2_raw$DrugName)
#min(length(drugid), length(drugs)) - length(intersect(drugs, drugid))

curationDrug_Temp <- curationDrug
curationDrug_Temp$gCSI.drugid <- tolower(gsub(badchars, "", curationDrug$gCSI.drugid))
druginfo_Temp <- tolower(gsub(badchars, "", druginfo))
y <- curationDrug[which(!curationDrug_Temp$gCSI.drugid %in% druginfo_Temp),]
x <- as.character(matchToIDTableDRUG(druginfo_Temp, curationDrug_Temp, "gCSI.drugid"))
druginfo <- data.frame("drugid"=c(druginfo, y$gCSI.drugid), "unique.id" = c(x, y$unique.drugid))
rownames(druginfo) <- druginfo$unique.id



celline.info <- unique(gCSI_V1_2_raw$CellLineName)
#ff <- gsub("HCC70 ", "HCC70", gCSI_V1_2_raw$CellLineName)
curationCell_Temp <- curationCell
curationCell_Temp$gCSI.cellid <- tolower(gsub(badchars, "", curationCell$gCSI.cellid))
celline.info_Temp <- tolower(gsub(badchars, "", celline.info))
celline.info[which(celline.info_Temp %in% celline.info_Temp[which(duplicated(celline.info_Temp))])]

x2 <- as.character(matchToIDTableCELL(celline.info_Temp, curationCell_Temp, "gCSI.cellid"))
y2 <- curationCell[which(!curationCell_Temp$gCSI.cellid %in% celline.info_Temp),]
celline.info <- data.frame("cellid"=c(celline.info, y2$gCSI.cellid), "unique.id" = c(x2, y2$unique.cellid))
rownames(celline.info) <- celline.info$unique.id

celldata <- read.csv("/pfs/getgcsi_2018/cellinfo_old.csv", na.strings=c("", " ", "NA"), row.names = 1)
celldata[,"CellLineName"] <- curationCell$gCSI.cellid[match(rownames(celldata),curationCell$unique.cellid)]
missing_cells <- celline.info$unique.id[which(!celline.info$unique.id %in% rownames(celldata))]
missingdf <- data.frame(matrix(NA, nrow=length(missing_cells),ncol = ncol(celldata)))
colnames(missingdf) <- colnames(celldata)
rownames(missingdf) <- missing_cells
missingdf$CellLineName <- celline.info$cellid[match(rownames(missingdf), celline.info$unique.id)]
missingdf$TissueMetaclass <- gCSI_V1_2_raw$PrimaryTissue[match(tolower(gsub(badchars, "", missing_cells)), tolower(gsub(badchars, "", gCSI_V1_2_raw$CellLineName)))]
missingdf$tissueid <- curationTissue$unique.tissueid[match(missing_cells,rownames(curationTissue))]
missingdf$unique.id <- rownames(missingdf)
celldata <- rbind(celldata,missingdf)
celline.info <- celldata
#gCSI_V1_2_raw$CellLineName <- gsub("HEC-1-B", "HEC-1",gCSI_V1_2_raw$CellLineName)
#gCSI_V1_2_raw$CellLineName <- gsub("KNS-81", "KNS-81-FD",gCSI_V1_2_raw$CellLineName)
u <- as.character(matchToIDTableCELL(tolower(gsub(badchars, "", unique(gCSI_V1_2_raw$CellLineName))), curationCell_Temp, "gCSI.cellid"))
#celldata$DoublingTime[row.names(celldata) %in% u] <- gCSI_V1_2_raw$doublingtime[which(isolated_cell %in% gCSI_V1_2_raw$CellLineName)]
isolated_cell <- celldata$CellLineName[row.names(celldata) %in% u]
celldata$DoublingTime[row.names(celldata) %in% u] <- gCSI_V1_2_raw$doublingtime[match(isolated_cell, gCSI_V1_2_raw$CellLineName)]



load("/pfs/recompGCSI2018/gCSI_auc_v1_2.RData")
colnames(gCSI_AUC)[4] <- "auc_recomputed"
colnames(gCSI_AUC)[5] <- "rep_no"
gCSI_AUC <- cbind(gCSI_AUC, "cellid"=NA, "drugid"=NA)
gCSI_AUC$drugid <- matchToIDTableDRUG(tolower(gsub(badchars, "", gCSI_AUC$DrugName)), curationDrug_Temp, "gCSI.drugid")
gCSI_AUC$cellid <- matchToIDTableCELL(tolower(gsub(badchars, "", gCSI_AUC$CellLineName)), curationCell_Temp, "gCSI.cellid") #tolower has to be used for "HCC-70 " spacing
rownames(gCSI_AUC) <- sprintf("%s_%s_%s", gCSI_AUC$cellid, gCSI_AUC$drugid, gCSI_AUC$experimentid)
load("/pfs/downloadgCSIProcessedSensitivity/gCSI_2017_sensitivity.RData")
old_experimentds <- which(sensitivityold$info$drugid %in% setdiff(rownames(druginfo), gCSI_AUC$drugid))
sensitivity_info <- rbind(gCSI_AUC[, c("cellid", "drugid", "experimentid")], cbind(sensitivityold$info[old_experimentds, c("cellid", "drugid")], "experimentid"=NA))
sensitivity_profiles <-  rbind(gCSI_AUC[,"auc_recomputed", drop=FALSE], sensitivityold$profiles[old_experimentds, "auc_recomputed", drop=FALSE])


gCSI_V1_2_metrics <- cbind(gCSI_v1_2_metric, "cellid"=NA, "drugid"=NA)
gCSI_V1_2_metrics$drugid <- curationDrug$unique.drugid[match(tolower(gsub(badchars, "", gCSI_V1_2_metrics$DrugName)), tolower(gsub(badchars, "", curationDrug$gCSI.drugid)))]
table(is.na(gCSI_V1_2_metrics$drugid))

gCSI_V1_2_metrics$cellid <- curationCell$unique.cellid[match(tolower(gsub(badchars, "", gCSI_V1_2_metrics$CellLineName)), tolower(gsub(badchars, "", curationCell$gCSI.cellid)))]
table(is.na(gCSI_V1_2_metrics$cellid))
rownames(gCSI_V1_2_metrics) <- sprintf("%s_%s_%s", gCSI_V1_2_metrics$cellid, gCSI_V1_2_metrics$drugid, gCSI_V1_2_metrics$ExperimentNumber)


sensitivity_profiles[rownames(gCSI_V1_2_metrics), "GR_AOC"] <- gCSI_V1_2_metrics$GR_AOC
sensitivity_profiles[rownames(gCSI_V1_2_metrics), "GReff"] <- gCSI_V1_2_metrics$GR_05uM_fit

load("/pfs/downloadgCSIProcessedRNA/gCSI_2017_molecprofile.RData")

gCSI@molecularProfiles$rnaseq$cellid <- gsub("U266B1", "U-266", gCSI@molecularProfiles$rnaseq$cellid) #compares molecular profile cell.id to unique.id of cell.info for checkPSet, so we need to change to U-266 (unique id of U266B1).
gCSI@molecularProfiles$mutation$cellid <- gsub("U266B1", "U-266", gCSI@molecularProfiles$mutation$cellid)
gCSI@molecularProfiles$mutation$cellid <- gsub("HEC1B", "HEC-1", gCSI@molecularProfiles$mutation$cellid)

gCSI_V1_2 <- PharmacoSet(molecularProfiles=gCSI@molecularProfiles,
                       name="gCSI",
                       cell=celline.info,
                       drug=druginfo,
                       sensitivityInfo=sensitivity_info,
                       sensitivityRaw=NULL,
                       sensitivityProfiles=sensitivity_profiles,
                       sensitivityN=NULL,
                       curationCell=curationCell,
                       curationDrug=curationDrug,
                       curationTissue=curationTissue,
                       datasetType="sensitivity")

save(gCSI_V1_2, file="/pfs/out/gCSI_2018.RData")
    
return (gCSI_V1_2)
    
    
}


library(PharmacoGxPrivate)
library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)

getgCSI(verbose=FALSE, nthread=1)
