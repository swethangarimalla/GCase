install.packages("openxlsx", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("ggpubr", repos = "http://cran.us.r-project.org")
install.packages("gridExtra", repos = "http://cran.us.r-project.org")

library(openxlsx)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(stats)
library(dplyr)
library(tidyr)
library(arsenal)
library(stringr)

#-------------------Functions----------------
roundUp <- function(x,to=10){to*(x%/%to + as.logical(x%%to))}
sem <- function(x) sd(x)/sqrt(length(x))
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

tukey_hsd <- function(x, ...){
  comparison <- comparison2 <- adj.p.value <- p.adj <-
    term <- group1 <- group2 <- NULL
  res <- TukeyHSD(x, ...) %>%
    broom::tidy() %>%
    mutate(comparison2 = comparison) %>%
    separate(comparison2, into= c("group2", "group1"), sep = "-") %>%
    rename(p.adj = adj.p.value) %>%
    mutate(p.adj = signif(p.adj, 2)) %>%
    select(term, group1, group2, everything())
  res
}
#---------------------------------------------

dir<-"/Users/swethagarimalla/Desktop/RCode/Bioanalytics/Experiments/PRV-2018-007/GCase/Raw Data/"
setwd(dir) #set directory as current working directory

BSAfiles <- list.files(pattern = '*BSA*') #Get all BSA files
GCfiles <- list.files(pattern = '*Gcase*')#Get all GC files 

BSA<-data.frame();
GC<-data.frame();

for(i in 1:length(BSAfiles)){ 
  print(i)
  #Load workbooks 
  bwkbk <- loadWorkbook(BSAfiles[i])
  gwkbk <- loadWorkbook(GCfiles[i])
  
  #find and read Dilution sheets
  bsheet <- grep("Dilution", names(bwkbk), value = TRUE)
  gsheet <- grep("Dilution", names(gwkbk), value = TRUE)
  
  bfile <- as.data.frame(read.xlsx(bwkbk,
                                   sheet = bsheet, 
                                   startRow = 13,
                                   colNames = T,
                                   skipEmptyRows = TRUE))
  gfile <- as.data.frame(read.xlsx(gwkbk,
                                   sheet = gsheet, 
                                   startRow = 13,
                                   colNames = T,
                                   skipEmptyRows = TRUE))
  
  #bind into one data frame
  BSA<-rbind(BSA, bfile)
  GC<-rbind(GC, gfile)
  
  
} 

colnames(BSA) <- paste0("B_", colnames(BSA)) #Add B_ and G_ to all the 
colnames(GC) <- paste0("G_", colnames(GC))   #colnames so we don't
#have duplicate colnames

BSA <- BSA[which(BSA$B_Well != " "),]        # Remove all empty rows
GC <- GC[which(GC$G_Well != " "),]


BSA$B_Abs <- as.numeric(BSA$B_Abs)          #make sure measurements are
GC$G_RFU<-as.numeric(GC$G_RFU)              #numeric

##Grubbs for BSA
BSA$B_Absmean <- ave(BSA$B_Abs, BSA$B_Sample, FUN = function(x) mean(x, na.rm=T)) #the mean including 0s
BSA$B_Abssd <- ave(BSA$B_Abs,BSA$B_Sample, FUN = function(x) sd(x, na.rm=T)) #the sd including 0s
BSA$B_Gvalue <- (BSA$B_Abs - BSA$B_Absmean)/BSA$B_Abssd #the G-value = (x-mean(x)/sd(x))

##if the G-value > 1.15 or Gvalue < - 1.15 for 3 replicates, it is an outlier and we do not consider the Quantity value
##for the downstream analysis so we assign it a value of NA
BSA$BSAAbs<-BSA$B_Abs
BSA[which(abs(BSA$B_Gvalue) >= 1.15), c("BSAAbs")] <- NA

GC$G_RFU[which((GC$G_RFU/GC$G_Dilution.factor) < 2.02)] <- 2.02 * GC$G_Dilution.factor[which((GC$G_RFU/GC$G_Dilution.factor) < 2.02)]
GC$G_RFU[which(is.nan(GC$G_RFU))] <- 2.02 * GC$G_Dilution.factor[which(is.nan(GC$G_RFU) < 2.02)]
GC$G_RFU[which(is.na(GC$G_RFU))] <- 2.02 * GC$G_Dilution.factor[which(is.nan(GC$G_RFU) < 2.02)]

##Grubbs for GC
GC$G_RFUmean <- ave(GC$G_RFU, GC$G_Sample, FUN = function(x) mean(x, na.rm=T)) #the mean including 0s
GC$G_RFUsd <- ave(GC$G_RFU,GC$G_Sample, FUN = function(x) sd(x, na.rm=T)) #the sd including 0s
GC$G_Gvalue <- (GC$G_RFU - GC$G_RFUmean)/GC$G_RFUsd #the G-value = (x-mean(x)/sd(x))

##if the G-value > 1.15 or Gvalue < - 1.15 for 3 replicates, it is an outlier and we do not consider the Quantity value
##for the downstream analysis so we assign it a value of NA
GC$GCRFU<-GC$G_RFU
GC[which(abs(GC$G_Gvalue) >= 1.15), c("GCRFU")] <- NA

GC2<-GC[order(GC$G_Sample),]
BSA2<-BSA[order(BSA$B_Sample),]
merged <- cbind(BSA2, GC2)

if(length(merged$B_Sample[!(merged$B_Sample %in% merged$G_Sample)]) == 0) { merged$Sample <- merged$B_Sample}
merged$unitsPerMG <- as.numeric(merged$GCRFU)/as.numeric(merged$BSAAbs)
merged$Sample <- merged$B_Sample
merged$mean_unitsPerMG <- ave(merged$unitsPerMG, merged$Sample, FUN = function(x) mean(x, na.rm=T))

merged$Sample[which(merged$Sample == "Empty")] <- "Empty-Empty"
merged$Mouse <-matrix(unlist(strsplit(merged$Sample, "-")), 
                      ncol = 2, 
                      byrow=TRUE)[,1]
merged$Tissue <-matrix(unlist(strsplit(merged$Sample, "-")), 
                       ncol = 2, 
                       byrow=TRUE)[,2]

merged <- merged[, c("Mouse", "Tissue", "Sample", "B_Abs", "BSAAbs", "B_Gvalue", 
                     "G_RFU", "GCRFU", "G_Gvalue", "unitsPerMG", "mean_unitsPerMG")]
tfile <- grep("Tissue Manifest", list.files(), value = TRUE)
tInventory <-loadWorkbook(tfile)
tInventory <- as.data.frame(read.xlsx(tInventory,
                                      sheet = grep("Inventory", getSheetNames(tfile), value = TRUE),
                                      colNames = T))
tInventory<- unique(tInventory[, c(grep("ID", colnames(tInventory)), 
                                   grep("Group", colnames(tInventory)),
                                   grep("Gender", colnames(tInventory)),
                                   grep("Tissue", colnames(tInventory)))])

colnames(tInventory) <- c("Sample.ID","PGID", "Group", "Gender", "TissueType")

tInventory$Sample.ID <- trim(tInventory$Sample.ID)
tInventory$TissueType <- trim(tInventory$TissueType)

tInventory$TissueType[which(tInventory$TissueType == "Rest of Cortex")] <- "Cortex"

tInventory$Sample <- paste0(tInventory$Sample.ID, "-", tInventory$TissueType)
merged<-unique(merge(merged, tInventory, 
                     by = "Sample",
                     all.x = TRUE))

merged <- merged[!(merged$Tissue == "Empty"),]

merged$Group<-as.character(merged$Group)
merged$Group[which(merged$Group == "1")] <- "PBS + Excipient" 
merged$Group[which(merged$Group == "2")] <- "CBE + Excipient"
merged$Group[which(merged$Group == "3")] <- "CBE + 3.2E10 vg PR001A"


#--------Plot BoxPlots by Tissue and Treatment---------
## Plot results
setwd("..")
outdir<-paste0("AutomatedOutput", Sys.time())
dir.create(outdir); #create dir/plotsByTissue to hold plots
setwd(outdir); #cd into dir/plotsByTissue


for (tissue in levels(factor(merged$Tissue))){
  byTissue<-merged[which(grepl(trim(tissue), merged$Tissue), arr.ind = T),]
  byTissue <- byTissue[which(!is.na(byTissue$mean_unitsPerMG)),]
  byTissue$Group<-as.character(byTissue$Group)
  
  byTissue<-unique(byTissue[,c("Group", "Tissue","mean_unitsPerMG")])
  
  gd_byTissue <- byTissue %>%
    group_by(Group) %>%
    summarise(mean = mean(mean_unitsPerMG), 
              sem = sem(mean_unitsPerMG))
  
  
  byTissue <- merge(byTissue, 
                    gd_byTissue, 
                    by = c("Group"), 
                    all.x = TRUE)
  
  byTissue$Group<-factor(byTissue$Group, levels=c("PBS + Excipient",
                                                  "CBE + Excipient", 
                                                  "CBE + 3.2E10 vg PR001A"))
  
  
  Combos <- aov(mean_unitsPerMG ~ Group, data = byTissue) %>% 
    tukey_hsd(... = "Group")
  Combos <- Combos[which(Combos$p.adj <= 0.05),]
  Combos <- rbind(Combos[which(Combos$group1 == "CBE + Excipient"),],
                  Combos[which(Combos$group2 == "CBE + Excipient"),])
  
  from <-3300#roundUp(max((byTissue$mean + byTissue$sem), na.rm = T))
  num <- dim(Combos)[1] 
  y.position <-seq(from = from, by = from/10, length.out = num)
  
  ggplot(byTissue, aes(x=Group, y=mean)) + #, fill=Group
    geom_bar(stat="identity", aes(fill=Group), 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean - sem, ymax=mean + sem), width=.2,
                  position=position_dodge(.9)) +
    labs(title=paste0("PRV-2018-007: GCase Activity Enzyme - ", tissue), x="Dose Group", y = "Effective Enzyme Activity (units/mg)")+
    theme_classic() +
    scale_fill_manual(values=c("grey", "red", 
                               "skyblue2", "skyblue3", "skyblue4", "darkblue")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_pvalue_manual( data = Combos, label = "p.adj",
                        y.position = y.position) +
    stat_compare_means(method = "anova",
                       label.x.npc = "right", 
                       label.y.npc = "top") + 
    theme(axis.text=element_text(size = 14,face="bold"),
          axis.title=element_text(size = 14,face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position="none") +
    #scale_y_continuous(expand=c(0,0), limits = c(0, max(merged$mean + 4*merged$sem))) +
    scale_x_discrete(labels = c( bquote("PBS + Excipient"), 
                                 bquote("CBE + Excipient"),
                                 bquote("CBE + 3.2 x " * 10^10 * " vg PR001A")))  
  
  ggsave(paste0("PRV_2018_007_GCase_", tissue, ".jpeg"))
  
}

gd_merged <- merged %>%
  group_by(Group, Tissue) %>%
  summarise(mean = mean(mean_unitsPerMG), 
            sem = sem(mean_unitsPerMG))
merged <- merge(merged, 
                gd_merged, 
                by = c("Group", "Tissue"), 
                all.x = TRUE)

## Create a new workbook
wkbk <- createWorkbook()
addWorksheet(wb = wkbk, "Rscript Output Data Table", gridLines = TRUE) #add a worksheet to workbook. 
writeDataTable(wkbk, "Rscript Output Data Table", merged)

unique_merged<-unique(merged[,c("Mouse", "Group","mean_unitsPerMG")])
addWorksheet(wb = wkbk, "Mean Units Per Sample", gridLines = TRUE) #add a worksheet to workbook. 
writeDataTable(wkbk, "Mean Units Per Sample", unique_merged)
## Save workbook to working directory
saveWorkbook(wkbk, file = paste0("RscriptEnzymeAssayOutput", Sys.Date(), ".xlsx"), overwrite = TRUE)


install.packages("devtools")
library(devtools)
install_github("Yue-Jiang/pzfx")
library(pzfx)

class(unique_merged$Mouse) <- "numeric"
attach(unique_merged)
unique_merged<-unique_merged[order(Group),] 
write_pzfx(unique_merged, "MeanUnitsPerMG_007.pzfx", row_names=TRUE)
out_df <- read_pzfx("MeanUnitsPerMG_007.pzfx", table=1)
head(out_df)

rFile <- "/Users/swethagarimalla/Desktop/RCode/Bioanalytics/Experiments/Scripts/GCase/GCase_v6_03142019.R"
newFile <- str_split(rFile, pattern = "/")
newFile <- unlist(newFile)[length(unlist(newFile))]
newFile <- paste0(getwd(), "/", newFile, "-", Sys.time(),".R")
file.copy(from = rFile, to = newFile)





