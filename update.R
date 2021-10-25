

rm(list=ls()) 
library(ggplot2)
library(vegan)
library(reshape2)
library(Hmisc)  
library(plotrix)
library(phyloseq)
library(bioDist) 
library(igraph)
library(edgeR)
library(formatR)
library(corrplot)
library(statnet)
library(circlize)
library(plyr)
library(agricolae)
library(lattice)
library(latticeExtra)
library(ggrepel)
library(dplyr)
library(ggcor)
        
##### Import Data #####
otu <- read.table("otutab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
otu <- as.matrix(otu)

design <- read.table("metadata.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
design$Concentration<-factor(design$Concentration,c("0","75","600"))
design$Site<-factor(design$Site,c("root","soil","film"))
design$Stage<-factor(design$Stage,c("Flowering","Maturity"))
design$Treatment1 <- factor(design$Treatment1,c("0-root","0-soil","75-soil","75-film","75-root","600-film","600-root","600-soil"))

tax <- read.table("taxonomy.txt", row.names=1, sep="\t", header=T ,stringsAsFactors=F,quote="")
tax$labels <- tax$Phylum
tax[ rownames(tax)[tax$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
tax[ rownames(tax)[tax$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
tax[ rownames(tax)[tax$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
tax[ rownames(tax)[tax$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(tax$labels)

##### Define sample types #####
floweringsamples <- rownames(design)[which(design$Stage == "Flowering")]
maturitysamples <- rownames(design)[which(design$Stage == "Maturity")]

edgeR_flowering <- DGEList(counts=flowering_otu, 
                           group=design_flowering$Concentration,
                           genes=tax_flowering)

edgeR_flowering <- calcNormFactors(edgeR_flowering)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_flowering <- cpm(edgeR_flowering, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_flowering <- otu_table(norm_flowering,taxa_are_rows=T)
phy_tax_flowering <-tax_table(as.matrix(tax_flowering))
phy_design_flowering <- sample_data(design_flowering)
physeq_norm_flowering <- phyloseq(phy_flowering,phy_tax_flowering,phy_design_flowering)
sample_data(physeq_norm_flowering)$Treatment1<- factor(sample_data(physeq_norm_flowering)$Treatment1,levels=c("0-root","0-soil","75-soil","75-film","75-root","600-film","600-root","600-soil"))

edgeR_Maturity <- DGEList(counts=Maturity_otu, 
                           group=design_Maturity$Concentration,
                           genes=tax_Maturity)

edgeR_Maturity <- calcNormFactors(edgeR_Maturity)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_Maturity <- cpm(edgeR_Maturity, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_Maturity <- otu_table(norm_Maturity,taxa_are_rows=T)
phy_tax_Maturity <-tax_table(as.matrix(tax_Maturity))
phy_design_Maturity <- sample_data(design_Maturity)
physeq_norm_Maturity <- phyloseq(phy_Maturity,phy_tax_Maturity,phy_design_Maturity)
sample_data(physeq_norm_Maturity)$Treatment1<- factor(sample_data(physeq_norm_Maturity)$Treatment1,levels=c("0-root","0-soil","75-soil","75-film","75-root","600-film","600-root","600-soil"))


## Create bray-curtis dissimiliartiy matrix
all_dis_Maturity <- vegdist(t(otu_table(physeq_norm_Maturity)),method="bray")
        
##### PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_Maturity <- adonis(all_dis_Maturity ~Site+Concentration+Treatment1*Treatment2, data=design_Maturity, permutations=9999)

paov_all_Maturity <- adonis(all_dis_Maturity ~Site+Concentration+Treatment1*Treatment2, data=design_Maturity, permutations=9999)
#####        

##### Beta diversity #####
pcoa_norm_flowering <- ordinate(physeq_norm_flowering,"PCoA","bray")
pcoa_all_flowering <- plot_ordination(physeq_norm_flowering, pcoa_norm_flowering, type="sites", color="Site", shape="Concentration")
pcoa_all_flowering <- pcoa_all_flowering+
          geom_point(size=4)+
          xlab(paste("PCo 1", paste("(",round(pcoa_norm_flowering$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
          ylab(paste("PCo 2", paste("(",round(pcoa_norm_flowering$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
          theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
          guides(color=guide_legend(nrow=1,byrow=TRUE))+
          guides(shape=guide_legend(nrow=1,byrow=TRUE))+
          theme(plot.title = element_text(face="bold", hjust = 0.5))+
          ggtitle("Flowering Community")

pcoa_norm_Maturity <- ordinate(physeq_norm_Maturity,"PCoA","bray")
pcoa_all_Maturity <- plot_ordination(physeq_norm_Maturity, pcoa_norm_Maturity, type="sites", color="Site", shape="Concentration")
pcoa_all_Maturity <- pcoa_all_Maturity+
  geom_point(size=4)+
  xlab(paste("PCo 1", paste("(",round(pcoa_norm_Maturity$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_norm_Maturity$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Maturity Community")
#####

## aov and LSD
PHYLUM_mat_flowering_aov <- read.csv("treatments.csv", header=T, row.names=1, stringsAsFactors=F, na.strings="NA",check.names=F)
PHYLUM_mat_flowering_aov$Treatment1 <-as.factor(PHYLUM_mat_flowering_aov$Treatment) 

for(i in PHYLUM_mat_flowering_aov[,4:13]) {
  fit <- aov(i~Treatment,data = PHYLUM_mat_flowering_aov)
  print(summary(fit))
  out <- LSD.test(fit,"Treatment")
  print(out$groups)
}

### Defining OTU colors by phylum (using the taxonomy file)
tax_flowering$cols <- tax_flowering$labels
table(tax_flowering$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_flowering, 1, mean) < 0.1)
low_count_phyla_flowering <- rownames(PHYLUM_mat_flowering)[sort(apply(PHYLUM_mat_flowering, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_flowering){
  tax_flowering[ rownames(tax_flowering)[tax_flowering$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_flowering$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_flowering <- rownames(PHYLUM_mat_flowering)[sort(apply(PHYLUM_mat_flowering, 1, mean), decr=T) > 1]
abundant_phyla_flowering
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Alphaproteobacteria" ], ]$cols <- "lightcoral"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Betaproteobacteria" ], ]$cols <- "palegreen3"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Gammaproteobacteria" ], ]$cols <- "palegreen1"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Deltaproteobacteria" ], ]$cols <- "palegreen4"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Actinobacteria" ], ]$cols <- "dodgerblue"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Firmicutes" ], ]$cols <- "plum1"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Acidobacteria" ], ]$cols <- "palegreen2"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Chloroflexi" ], ]$cols <- "steelblue1"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Candidatus_Saccharibacteria" ], ]$cols <- "steelblue4"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Verrucomicrobia" ], ]$cols <- "palevioletred3"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Gemmatimonadetes" ], ]$cols <- "peachpuff2"
tax_flowering[ rownames(tax_flowering)[tax_flowering$labels=="Nitrospirae" ], ]$cols <- "orange"

## collaps OTU colors to prepare Phylum level colors
label_cols_flowering <- tax_flowering[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_flowering <- ddply(label_cols_flowering, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_flowering) <- PHYLA_label_cols_flowering[,1]
PHYLA_label_cols_flowering <- PHYLA_label_cols_flowering[c(abundant_phyla_flowering, low_count_phyla_flowering),]
PHYLA_label_cols_flowering

## Legend for Phylum colors
PHYLA_label_cols_flowering_legend <- PHYLA_label_cols_flowering[1:14,]
PHYLA_label_cols_flowering_legend[14,1] <- "other"
rownames(PHYLA_label_cols_flowering_legend)[14] <- "other"
PHYLA_label_cols_flowering_legend
#####


##### Plot 
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylumar_flowering <- barplot(as.matrix(PHYLUM_mat_flowering_mean_treat), col=PHYLA_label_cols_flowering[rownames(PHYLUM_mat_flowering_mean_treat),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylumar_flowering, labels=colnames(PHYLUM_mat_flowering_mean_treat), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Flowering Community")
legend(9.5, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_flowering_legend$labels), 
       fill=rev(PHYLA_label_cols_flowering_legend$cols), 
       border=rev(PHYLA_label_cols_flowering_legend$cols) )

par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylumar_Maturity <- barplot(as.matrix(PHYLUM_mat_Maturity_mean_treat), col=PHYLA_label_cols_Maturity[rownames(PHYLUM_mat_Maturity_mean_treat),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylumar_Maturity, labels=colnames(PHYLUM_mat_Maturity_mean_treat), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Maturity Community")
legend(9.5, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_Maturity_legend$labels), 
       fill=rev(PHYLA_label_cols_Maturity_legend$cols), 
       border=rev(PHYLA_label_cols_Maturity_legend$cols) )
#####



## Alpha Diversity Index
cp_flowering_otu <- as.data.frame(t(rrarefy(t(flowering_otu), min(colSums(flowering_otu)))))
cp_flowering_otu_shannon <- vegan::diversity(t(cp_flowering_otu), index = "shannon")
design_flowering$Shannon <- cp_flowering_otu_shannon
cp_flowering_otu_simpson <- vegan::diversity(t(cp_flowering_otu), index = "simpson")
design_flowering$Simpson <- cp_flowering_otu_simpson
chao1 <- t(estimateR(t(cp_flowering_otu)))
chao1 <- as.data.frame(chao1)
design_flowering$chao1 <- chao1$chao1
design_flowering$ACE <- chao1$ACE

cp_maturity_otu <- as.data.frame(t(rrarefy(t(maturity_otu), min(colSums(maturity_otu)))))
cp_maturity_otu_shannon <- vegan::diversity(t(cp_maturity_otu), index = "shannon")
design_maturity$Shannon <- cp_maturity_otu_shannon
cp_maturity_otu_simpson <- vegan::diversity(t(cp_maturity_otu), index = "simpson")
design_maturity$simpson <- cp_maturity_otu_simpsoner
chao1 <- t(estimateR(t(cp_maturity_otu)))
chao1 <- as.data.frame(chao1)
design_maturity$chao1 <- chao1$S.chao1
design_maturity$ACE <- chao1$S.ACE

for(i in design_flowering[,9:12]) {
  fit <- aov(i~Treatment1,data = design_flowering)
  print(summary(fit))
  out <- LSD.test(fit,"Treatment1", p.adj="none" )
  print(out$groups)
}

shannon_maturity <- ggplot(design_maturity, aes(x=factor(Site,levels = c("Bulk soil","Plastisphere compartment")), y=Shannon, fill =Concentration))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) + 
  scale_fill_manual(values=c("#ffffff","#5ba9de","b55e6e"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Site", y="Maturity shannon index")

#####


##### Defining  Bulk soil/Rhizosphere compartment/Plastisphere compartment bacteria communities for beta diversity analysis #####
## Apply sequence count threshold to bacteria soil community and TMM normalize counts
##### Define sample types #####
floweringsamples_soil <- rownames(design_flowering)[which(design_flowering$Site == "Bulk soil")]
floweringsamples_root <- rownames(design_flowering)[which(design_flowering$Site == "Rhizosphere compartment")]
floweringsamples_film <- rownames(design_flowering)[which(design_flowering$Site == "Plastisphere compartment")]

flowering_soil <- flowering_otu[, floweringsamples_soil ]
flowering_soil <- flowering_soil[rowSums(flowering_soil) > 0,]
dim(flowering_soil)
keep_flowering_soil <- which(rowSums(flowering_soil >= 2) >= 3)

flowering_soil <- flowering_soil[keep_flowering_soil,]

dim(flowering_soil)

flowering_soil_RA <- t(t(flowering_soil)/colSums(flowering_soil)) * 100
colSums(flowering_soil_RA)
nrow(flowering_soil_RA)

tax_soil_flowering <- tax_flowering[rownames(flowering_soil),]
design_flowering_soil <- droplevels(design_flowering[floweringsamples_soil,])

edgeR_flowering_soil <- DGEList(counts=flowering_soil, 
                                group=design_flowering_soil$Concentration,
                                genes=tax_soil_flowering)

edgeR_flowering_soil <- calcNormFactors(edgeR_flowering_soil)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_soil_flowering <- cpm(edgeR_flowering_soil, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_soil_flowering <- vegdist(t(otu_table(physeq_soil_norm_flowering)),method="bray")

##### B overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_soil_flowering <- adonis(all_dis_soil_flowering ~Concentration+Treatment1*Treatment3, data=design_flowering_soil, permutations=9999)

##### Supplementary Table S2: global PERMANOVA #####
paov_all_soil_flowering  


flowering_root <- flowering_otu[, floweringsamples_root ]
flowering_root <- flowering_root[rowSums(flowering_root) > 0,]
dim(flowering_root)
keep_flowering_root <- which(rowSums(flowering_root >= 2) >= 3)

flowering_root <- flowering_root[keep_flowering_root,]

dim(flowering_root)

flowering_root_RA <- t(t(flowering_root)/colSums(flowering_root)) * 100
colSums(flowering_root_RA)
nrow(flowering_root_RA)

tax_root_flowering <- tax_flowering[rownames(flowering_root),]
design_flowering_root <- droplevels(design_flowering[floweringsamples_root,])

edgeR_flowering_root <- DGEList(counts=flowering_root, 
                                group=design_flowering_root$Concentration,
                                genes=tax_root_flowering)

edgeR_flowering_root <- calcNormFactors(edgeR_flowering_root)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_root_flowering <- cpm(edgeR_flowering_root, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_root_flowering <- vegdist(t(otu_table(physeq_root_norm_flowering)),method="bray")

##### B overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_root_flowering <- adonis(all_dis_root_flowering ~Concentration+Treatment1*Treatment3, data=design_flowering_root, permutations=9999)

##### Supplementary Table S2: global PERMANOVA #####
paov_all_root_flowering  


flowering_film <-otu[, floweringsamples_film]

flowering_film <- flowering_film[rowSums(flowering_film) > 0,]
dim(flowering_film)
keep_flowering_film <- which(rowSums(flowering_film >= 2) >= 3)

flowering_film <- flowering_film[keep_flowering_film,]

dim(flowering_film)

flowering_film_RA <- t(t(flowering_film)/colSums(flowering_film)) * 100
colSums(flowering_film_RA)
nrow(flowering_film_RA)

tax_film_flowering <- tax_flowering[rownames(flowering_film),]
design_flowering_film <- droplevels(design_flowering[floweringsamples_film,])
design_flowering_film <- design[floweringsamples_film,]

edgeR_flowering_film <- DGEList(counts=flowering_film, 
                                group=design_flowering_film$Concentration,
                                genes=tax_film_flowering)

edgeR_flowering_film <- calcNormFactors(edgeR_flowering_film)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_film_flowering <- cpm(edgeR_flowering_film, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_film_flowering <- vegdist(t(otu_table(physeq_film_norm_flowering)),method="bray")

##### B overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_film_flowering <- adonis(all_dis_film_flowering ~Concentration, data=design_flowering_film, permutations=9999)

##### Supplementary Table S2: global PERMANOVA #####
paov_all_film_flowering  


## Apply sequence count threshold to bacteria soil community and TMM normalize counts
##### Define sample types #####
floweringsamples_ck <- rownames(design_flowering)[which(design_flowering$Concentration == "0")]
floweringsamples_low <- rownames(design_flowering)[which(design_flowering$Concentration == "75")]
floweringsamples_high <- rownames(design_flowering)[which(design_flowering$Concentration == "600")]

flowering_ck <- otu[, c(floweringsamples_ck,maturitysamples_soil_low)]
flowering_ck <- flowering_ck[rowSums(flowering_ck) > 0,]
dim(flowering_ck)
keep_flowering_ck <- which(rowSums(flowering_ck >= 2) >= 3)

flowering_ck <- flowering_ck[keep_flowering_ck,]

dim(flowering_ck)

flowering_ck_RA <- t(t(flowering_ck)/colSums(flowering_ck)) * 100
colSums(flowering_ck_RA)
nrow(flowering_ck_RA)

tax_ck_flowering <- tax_flowering[rownames(flowering_ck),]
design_flowering_ck <- droplevels(design_flowering[floweringsamples_ck,])
design_flowering_ck <- design[c(floweringsamples_ck,maturitysamples_soil_low),]

edgeR_flowering_ck <- DGEList(counts=flowering_ck, 
                              group=design_flowering_ck$Site,
                              genes=tax_ck_flowering)

edgeR_flowering_ck <- calcNormFactors(edgeR_flowering_ck)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_ck_flowering <- cpm(edgeR_flowering_ck, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_ck_flowering <- vegdist(t(otu_table(physeq_ck_norm_flowering)),method="bray")

flowering_low <- flowering_otu[, floweringsamples_low ]
flowering_low <- flowering_low[rowSums(flowering_low) > 0,]
dim(flowering_low)
keep_flowering_low <- which(rowSums(flowering_low >= 2) >= 3)

flowering_low <- flowering_low[keep_flowering_low,]

dim(flowering_low)

flowering_low_RA <- t(t(flowering_low)/colSums(flowering_low)) * 100
colSums(flowering_low_RA)
nrow(flowering_low_RA)

tax_low_flowering <- tax_flowering[rownames(flowering_low),]
design_flowering_low <- droplevels(design_flowering[floweringsamples_low,])

edgeR_flowering_low <- DGEList(counts=flowering_low, 
                               group=design_flowering_low$Concentration,
                               genes=tax_low_flowering)

edgeR_flowering_low <- calcNormFactors(edgeR_flowering_low)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_low_flowering <- cpm(edgeR_flowering_low, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_low_flowering <- vegdist(t(otu_table(physeq_low_norm_flowering)),method="bray")


flowering_high <- flowering_otu[, floweringsamples_high ]
flowering_high <- flowering_high[rowSums(flowering_high) > 0,]
dim(flowering_high)
keep_flowering_high <- which(rowSums(flowering_high >= 2) >= 3)

flowering_high <- flowering_high[keep_flowering_high,]

dim(flowering_high)

flowering_high_RA <- t(t(flowering_high)/colSums(flowering_high)) * 100
colSums(flowering_high_RA)
nrow(flowering_high_RA)

tax_high_flowering <- tax_flowering[rownames(flowering_high),]
design_flowering_high <- droplevels(design_flowering[floweringsamples_high,])

edgeR_flowering_high <- DGEList(counts=flowering_high, 
                                group=design_flowering_high$Concentration,
                                genes=tax_high_flowering)

edgeR_flowering_high <- calcNormFactors(edgeR_flowering_high)

## Get TMM normalized counts expressed as relative abundance counts per million 
norm_high_flowering <- cpm(edgeR_flowering_high, normalized.lib.sizes=T, log=F)

## Create bray-curtis dissimiliartiy matrix
all_dis_high_flowering <- vegdist(t(otu_table(physeq_high_norm_flowering)),method="bray")
#####     


## Define Soil, Root and Film enriched BACTERIA
edgeR_soil_film_enrich <- DGEList(counts=flowering_soil_film, 
                                  group=design_soil_film$Site)

edgeR_soil_film_enrich <- calcNormFactors(edgeR_soil_film_enrich)

otu_norm_soil_film_enrich <- cpm(edgeR_soil_film_enrich, normalized.lib.sizes=T, log=F)

model_mat_enrich_soil_film <- model.matrix(~Block+Site, data=design_soil_film)

dge_enrich_soil_film <- estimateGLMRobustDisp(edgeR_soil_film_enrich, design=model_mat_enrich_soil_film)

fit_enrich_soil_film <- glmFit(dge_enrich_soil_film, design=model_mat_enrich_soil_film)
lrt_enrich_soil_film <- glmLRT(fit_enrich_soil_film, coef="Sitefilm")
tt_enrich_soil_film <- topTags(lrt_enrich_soil_film, n=Inf, p.value=1)
head(tt_enrich_soil_film$table)

film_enrich_soil <- tt_enrich_soil_film$table[tt_enrich_soil_film$table$logFC < 0 & tt_enrich_soil_film$table$FDR < 0.05,]
soil_enrich_soil <- tt_enrich_soil_film$table[tt_enrich_soil_film$table$logFC > 0 & tt_enrich_soil_film$table$FDR < 0.05,]

### Soil&Film
tt_enrich_soil_film <- as.data.frame(tt_enrich_soil_film)
forMA_soil_film <- data.frame(2^tt_enrich_soil_film$logCPM,
                              tt_enrich_soil_film$logFC,
                              tt_enrich_soil_film$FDR<0.05)
colnames(forMA_soil_film) <- c("CPM", "logFC", "signif")

## define colors
forMA_soil_film$col[forMA_soil_film$signif==F] <- "dimgrey"
forMA_soil_film$col[forMA_soil_film$logFC>0 & forMA_soil_film$signif==T] <- "sienna4"
forMA_soil_film$col[forMA_soil_film$logFC<0 & forMA_soil_film$signif==T] <- "#DAA520"
# order for plotting colors
forMA_soil_film$ord[forMA_soil_film$col=="dimgrey"] <- 1
forMA_soil_film$ord[forMA_soil_film$col=="#DAA520"] <- 2
forMA_soil_film$ord[forMA_soil_film$col=="sienna4"] <- 3
forMA_soil_film <- forMA_soil_film[sort(forMA_soil_film$ord,ind=T,decr=F)$ix,]

## define pch
forMA_soil_film$pch[forMA_soil_film$signif==F] <- 1
forMA_soil_film$pch[forMA_soil_film$signif==T] <- 16
        
plot(forMA_soil_film$CPM, forMA_soil_film$logFC, log="x", ylim=c(-8.2,7.6),
             main="Soil&Film",
             ylab="log fold change", xlab="average abundance CPM",
             col=forMA_soil_film$col, cex=1, pch=forMA_soil_film$pch, las=1)


edgeR_root_film_enrich <- DGEList(counts=flowering_root_film, 
                                  group=design_root_film$Site)

edgeR_root_film_enrich <- calcNormFactors(edgeR_root_film_enrich)

otu_norm_root_film_enrich <- cpm(edgeR_root_film_enrich, normalized.lib.sizes=T, log=F)

model_mat_enrich_root_film <- model.matrix(~Block+Site, data=design_root_film)

dge_enrich_root_film <- estimateGLMRobustDisp(edgeR_root_film_enrich, design=model_mat_enrich_root_film)

fit_enrich_root_film <- glmFit(dge_enrich_root_film, design=model_mat_enrich_root_film)
lrt_enrich_root_film <- glmLRT(fit_enrich_root_film, coef="Sitefilm")
tt_enrich_root_film <- topTags(lrt_enrich_root_film, n=Inf, p.value=1)
head(tt_enrich_root_film$table)

film_enrich_root_root <- tt_enrich_root_film$table[tt_enrich_root_film$table$logFC < 0 & tt_enrich_root_film$table$FDR < 0.05,]
root_enrich_root_root <- tt_enrich_root_film$table[tt_enrich_root_film$table$logFC > 0 & tt_enrich_root_film$table$FDR < 0.05,]

### Soil&Film
tt_enrich_root_film <- as.data.frame(tt_enrich_root_film)
forMA_root_film <- data.frame(2^tt_enrich_root_film$logCPM,
                              tt_enrich_root_film$logFC,
                              tt_enrich_root_film$FDR<0.05)
colnames(forMA_root_film) <- c("CPM", "logFC", "signif")

## define colors
forMA_root_film$col[forMA_root_film$signif==F] <- "dimgrey"
forMA_root_film$col[forMA_root_film$logFC>0 & forMA_root_film$signif==T] <- "sienna4"
forMA_root_film$col[forMA_root_film$logFC<0 & forMA_root_film$signif==T] <- "forestgreen"
# order for plotting colors
forMA_root_film$ord[forMA_root_film$col=="dimgrey"] <- 1
forMA_root_film$ord[forMA_root_film$col=="forestgreen"] <- 2
forMA_root_film$ord[forMA_root_film$col=="sienna4"] <- 3
forMA_root_film <- forMA_root_film[sort(forMA_root_film$ord,ind=T,decr=F)$ix,]

## define pch
forMA_root_film$pch[forMA_root_film$signif==F] <- 1
forMA_root_film$pch[forMA_root_film$signif==T] <- 16

plot(forMA_root_film$CPM, forMA_root_film$logFC, log="x", ylim=c(-7.6,7.6),
     main="Root&Film",
     ylab="log fold change", xlab="average abundance CPM",
     col=forMA_root_film$col, cex=1, pch=forMA_root_film$pch, las=1)
text(rep(150,150), c(-7,7), label=c("root", "film"), adj=0, cex=1.5, col=c("forestgreen","sienna4"))
dev.off()
#####        



##### Identifiying OTUs with indicator species analysis #####
## Define indicator species for soil bacteria community.
set.seed(8046)

indicatorsp_soil_flowering <- multipatt(indic_root_flowering,indic_root_groups_flowering,func = "r.g",control=how(nperm=9999))
CK_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.0 == 1 & indic_soil_df_flowering$p.value < 0.05),])
LOW_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.75 == 1 & indic_soil_df_flowering$p.value < 0.05),])
HIGH_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.600 == 1 & indic_soil_df_flowering$p.value < 0.05),])

soil_r_values_flowering <- rbind(CK_soil_flowering,LOW_soil_flowering,HIGH_soil_flowering)
colnames(soil_r_values_flowering)[1:3] <-c("CK","LOW","HIGH")

## Range of correlation coefficients
range(soil_r_values_flowering[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(soil_r_values_flowering)))

## Proportion of soil bacteria OTUs responding to cropping system
length(unique(rownames(soil_r_values_flowering)))/nrow(norm_soil_flowering)

## Proportion of soil bacteria sequences responding to cropping system
soil_flowering_ra <- t(t(flowering_soil)/colSums(flowering_soil)) * 100
sum(colSums(soil_flowering_ra[unique(rownames(soil_r_values_flowering)),]))/sum(colSums(soil_flowering_ra))

## Management sensitive bulk soil bacteria OTus
indic_flowering_soil <- unique(rownames(soil_r_values_flowering))
sum(colSums(soil_flowering_ra[indic_flowering_soil,]))/sum(colSums(soil_flowering_ra))*100

## Define soil samples by management system
CK_samples_soil <- rownames(design_flowering_soil[design_flowering_soil$Concentration=="0",])
LOW_samples_soil <- rownames(design_flowering_soil[design_flowering_soil$Concentration=="75",])
HIGH_samples_soil <- rownames(design_flowering_soil[design_flowering_soil$Concentration=="600",])

## Soil bacteria: Calculate percentage of sequences classified into each phylum
flowering_soil_mso <- flowering_soil[unique(rownames(soil_r_values_flowering)),]
flowering_soil_mso_ra <- t(t(flowering_soil_mso)/colSums(flowering_soil_mso)) * 100

PHYLAnames_soil_mso_flowering <- names(sort(table(tax_soil_flowering[indic_flowering_soil,]$Phylum),decr=T))
length(PHYLAnames_soil_mso_flowering)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(flowering_soil_mso_ra)
for (i in PHYLAnames_soil_mso_flowering){
  x <- array(colSums(flowering_soil_mso_ra[rownames(tax_soil_flowering[indic_flowering_soil,])[which(tax_soil_flowering[indic_flowering_soil,]$Phylum == paste(i))], ,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_soil_mso_flowering)
colnames(y) <- paste(colnames(flowering_soil_mso_ra))
PHYLUM_mat_soil_flowering_mso <- y
PHYLUM_mat_soil_flowering_mso_mean <- sort(apply(PHYLUM_mat_soil_flowering_mso,1,mean),decr=T)
PHYLUM_mat_soil_flowering_mso <- PHYLUM_mat_soil_flowering_mso[names(PHYLUM_mat_soil_flowering_mso_mean),]

## Get sequences abundances by phylum
PHYLUM_mat_soil_flowering_mso_mean

## Phylum matrix for heatmap of phyla abundances across 
soil_flowering_mso_ra <- norm_soil_flowering[unique(rownames(soil_r_values_flowering)),]

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(soil_flowering_mso_ra)
for (i in PHYLAnames_soil_mso_flowering){
  x <- array(colSums(soil_flowering_mso_ra[rownames(tax_soil_flowering[indic_flowering_soil,])[which(tax_soil_flowering[indic_flowering_soil,]$Phylum == paste(i))], ,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_soil_mso_flowering)
colnames(y) <- paste(colnames(soil_flowering_mso_ra))
PHYLUM_mat_soil_flowering_mso <- y
PHYLUM_mat_soil_flowering_mso_mean <- sort(apply(PHYLUM_mat_soil_flowering_mso,1,mean),decr=T)
PHYLUM_mat_soil_flowering_mso <- PHYLUM_mat_soil_flowering_mso[names(PHYLUM_mat_soil_flowering_mso_mean),]

## Make matrix of phyla abundances by cropping system
soil_flowering_mso_phylum <-as.matrix(cbind(`CK`=apply(PHYLUM_mat_soil_flowering_mso[,CK_samples_soil],1,mean),
                                            `LOW`=apply(PHYLUM_mat_soil_flowering_mso[,LOW_samples_soil],1,mean),
                                            `HIGH`=apply(PHYLUM_mat_soil_flowering_mso[,HIGH_samples_soil],1,mean)))

## Identify indicator species in bulk soil bacteria communities
indic_root_flowering <- as.data.frame(t(norm_root_flowering))
indic_root_groups_flowering <- design_flowering_root$Concentration
length(unique(indic_root_groups_flowering))

pdf(paste0("soil_flowering_phylum_heartmap.pdf"),height=10, width=10)
par(las=1)
colors <- maPalette(l="lightgray",m="lightskyblue",h="orangered3", k=60)
heatmap.2(log2(soil_flowering_mso_phylum+1), col=colors, Rowv=F, Colv=F, margins=c(5,12),scale="none", trace="none"
          ,density.info="none",key.title=NA,srtCol=45,
          labRow=rownames(soil_flowering_mso_phylum),
          dendrogram="none",key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))
dev.off()

## flowering Soil management sensitive OTU abundances heatmap
soil_flowering_mso <-as.matrix(cbind(`CK`=apply(norm_soil_flowering[indic_flowering_soil,CK_samples_soil],1,mean),
                                     `LOW`=apply(norm_soil_flowering[indic_flowering_soil,LOW_samples_soil],1,mean),
                                     `HIGH`=apply(norm_soil_flowering[indic_flowering_soil,HIGH_samples_soil],1,mean)))

table(tax_soil_flowering[indic_flowering_soil,]$Phylum)
soil_flowering_sidecol <- tax_soil_flowering[indic_flowering_soil,]$Phylum
names(soil_flowering_sidecol) <- soil_flowering_sidecol

soil_flowering_sidecol[soil_flowering_sidecol == "Proteobacteria"] <- "seashell4"
soil_flowering_sidecol[soil_flowering_sidecol == "Actinobacteria"] <- "red"
soil_flowering_sidecol[soil_flowering_sidecol == "Bacteroidetes"] <- "darkslateblue"
soil_flowering_sidecol[soil_flowering_sidecol == "Candidatus_Saccharibacteria"] <- "tan4"
soil_flowering_sidecol[soil_flowering_sidecol == "Acidobacteria"] <- "blue"
soil_flowering_sidecol[soil_flowering_sidecol == "Chloroflexi"] <- "violetred4"
soil_flowering_sidecol[soil_flowering_sidecol == "Deinococcus-Thermus"] <- "orchid"
soil_flowering_sidecol[soil_flowering_sidecol == "Gemmatimonadetes"] <- "chartreuse"
soil_flowering_sidecol[soil_flowering_sidecol == "Firmicutes"] <- "orange"
soil_flowering_sidecol[soil_flowering_sidecol == "Verrucomicrobia"] <- "green"
soil_flowering_sidecol[soil_flowering_sidecol == "BRC1"] <- "yellow"

pdf(paste0("soil_OTU_flowering_phylum_heartmap.pdf"),height=10, width=10)
par(las=1)
colors <- maPalette(l="lightgray",m="lightskyblue",h="orangered3", k=60)
heatmap.2(log2(soil_flowering_mso+1), col=colors, Rowv=F, Colv=F, margins=c(5,12), scale="none", trace="none",
          density.info="none", key.title=NA, srtCol=45,
          labRow=paste(rownames(soil_flowering_mso),tax_soil_flowering[indic_flowering_soil,]$Family), RowSideColors=soil_flowering_sidecol,cexRow = 0.3,
          dendrogram="none",key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))
dev.off()

soilflowering_otu <-c("OTU2383","OTU1932","OTU834","OTU1547","OTU595","OTU86","OTU241","OTU97")

pdf(paste0("soilflowering_otuoxplot.pdf"),encoding="MacRoman",height=10,width=10)
par(mfrow=c(3,3))
for(i in soilflowering_otu)
  plotOTU(sample_table = design_flowering_soil, count_matrix = norm_soil_flowering,
          OTU_id = i, sample_factor = "Concentration",
          main = paste(i, tax_soil_flowering[i,"Family"]), ylab = "CPM",
          type = "box")
dev.off()
#####

  
        
## NETWORK
## soil
soil_flowering_occor <- rcorr(t(norm_soil_flowering),type=c("spearman"))
        
## Create data frame of co-occurring OTUs
soil_flowering_cor_df <- CorrDF(soil_flowering_occor$r,soil_flowering_occor$P)
soil_flowering_cor_df$padj <- p.adjust(soil_flowering_cor_df$p, method = "none")
        
## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
soil_flowering_cor_df_padj <- soil_flowering_cor_df[which(soil_flowering_cor_df$cor > 0.7),]
soil_flowering_cor_df_padj <- soil_flowering_cor_df_padj[which(soil_flowering_cor_df_padj$padj < 0.001),]
        
## Create co-occurrence network with igraph
flowering_soil_net <- graph_from_data_frame(soil_flowering_cor_df_padj,direct=F, vertices = flowering_soil_nodeattrib)

## Network properties ##
length(V(flowering_soil_net))

length(E(flowering_soil_net))

transitivity(flowering_soil_net)

graph.density(flowering_soil_net)

flowering_soil_net_cfg <- cluster_fast_greedy(as.undirected(flowering_soil_net))
flowering_soil_net_modules <- sort(table(membership(flowering_soil_net_cfg)),decr=T)

flowering_soil_net_betweenness <- betweenness(flowering_soil_net,normalized = T)
flowering_soil_net_degree <- sort(degree(flowering_soil_net),decreasing = T)
        
##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_soil <- layout_(flowering_soil_net,with_fr(niter=9999, grid="nogrid"))
        
## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_soil) <-  NULL
        
plot(flowering_soil_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_soil)


## root
root_flowering_occor <- rcorr(t(norm_root_flowering),type=c("spearman"))

## Create data frame of co-occurring OTUs
root_flowering_cor_df <- CorrDF(root_flowering_occor$r,root_flowering_occor$P)
root_flowering_cor_df$padj <- p.adjust(root_flowering_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
root_flowering_cor_df_padj <- root_flowering_cor_df[which(root_flowering_cor_df$cor > 0.7),]
root_flowering_cor_df_padj <- root_flowering_cor_df_padj[which(root_flowering_cor_df_padj$padj < 0.001),]

## Create co-occurrence network with igraph
flowering_root_net <- graph_from_data_frame(root_flowering_cor_df_padj,direct=F, vertices = flowering_root_nodeattrib)

## Network properties ##
length(V(flowering_root_net))

length(E(flowering_root_net))

transitivity(flowering_root_net)

graph.density(flowering_root_net)

flowering_root_net_cfg <- cluster_fast_greedy(as.undirected(flowering_root_net))
flowering_root_net_modules <- sort(table(membership(flowering_root_net_cfg)),decr=T)

flowering_root_net_betweenness <- betweenness(flowering_root_net,normalized = T)
flowering_root_net_degree <- sort(degree(flowering_root_net),decreasing = T)

##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_root <- layout_(flowering_root_net,with_fr(niter=9999, grid="nogrid"))

## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_root) <-  NULL

plot(flowering_root_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_root)


## film
film_flowering_occor <- rcorr(t(norm_film_flowering),type=c("spearman"))

## Create data frame of co-occurring OTUs
film_flowering_cor_df <- CorrDF(film_flowering_occor$r,film_flowering_occor$P)
film_flowering_cor_df$padj <- p.adjust(film_flowering_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
film_flowering_cor_df_padj <- film_flowering_cor_df[which(film_flowering_cor_df$cor > 0.7),]
film_flowering_cor_df_padj <- film_flowering_cor_df_padj[which(film_flowering_cor_df_padj$padj < 0.001),]

## Create co-occurrence network with igraph
flowering_film_net <- graph_from_data_frame(film_flowering_cor_df_padj,direct=F, vertices = flowering_film_nodeattrib)

## Network properties ##
length(V(flowering_film_net))

length(E(flowering_film_net))

transitivity(flowering_film_net)

graph.density(flowering_film_net)

flowering_film_net_cfg <- cluster_fast_greedy(as.undirected(flowering_film_net))
flowering_film_net_modules <- sort(table(membership(flowering_film_net_cfg)),decr=T)

flowering_film_net_betweenness <- betweenness(flowering_film_net,normalized = T)
flowering_film_net_degree <- sort(degree(flowering_film_net),decreasing = T)

##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_film <- layout_(flowering_film_net,with_fr(niter=9999, grid="nogrid"))

## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_film) <-  NULL

plot(flowering_film_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_film)


## CK
CK_flowering_occor <- rcorr(t(norm_CK_flowering),type=c("spearman"))

## Create data frame of co-occurring OTUs
CK_flowering_cor_df <- CorrDF(CK_flowering_occor$r,CK_flowering_occor$P)
CK_flowering_cor_df$padj <- p.adjust(CK_flowering_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
CK_flowering_cor_df_padj <- CK_flowering_cor_df[which(CK_flowering_cor_df$cor > 0.7),]
CK_flowering_cor_df_padj <- CK_flowering_cor_df_padj[which(CK_flowering_cor_df_padj$padj < 0.001),]

## Create co-occurrence network with igraph
flowering_CK_net <- graph_from_data_frame(CK_flowering_cor_df_padj,direct=F, vertices = flowering_CK_nodeattrib)

## Network properties ##
length(V(flowering_CK_net))

length(E(flowering_CK_net))

transitivity(flowering_CK_net)

graph.density(flowering_CK_net)

flowering_CK_net_cfg <- cluster_fast_greedy(as.undirected(flowering_CK_net))
flowering_CK_net_modules <- sort(table(membership(flowering_CK_net_cfg)),decr=T)

flowering_CK_net_betweenness <- betweenness(flowering_CK_net,normalized = T)
flowering_CK_net_degree <- sort(degree(flowering_CK_net),decreasing = T)

##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_CK <- layout_(flowering_CK_net,with_fr(niter=9999, grid="nogrid"))

## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_CK) <-  NULL

plot(flowering_CK_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_CK)


## LOW
low_flowering_occor <- rcorr(t(norm_low_flowering),type=c("spearman"))

## Create data frame of co-occurring OTUs
low_flowering_cor_df <- CorrDF(low_flowering_occor$r,low_flowering_occor$P)
low_flowering_cor_df$padj <- p.adjust(low_flowering_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
low_flowering_cor_df_padj <- low_flowering_cor_df[which(low_flowering_cor_df$cor > 0.7),]
low_flowering_cor_df_padj <- low_flowering_cor_df_padj[which(low_flowering_cor_df_padj$padj < 0.001),]

## Create co-occurrence network with igraph
flowering_low_net <- graph_from_data_frame(low_flowering_cor_df_padj,direct=F, vertices = flowering_low_nodeattrib)

## Network properties ##
length(V(flowering_low_net))

length(E(flowering_low_net))

transitivity(flowering_low_net)

graph.density(flowering_low_net)

flowering_low_net_cfg <- cluster_fast_greedy(as.undirected(flowering_low_net))
flowering_low_net_modules <- sort(table(membership(flowering_low_net_cfg)),decr=T)

flowering_low_net_betweenness <- betweenness(flowering_low_net,normalized = T)
flowering_low_net_degree <- sort(degree(flowering_low_net),decreasing = T)

##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_low <- layout_(flowering_low_net,with_fr(niter=9999, grid="nogrid"))

## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_low) <-  NULL

plot(flowering_low_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_low)


## high
high_flowering_occor <- rcorr(t(norm_high_flowering),type=c("spearman"))

## Create data frame of co-occurring OTUs
high_flowering_cor_df <- CorrDF(high_flowering_occor$r,high_flowering_occor$P)
high_flowering_cor_df$padj <- p.adjust(high_flowering_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
high_flowering_cor_df_padj <- high_flowering_cor_df[which(high_flowering_cor_df$cor > 0.7),]
high_flowering_cor_df_padj <- high_flowering_cor_df_padj[which(high_flowering_cor_df_padj$padj < 0.001),]

## Create co-occurrence network with igraph
flowering_high_net <- graph_from_data_frame(high_flowering_cor_df_padj,direct=F, vertices = flowering_high_nodeattrib)

## Network properties ##
length(V(flowering_high_net))

length(E(flowering_high_net))

transitivity(flowering_high_net)

graph.density(flowering_high_net)

flowering_high_net_cfg <- cluster_fast_greedy(as.undirected(flowering_high_net))
flowering_high_net_modules <- sort(table(membership(flowering_high_net_cfg)),decr=T)

flowering_high_net_betweenness <- betweenness(flowering_high_net,normalized = T)
flowering_high_net_degree <- sort(degree(flowering_high_net),decreasing = T)

##### Individual Co-occurence networks #####
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_flowering_high <- layout_(flowering_high_net,with_fr(niter=9999, grid="nogrid"))

## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_high) <-  NULL

plot(flowering_high_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_high)


## Properties of networks
network_properties <- read.csv("network.csv", header=T,stringsAsFactors=F, na.strings="NA",check.names=F)

network_properties$Community <-factor(network_properties$Community,c("Soil","Root","Film"))

barplot_otu <- barplot(network_properties$nodes, names.arg=network_properties$Concentration, 
                       las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                       main="nodes", col=CS_cols, border=F)

barplot_Connections <- barplot(network_properties$edges, names.arg=network_properties$Concentration, 
                               las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                               main="edges", col=CS_cols, border=F)

barplot_Density <- barplot(network_properties$Density, names.arg=network_properties$Concentration, 
                           las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                           main="Density", col=CS_cols, border=F)

barplot_Avg <- barplot(network_properties$`Avg. connectivity`, names.arg=network_properties$Concentration, 
                       las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                       main="Avg. connectivity", col=CS_cols, border=F)

barplot_Modules <- barplot(network_properties$Modules, names.arg=network_properties$Concentration, 
                           las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                           main="Modules", col=CS_cols, border=F)

barplot_Clustering <- barplot(network_properties$`Clustering coefficient`, names.arg=network_properties$Concentration, 
                              las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                              main="Clustering coefficient", col=CS_cols, border=F)

barplot_Modules <- barplot(network_properties$Modularity, names.arg=network_properties$Concentration, 
                           las=2, ylab="cumulative relative abundance", cex.lab=1, cex.axis=.7, cex.names=.7,
                           main="Modularity")


#####
        


##### circlize #####
PHYLUM_flowering_chord <- chordDiagram(PHYLUM_flowering,
                                       directional = TRUE,
                                       diffHeight = 0.06, 
                                       reduce = 0,
                                       transparency = 0.5)

##### weather #####
a <- read.csv("weather.csv",row.names = 1,header = T)

###2014
weather2014 <- a[which(a$Year=="2014"),]

weather2014_rain <- xyplot(Rainfall~ day,weather2014, type = "h" , lwd=2, col="black",ylim = c(0,100))
weather2014_temp <- xyplot(Maximum+Average+Minimum ~ day,weather2014, type = "l" , col=c( "#FB6542","#ffbb00","#375E97"), lwd=2.8,ylim = c(-40,40),cex=1)
weather2014_plot <- doubleYScale(weather2014_rain, weather2014_temp, add.ylab2 = TRUE, use.style=FALSE )

###2015
weather2015 <- a[which(a$Year=="2015"),]

weather2015_rain <- xyplot(Rainfall~ day,weather2015, type = "h" , lwd=2, col="black",ylim = c(0,100))
weather2015_temp <- xyplot(Maximum+Average+Minimum ~ day,weather2015, type = "l" , col=c("#FB6542","#ffbb00","#375E97"), lwd=2.8,ylim = c(-40,40),cex=1)
weather2015_plot <- doubleYScale(weather2015_rain, weather2015_temp, add.ylab2 = TRUE, use.style=FALSE )

###2016
weather2016 <- a[which(a$Year=="2016"),]

weather2016_rain <- xyplot(Rainfall~ day,weather2016, type = "h" , lwd=2, col="black",ylim = c(0,100))
weather2016_temp <- xyplot(Maximum+Average+Minimum ~ day,weather2016, type = "l" , col=c("#FB6542","#ffbb00","#375E97"), lwd=2.8,ylim = c(-40,40),cex=1)
weather2016_plot <- doubleYScale(weather2016_rain, weather2016_temp, add.ylab2 = TRUE, use.style=FALSE )

###2017
weather2017 <- a[which(a$Year=="2017"),]

weather2017_rain <- xyplot(Rainfall~ day,weather2017, type = "h" , lwd=2, col="black",ylim = c(0,100))
weather2017_temp <- xyplot(Maximum+Average+Minimum ~ day,weather2017, type = "l" , col=c("#FB6542","#ffbb00","#375E97"), lwd=2.8,ylim = c(-40,40),cex=1)
weather2017_plot <- doubleYScale(weather2017_rain, weather2017_temp, add.ylab2 = TRUE, use.style=FALSE )

###2018
weather2018 <- a[which(a$Year=="2018"),]

weather2018_rain <- xyplot(Rainfall~ day,weather2018, type = "h" , lwd=2, col="black",ylim = c(0,100))
weather2018_temp <- xyplot(Maximum+Average+Minimum ~ day,weather2018, type = "l" , col=c("#FB6542","#ffbb00","#375E97"), lwd=2.8,ylim = c(-40,40),cex=1)
weather2018_plot <- doubleYScale(weather2018_rain, weather2018_temp, add.ylab2 = TRUE, use.style=FALSE )

grid.newpage()
grid.arrange(weather2014_plot,weather2015_plot,weather2016_plot,weather2017_plot,weather2018_plot,ncol = 2)
##### 

##### Env #####
###  RDA
env <- read.csv("env.csv", header=T, row.names=1, stringsAsFactors=F, na.strings="NA",check.names=F)
env <- env[,-c(1:2)]
env_log <- log10(env)

design_flowering_rda <- design_flowering[-which(design_flowering$Site=="film"),]
flowering_otu_rda <- Family_mat_flowering[,rownames(design_flowering_rda)]
flowering_otu_rda <- decostand(flowering_otu_rda,method = "hellinger")
flowering_otu_rda <- t(flowering_otu_rda)
flowering_otu_rda <- as.data.frame(flowering_otu_rda, scale = FALSE)

flowering_rda <- rda(flowering_otu_rda~.,env_log)
flowering_rda_ill <- summary(flowering_rda)

vif.cca(flowering_rda)
anova(flowering_rda, by = "term", permutations=999)

RsquareAdj(rda(flowering_otu_rda~PH,data=env_log))$adj.r.squared
RsquareAdj(rda(flowering_otu_rda~DOC+NH+NO+P+MBC+MBN+MBP+PH,data=env_log))$adj.r.squared

st=as.data.frame(flowering_rda_ill$sites)[1:2]
bi=as.data.frame(flowering_rda_ill$biplot)[1:2]*3

rda_species <- flowering_rda_ill[["species"]]
rda_species <- as.data.frame(rda_species)
rda_species <- rda_species[order(rda_species[,1],rda_species[2]),]
top10 <- c("OTU57","OTU575","OTU392","OTU206","OTU55","OTU115","OTU384","OTU371","OTU423","OTU338")
rda_top10 <- rda_species[top10,]
rda_top10 <- as.data.frame(rda_top10)

ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = design_flowering_rda$Treatment1,size=4))+
  geom_point(data = rda_species,aes(RDA1,RDA2),size=3)+
  geom_segment(data = bi,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = bi,aes(RDA1,RDA2,label=row.names(bi)))+
  labs(x="RDA1 11.98%",y="RDA2 7.3%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())

###  Mental test
set_scale()

dist.env <- dist(env_log, method = 'euclidean')
dist.abund <- vegdist(flowering_otu_rda, method = 'bray')
abund_env <- mantel(dist.abund, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_env

mantel <- mantel_test(flowering_otu_rda,env_log,
                      spec.select = list(ck_root = 1:3,
                                         ck_soil = 4:6,
                                         low_root = 7:9,
                                         low_soil = 10:12,
                                         high_root = 13:15,
                                         high_soil = 16:18)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

quickcor(env_log, type = "lower") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
#####
