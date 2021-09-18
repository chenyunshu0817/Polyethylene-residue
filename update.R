

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
design <- read.table("metadata.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
tax <- read.table("taxonomy.txt", row.names=1, sep="\t", header=T ,stringsAsFactors=F,quote="")

##### Define sample types #####
floweringsamples <- rownames(design)[which(design$Stage == "Flowering")]

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

## Create bray-curtis dissimiliartiy matrix
all_dis_flowering <- vegdist(t(otu_table(physeq_norm_flowering)),method="bray")
        
##### PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_flowering <- adonis(all_dis_flowering ~Site+Concentration+Treatment1*Treatment2, data=design_flowering, permutations=9999)
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


##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylumar_flowering <- barplot(as.matrix(PHYLUM_mat_flowering_mean_treat), col=PHYLA_label_cols_flowering[rownames(PHYLUM_mat_flowering_mean_treat),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylumar_flowering, labels=colnames(PHYLUM_mat_flowering_mean_treat), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Flowering Community")
legend(9.5, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_flowering_legend$labels), 
       fill=rev(PHYLA_label_cols_flowering_legend$cols), 
       border=rev(PHYLA_label_cols_flowering_legend$cols) )

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

for(i in design_flowering[,9:12]) {
  fit <- aov(i~Treatment1,data = design_flowering)
  print(summary(fit))
  out <- LSD.test(fit,"Treatment1", p.adj="none" )
  print(out$groups)
}
#####

## Define Soil, Root and Film enriched BACTERIA
film_enrich_soil <- tt_enrich_soil_film$table[tt_enrich_soil_film$table$logFC < 0 & tt_enrich_soil_film$table$FDR < 0.05,]
soil_enrich_soil <- tt_enrich_soil_film$table[tt_enrich_soil_film$table$logFC > 0 & tt_enrich_soil_film$table$FDR < 0.05,]
        
plot(forMA_soil_film$CPM, forMA_soil_film$logFC, log="x", ylim=c(-8.2,7.6),
             main="Soil&Film",
             ylab="log fold change", xlab="average abundance CPM",
             col=forMA_soil_film$col, cex=1, pch=forMA_soil_film$pch, las=1)
#####        

##### Identifiying OTUs with indicator species analysis #####
## Define indicator species for soil bacteria community.
set.seed(8046)

indicatorsp_soil_flowering <- multipatt(indic_root_flowering,indic_root_groups_flowering,func = "r.g",control=how(nperm=9999))
CK_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.0 == 1 & indic_soil_df_flowering$p.value < 0.05),])
LOW_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.75 == 1 & indic_soil_df_flowering$p.value < 0.05),])
HIGH_soil_flowering <- as.matrix(indic_soil_df_flowering[which(indic_soil_df_flowering$s.600 == 1 & indic_soil_df_flowering$p.value < 0.05),])

pdf(paste0("soil_flowering_phylum_heartmap.pdf"),height=10, width=10)
par(las=1)
colors <- maPalette(l="lightgray",m="lightskyblue",h="orangered3", k=60)
heatmap.2(log2(soil_flowering_mso_phylum+1), col=colors, Rowv=F, Colv=F, margins=c(5,12),scale="none", trace="none"
          ,density.info="none",key.title=NA,srtCol=45,
          labRow=rownames(soil_flowering_mso_phylum),
          dendrogram="none",key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))
dev.off()

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
coords_flowering_film <- layout_(flowering_film_net,with_fr(niter=9999, grid="nogrid"))
        
## Import pre-calculated FR layout coordinates to save time 
dimnames(coords_flowering_film) <-  NULL
        
plot(flowering_soil_net, vertex.label=NA, vertex.size=3, layout=coords_flowering_soil)
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
