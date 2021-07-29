library(dplyr)
library(dada2)
library(phyloseq)
set.seed(1234)

path <- "cleandata/" 
outpath <- "Result/"
files.list <- read.csv("sample-metadata.tsv", sep="\t",row.names=1, header=T)
#################
list.files(path)

fnFs <- paste0(path, rownames(files.list),"_clean_1.fq")
fnRs <-paste0(path, rownames(files.list),"_clean_2.fq")

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

filtFs <- file.path(outpath, "02_filtered", paste0(rownames(files.list), "_F_filt.fq.gz"))
filtRs <- file.path(outpath, "02_filtered", paste0(rownames(files.list), "_R_filt.fq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,240),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=T)
head(out)
#
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#ASV
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#
getN <- function(x) sum(getUniques(x))
#
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100, 1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","total_perc_reads_lost")
rownames(track) <- rownames(files.list)
head(track)


taxa <- assignTaxonomy(seqtab.nochim, "gg_13_8_train_set_97.fa.gz", multithread=T, tryRC=T)
taxa.print <- taxa 

# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- rownames(files.list)
write.table(asv_tab, "ASVs_counts.txt", sep="\t", quote=F)

asv_tax <- taxa.print
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.txt", sep="\t", quote=F)


########################
### diversity function
alpha <- function(x, tree = NULL, base = exp(1)) {
  
  Richness <- specnumber(x)
  Shannon <- diversity(x, index = 'shannon')
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 
  Pielou <- Shannon / log(Richness)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, goods_coverage)
  
  if (!is.null(tree)) {
    
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    
    names(PD_whole_tree) <- 'PD_whole_tree'
    
    result <- cbind(result, PD_whole_tree)
    
  }
  
  result
  
}


##############
library(vegan)
library(ade4)
library(ggpubr)
library(ape)

beta <- function(data = dat.data, pheno = dat.pheno, Group=Group){
  datNorm <-decostand(data,"hell")
  datNorm.bray <- vegdist(datNorm, method="bray")
  
  is.euclid(datNorm.bray)
  is.euclid(sqrt(datNorm.bray))
  datNorm.braysq <- sqrt(datNorm.bray)
  
  mds <- cmdscale(datNorm.braysq, k=2, eig=TRUE)
  mds_point <- data.frame(mds$points)
  
  
  ####
  mds_pheno <- cbind(mds_point,pheno);
  
  return (list (datNorm.braysq = datNorm.braysq, mds_pheno = mds_pheno, eig = mds$eig));
}

filter_zero <- function(species.pheno = species.pheno) {species.pheno.filter <- species.pheno  %>% rownames_to_column() %>% 
  dplyr::mutate(sum_spe =  rowSums(.[2:dim(species.pheno)[2]] ) ) %>% 
  column_to_rownames() %>% subset(sum_spe > 0)  %>% dplyr::select (-sum_spe) %>% t() %>% as.data.frame()#总丰度>0.001%

species.filter.sample <- species.pheno.filter %>% rownames_to_column() %>% 
  dplyr::mutate(sum_id =  rowSums(.[2:dim(species.pheno.filter)[2]] ) ) %>% 
  column_to_rownames() %>% subset(sum_id > 0)  %>% dplyr::select (-sum_id)

return(species.filter.sample)
}
pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
  
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    x <- dat.data[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] 
    y <- factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    
    datNorm <-decostand(x,"hell")
    datNorm.bray <- vegdist(datNorm, method="bray")
    
    datNorm.braysq <- sqrt(datNorm.bray)
    set.seed(1234)
    ad = adonis(datNorm.braysq ~ y );
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  
}

####################
asv_table <- cbind(asv_tax, asv_tab)
asv_genus <- data.frame(Names = paste0(asv_tax[,5],";",asv_tax[,6]), asv_tab)
asv_genus.count <- aggregate(. ~ Names, data = asv_genus, sum)
asv_genus.relative <- decostand(asv_genus.count[,-1], 'total',2) %>%  `rownames<-`(asv_genus.count[,1]) %>%
                      rownames_to_column(var = "rowname") %>%
                      filter(!str_detect(rowname, ";NA")) %>%
                      filter(!str_detect(rowname, ";g__$")) %>%
                      column_to_rownames(var = "rowname")
write.table(asv_genus.relative, "ASVs_relactive_abundance.txt", sep="\t", quote=F)

################
# alpha diversity
asv_genus.alpha <- asv_genus.relative %>% t() %>% as.data.frame()
alpha.result <- alpha(asv_genus.alpha)
alpha.plot <- data.frame(alpha.result, files.list)

chose <- subset(alpha.plot, Source=="Fecal")
aov1 <- aov(Shannon~Transfer, chose)
aov.result <- summary(aov1)
aov.result
tukey = TukeyHSD(x=aov1, 'Transfer', conf.level=0.95)
tukey

alpha.plot$Transfer <- factor(alpha.plot$Transfer, levels =c("Normal","I_III","IV"),ordered = T)
pdf("Shannon_fecal.pdf",width=3,height=2.5)
p <- ggplot(data = alpha.plot, aes(x = Transfer, y = as.numeric(Shannon), fill = Transfer)) +
      stat_boxplot(geom=('errorbar'), width=.4) + geom_boxplot() +
      scale_fill_manual( values = c("#7EC0EE","#1C86EE","#0000FF")) +
      geom_signif ( comparisons = list(c("Normal","I_III"), c("I_III","IV"), c("Normal","IV")), 
                map_signif_level =T, test = t.test, step_increase = 0.1)#map_signif_level = F,则显示p值
p  +  labs( y = "Shannon", color = NULL) + theme_bw()+theme(panel.grid=element_blank(), legend.position = "none")

dev.off()

#################
##Beta diversity

asv_genus.relative[asv_genus.relative < .0001] <- 0

dat <- asv_genus.relative %>% dplyr::select(rownames(files.list)) %>% filter_zero() %>% 
            as.data.frame() %>% cbind(files.list)

dat.compare <- dat %>% subset(Source %in% "Sputum") ###
dat.data <- dat.compare %>% dplyr::select(-c("Source","Brain","Disease","Stage","Group","Transfer"));
dat.pheno <- dat.compare  %>% dplyr::select(Group = Disease);###
dat.beta <- beta(data = dat.data, pheno = dat.pheno);

set.seed(1234)
adonis(dat.beta$datNorm.braysq~Group, dat.pheno)

##############
#sputum  "#ffbcbc","#ff7f7f","#FF0000"
#fecal "#7EC0EE","#1C86EE","#0000FF"
dat.beta$mds_pheno$Group <- factor(dat.beta$mds_pheno$Group, levels =c("Normal","I_III","IV"),ordered = T)
p <- ggscatter(dat.beta$mds_pheno, x= "X1", y = "X2", 
               color = "Group",shape = "Group",size=3, palette = c("#ffbcbc","#7EC0EE"),   
               ellipse = TRUE,  
               ellipse.level = 0.65
                ) +
      labs(title="Sputum CONvsNSCLC", 
          x = paste("PCoA 1 (", format(100*dat.beta$eig[1]/sum(dat.beta$eig), digits = 4), "%)",sep = ""), 
          y = paste("PCoA 2 (", format(100*dat.beta$eig[2]/sum(dat.beta$eig), digits = 4), "%)",sep = "")) +
      annotate("text", x=0.2, y=0.3, label="Pvalue = 0.001***\nR^2 = 0.36",fontface="italic", family="serif", colour="black", size=4) +
      scale_shape_manual(values=c(16,15,17))
pp <- p + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA))+
      theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10), legend.position = c(0.9,0.1),legend.background=element_rect(
      fill = "white",
      colour = "black", 
      size = .5 ) )
pp
pdf("pic/PCoA_Sputum_disease.pdf",width=5,height=5)
print(pp)
dev.off()

pairwise.adonis(dat.data, dat.pheno$Group,  p.adjust.m= "fdr")
###################
#### Top 10 Genus
library(tidyr)
dat.tidy <- asv_genus.relative %>% rownames_to_column(var = "Genus") %>%  gather(key= "SampleID", value = "abundance", -Genus) 
match <- match(dat.tidy$SampleID, rownames(files.list))
dat.tidy.meta <- cbind(dat.tidy, files.list[match,])


#group_by Disease / Transfer
tidy.disease <- dat.tidy.meta %>% group_by( Source, Disease, Genus) %>% summarise( Mean = mean(abundance))  %>% 
  mutate(r.mean=rank(-abs(Mean)))
genera.top10 <- tidy.disease %>% subset(r.mean < 11) 

pdf("top10_disease.pdf")
tidy.disease %>%
  subset(Genus %in% unique(genera.top10$Genus) ) %>%
  # highlight feature rank for the top 20 features
  mutate(r.mean=case_when(r.mean > 10~NA_real_, TRUE~r.mean)) %>% 
  mutate(Groups =paste0(Source,"_" ,Disease)) %>%
  ggplot(aes(y=Genus, x=Groups, fill=Mean)) +
  geom_tile() +
    scale_fill_gradientn(colours=rev(
    c('#BE5400', "#FFA300", "#EFC06E", 'white'))) +
  theme_minimal() +
  geom_text(aes(label=r.mean), col='black', size= 6) +
  theme(panel.grid = element_blank()) +
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=8))
dev.off()

#####
library(ggalluvial)
library(RColorBrewer)
fecal_trans.top10 <- genera.top10 %>% subset(Source=="Fecal") %>% as.data.frame()
colors <- getPalette(length(unique(fecal_trans.top10$Genus)))
fecal_trans.top10$Transfer <- factor(fecal_trans.top10$Transfer, levels = c("Normal", "I_III","IV"), order=T) 


p <- ggplot(fecal_trans.top10, aes(x=Transfer, y=Mean, fill=Genus, alluvium = Genus))+
      geom_bar(stat="identity", width = 0.5)+
      scale_fill_manual(values = colors)+
      geom_alluvium(aes(fill = Genus))+
      labs(title = "Fecal transfer",y = "Genus mean of relactive abundance")
pp <- p + theme_bw() + theme(panel.grid=element_blank(), legend.position = "right")
pp
pdf("Sputum_transfer_top10.pdf", height = 5)
print(pp)
dev.off()


############################
#lefse analysis
library(microbiomeMarker)
#seqtab.nochim
seq_tab <- readRDS(system.file("extdata", "dada2_seqtab.rds",
                               package= "microbiomeMarker"))#dada2::removeBimeraDenovo()
#taxa
tax_tab <- readRDS(system.file("extdata", "dada2_taxtab.rds",
                               package= "microbiomeMarker"))#dada2::assignTaxonomy() or dada2::addSpecies()
#files.list
sam_tab <- read.table(system.file("extdata", "dada2_samdata.txt",
                                  package= "microbiomeMarker"), sep = "\t", header = TRUE, row.names = 1)

rownames(seqtab.nochim) <- rownames(files.list)
taxa_nospecies <- taxa[,-7]
ps <- import_dada2(seq_tab = seqtab.nochim, 
                   tax_tab = taxa_nospecies, 
                   sam_tab = files.list)

sputum.ps <- phyloseq::subset_samples(
  ps,
  Source %in% c("Sputum")
)


mm <- lefse(
  sputum.ps,
  class="Disease",
  subclass = NULL,
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multicls_strat = T,
  sample_min = 10,
  only_same_subcls = FALSE,
  curv = FALSE
)

plot_ef_bar(mm, label_level = 1) +
  scale_fill_manual(values = c("Normal" = "blue", "NSCLC" = "red"))

library(data.table)

feature <-  marker_table(mm)[grep(pattern = "\\|g__\\S", perl = T, marker_table(mm)$feature),]
feature$feature <- unlist(lapply(strsplit(as.character(feature$feature), split="__", fixed=T), 
                                  function(data){y <- paste(data[6],"__",data[7], sep="") }))

p <- ggplot(data=feature, mapping=aes(x=feature,y= lda,fill= enrich_group))+
      geom_bar(stat="identity",size=1,width=0.7)+
      coord_flip()+
      scale_x_discrete(limits=rev(unique(feature$feature)))+
      guides(fill=guide_legend(reverse=F))+
      labs(y="LDA score (log10)", x="Enrichment genus") + scale_fill_manual(values = c("#3CB371B2","#DAA520B2")) 
genus.plot <- p + theme_set(theme_classic()) + theme(panel.grid.major=element_line(colour=NA)) +
      theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))+
      theme(axis.title.x=element_text(size=9,face="bold"), axis.title.y=element_text(size=9,face="bold")) +
      theme(legend.title = element_text(size=9, face="bold"),legend.position = "top") +
      theme(legend.text = element_text(size =6,face="bold"))
genus.plot







