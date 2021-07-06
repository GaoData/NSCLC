library(dplyr)
library(dada2)
library(phyloseq)
set.seed(1234)
setwd("/mnt/raid5/nagao/Data/Xiehe/qiime2/New_202103")
path <- "/mnt/raid5/nagao/Data/Xiehe/qiime2/01_cleandata/" 
outpath <- "/mnt/raid5/nagao/Data/Xiehe/qiime2/New_202103"
files.list <- read.csv("sample-metadata.tsv", sep="\t",row.names=1, header=T)
#################
list.files(path)
# 返回测序正向文件完整文件名
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnFs <- paste0(path, rownames(files.list),"_clean_1.fq")

# 返回测序反向文件完整文件名
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs <-paste0(path, rownames(files.list),"_clean_2.fq")

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# basename(fnFs)提取文件名(不要目录)
# strsplit按`_`进行分割字符，返回列表
# sapply批量操作，这里批量提取列表中第一个元素，即样本名
# 提取文件名中`_`分隔的第一个单词作为样品名
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#这里用files.list代替sample.names

#绘制前4个样本的质量示意图
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

# Place filtered files in filtered/ subdirectory
# 将过滤后的文件存于filtered子目录，设置输出文件名

filtFs <- file.path(outpath, "02_filtered", paste0(rownames(files.list), "_F_filt.fq.gz"))
filtRs <- file.path(outpath, "02_filtered", paste0(rownames(files.list), "_R_filt.fq.gz"))
# 过滤文件输出，输出和参数，统计结果保存于Out
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,240),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=T)
head(out)
#错误率
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#去除重复序列
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#基于错误模型进一步质控
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#序列拼接
#默认情况下，仅当正向和反向序列重叠至少12个碱基并且在重叠区域中彼此相同时才输出合并序列
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#生成ASV表
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 查看序列长度分布
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#去除嵌合体
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#求和函数
getN <- function(x) sum(getUniques(x))
#合并各样本分步数据量
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100, 1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#修改列名
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","total_perc_reads_lost")
#行名修改为样本名
rownames(track) <- rownames(files.list)
#统计结果预览
head(track)


#taxa <- assignTaxonomy(seqtab.nochim, "/mnt/raid1/data/databases/qiime/silva_nr99_v138.1_train_set.fa.gz", multithread=T, tryRC=T)
taxa <- assignTaxonomy(seqtab.nochim, "/mnt/raid1/data/databases/qiime/gg_13_8_train_set_97.fa.gz", multithread=T, tryRC=T)
#taxa <- addSpecies(taxa,  "/mnt/raid1/data/databases/qiime/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa 
# 另存物种注释变量，去除序列名，只显示物种信息
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


library(vegan)
#导入的数据集中每行为一个样本，每例为一个物种或者OTU，且导入数据应为原始注释结果
# 行名为sample 列名为feature

otu = t(asv_tab)
#统计每个样本的物种数或OTU数目
S <- specnumber(otu) 
#确定在样本中所包含的最小的reads条数
raremax <- min(rowSums(otu)) #3,114
#计算一定reads数下预期的物种丰富度
Srare <- rarefy(otu, raremax)
#预期物种丰富度与实际物种丰富度的分布
plot(S, Srare, xlab = "Observed No. of Genus", ylab = "Rarefied No. of Genus")
pdf("rarecurve.pdf")
rarecurve(otu, step = 500, col = rainbow(10),xlab = "ASVs Number", ylab = "Genus Richness",label = F)
dev.off()

pdf("rarecurve_3114.pdf")
rarecurve(otu, step = 500, sample=raremax, col = rainbow(10),xlab = "ASVs Number", ylab = "Genus Richness",label = F)
dev.off()

pdf("rarecurve_5235.pdf")
rarecurve(otu, step = 500, sample=5235, col = rainbow(10),xlab = "ASVs Number", ylab = "Genus Richness",label = F)
dev.off()
##################

library(dplyr)
asv_tab <- read.table( "ASVs_counts.txt",sep="\t",header=T,row.names = 1)
asv_tab_filt <- asv_tab %>% dplyr::select(!YDQ_F)
otu1 = otu_table(asv_tab_filt, taxa_are_rows = T)
phyloseq = phyloseq(otu1)

set.seed(1234)
#这种方法会自动去除一些低丰度的otu
rare.data = rarefy_even_depth(phyloseq,replace = TRUE)
sample_sums(phyloseq)
sample_sums(rare.data)

#提取抽平后的otu表格 
rare.otu = rare.data@.Data %>% as.data.frame()
rare.tax <- asv_tax %>% subset(rownames(asv_tax) %in% rownames(rare.otu))
rare.table <- rare.otu %>% subset(rownames(rare.otu) %in% rownames(asv_tax)) %>% bind_cols(rare.tax)
rare_genus <- data.frame(Names = paste0(asv_tax[,5],";",asv_tax[,6]), rare.otu)
asv_genus.count <- aggregate(. ~ Names, data = asv_genus, sum)
asv_genus.relative <- decostand(asv_genus.count[,-1], 'total',2) %>%  `rownames<-`(asv_genus.count[,1]) %>%
  rownames_to_column(var = "rowname") %>%
  filter(!str_detect(rowname, ";g__$")) %>%
  column_to_rownames(var = "rowname")

########################

alpha <- function(x, tree = NULL, base = exp(1)) {
  
  Richness <- specnumber(x)
  Shannon <- diversity(x, index = 'shannon')
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
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

test = rarefy_even_depth(phyloseq,replace = F, sample.size = 5235, rngseed = 1234)
rare.otu = test@.Data %>% as.data.frame() %>% filter_zero() %>% as.data.frame
 
alpha.result <- alpha(rare.otu)
rare.pheno <- files.list %>% subset(rownames(files.list) %in% rownames(alpha.result))

alpha.plot <- data.frame(alpha.result, rare.pheno)


chose <- subset(alpha.plot, Source=="Sputum")
aov1 <- aov(Shannon~Transfer, chose)
aov.result <- summary(aov1)
aov.result
tukey = TukeyHSD(x=aov1, 'Transfer', conf.level=0.95)
tukey

tukey = as.data.frame(tukey$Transfer)
tukey$pair = rownames(tukey)
p <- ggplot(tukey, aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1), label=c("p<0.01","p<0.05","Non-Sig")))) +
      theme_bw(base_size = 16)+
      geom_hline(yintercept=0, lty="11", colour="grey30",size = 1) +
      geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2,size = 1) +
      geom_point(aes(pair, diff),size = 2) +
      labs(colour="")+
      theme(axis.text.x = element_text(size = 14))


dat.asv <- data.frame(rare.otu, rare.pheno)
dat.compare <- dat.asv %>% subset(Source %in% "Sputum") ###选定比较组数据
dat.data <- dat.compare %>% dplyr::select(-c("Source","Brain","Disease","Stage","Group","Transfer"));
dat.pheno <- dat.compare  %>% dplyr::select(Group = Transfer);###改比较组
dat.beta <- beta(data = dat.data, pheno = dat.pheno);

set.seed(1234)
adonis(dat.beta$datNorm.braysq~Group, dat.pheno)



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



asv_genus.relative[asv_genus.relative < .0001] <- 0

dat <- asv_genus.relative %>% dplyr::select(rownames(files.list)) %>% filter_zero() %>% 
            as.data.frame() %>% cbind(files.list)

dat.compare <- dat %>% subset(Source %in% "Sputum") ###选定比较组数据
dat.data <- dat.compare %>% dplyr::select(-c("Source","Brain","Disease","Stage","Group","Transfer"));
dat.pheno <- dat.compare  %>% dplyr::select(Group = Disease);###改比较组
dat.beta <- beta(data = dat.data, pheno = dat.pheno);

set.seed(1234)
adonis(dat.beta$datNorm.braysq~Group, dat.pheno)

##############
#sputum  "#ffbcbc","#ff7f7f","#FF0000"
#fecal "#7EC0EE","#1C86EE","#0000FF"
dat.beta$mds_pheno$Group <- factor(dat.beta$mds_pheno$Group, levels =c("Normal","I_III","IV"),ordered = T)
p <- ggscatter(dat.beta$mds_pheno, x= "X1", y = "X2", 
               color = "Group",shape = "Group",size=3, palette = c("#ffbcbc","#7EC0EE"),  ####分组颜色 
               ellipse = TRUE,  
               ellipse.level = 0.65#,
               #ggtheme = theme_minimal()
                ) +
      labs(title="Sputum CONvsNSCLC", #####图片标题
          x = paste("PCoA 1 (", format(100*dat.beta$eig[1]/sum(dat.beta$eig), digits = 4), "%)",sep = ""), 
          y = paste("PCoA 2 (", format(100*dat.beta$eig[2]/sum(dat.beta$eig), digits = 4), "%)",sep = "")) +
      annotate("text", x=0.2, y=0.3, label="Pvalue = 0.001***\nR^2 = 0.36",fontface="italic", family="serif", colour="black", size=4) +
      scale_shape_manual(values=c(16,15,17))
pp <- p + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA))+
      theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10), legend.position = c(0.9,0.1),legend.background=element_rect(
      fill = "white", # 填充色
      colour = "black", # 框线色
      size = .5 ) )
pp
pdf("pic/PCoA_Sputum_disease.pdf",width=5,height=5)
print(pp)
dev.off()

pairwise.adonis(dat.data, dat.pheno$Group,  p.adjust.m= "fdr")
###################
library(tidyr)
# gather()命令转换说明：
# gather（data=数据框名，key="key名"，value="value名"，要转换的列1，列2，列3）
# gene_exp_tidy <- gather(data = gene_exp, key = "sample_name", value = "expression", Sample1, Sample2, Sample3)
# 在指定要转换的列时，也可不用列名，直接指定列的编号即可
# gene_exp_tidy <- gather(data = gene_exp, key = "sample_name", value = "expression", 2:4)
#  在指定要转换的列时，也可指定不需转换的列，其他列参与转换
# gene_exp_tidy <- gather(data = gene_exp, key = "sample_name", value = "expression", -GeneId)

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
  #scale_fill_gradient(low = "#EFC06E",
   #                    mid =  "#FFA300",
    #                   high = "#BE5400", na.value = NA) +
  scale_fill_gradientn(colours=rev(
    c('#BE5400', "#FFA300", "#EFC06E", 'white'))) +
  theme_minimal() +
  geom_text(aes(label=r.mean), col='black', size= 6) +
  theme(panel.grid = element_blank()) +
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=8))
dev.off()


ggplot(as.data.frame(Titanic),
       aes(y = Freq,
           axis1 = Survived, axis2 = Sex, axis3 = Class)) +
  geom_alluvium(aes(fill = Class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Survived", "Sex", "Class")) +
  coord_flip() +
  ggtitle("Titanic survival by class and sex")

#####
library(ggalluvial)
library(RColorBrewer)
fecal_trans.top10 <- genera.top10 %>% subset(Source=="Fecal") %>% as.data.frame()
colors <- getPalette(length(unique(fecal_trans.top10$Genus)))
fecal_trans.top10$Transfer <- factor(fecal_trans.top10$Transfer, levels = c("Normal", "I_III","IV"), order=T) 


p <- ggplot(fecal_trans.top10, aes(x=Transfer, y=Mean, fill=Genus, alluvium = Genus))+
      geom_bar(stat="identity", width = 0.5)+
      #theme_bw()+
      #guides(fill=F)+
      scale_fill_manual(values = colors)+
      geom_alluvium(aes(fill = Genus))+
      #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      labs(title = "Fecal transfer",y = "Genus mean of relactive abundance")
pp <- p + theme_bw() + theme(panel.grid=element_blank(), legend.position = "right")
pp
pdf("Sputum_transfer_top10.pdf", height = 5)
print(pp)
dev.off()


############################
#lefse
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
  #transform = c("identity", "log10", "log10p"),
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multicls_strat = T,
  #correct = c("0", "1", "2"),
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

####STAMP
two_group_welch <- test_two_groups(
  enterotypes_arumugam, 
  group = "Gender", 
  method = "welch.test"#welch test, t test and white test.
)

multiple_group_anova <-  test_multiple_groups(
  ps,
  group = "Enterotype", 
  method = "anova"#anova and kruskal test
)
pht <- posthoc_test(ps, group = "Enterotype")
markers <- marker_table(multiple_group_anova)$feature
plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")





