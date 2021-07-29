library(SIAMCAT)
library(dplyr)
library(randomForest)

compare.groups <- list("I_III.IV" = c("I_III","IV"),"CON.IV"= c("Normal","IV"),"CON.I_III"=c("Normal","I_III"),"CON.NSCLC"=c("Normal","NSCLC"));
names <- c("I_III.IV","CON.IV","CON.I_III")
#siamcat function
siamcat.cross <- function(x,y) {
  siamcat(feat = x, meta = y,
          label='Transfer', case=compare.groups[[d]][2]   # change The compare groups;
  ) %>% 
    filter.features( filter.method = 'abundance',
                     cutoff = 0.001,
                     rm.unmapped = T,
                     verbose=2
    ) %>%  
    normalize.features(norm.method = "log.std",
                       norm.param = list(log.n0 = 1e-05, sd.min.q = 0.1),
                       verbose = 2
    ) %>% 
    create.data.split( num.folds = 10,
                       num.resample = 10
    ) %>% 
    train.model( method = "randomForest" 
    ) %>%
    make.predictions(normalize.holdout = F) %>% evaluate.predictions()
}
siamcat.deg <- function(x,y) {
  siamcat(feat = x, meta = y,
          label='Transfer', case=compare.groups[[d]][2]   # change The compare groups;
  ) %>% 
    filter.features( filter.method = 'abundance',
                     cutoff = 0.001,
                     rm.unmapped = TRUE,
                     verbose=2
    ) %>%  
    normalize.features(norm.method = "log.std",
                       norm.param = list(log.n0 = 1e-05, sd.min.q = 0.1),
                       verbose = 2
    ) %>%
    check.associations(
      sort.by = 'fc',
      alpha = 0.1,
      mult.corr = "fdr",
      detect.lim = 10 ^-6,
      plot.type = "box",
      panels = c("fc", "prevalence", "auroc"),
      max.show = 20,
      feature.type = "normalized"
    )
}
###
# machine learning
pheno <- read.csv("sample-metadata.tsv", sep="\t",row.names=1, header=T)

f.dat <- asv_genus.relative %>% dplyr::select( paste0(rownames(pheno),"_F")) 
rownames(f.dat) <- paste0("F:", rownames(f.dat))
colnames(f.dat) <- rownames(pheno)
s.dat <- asv_genus.relative %>% dplyr::select( paste0(rownames(pheno),"_S")) 
rownames(s.dat) <- paste0("S:", rownames(s.dat))
colnames(s.dat) <- rownames(pheno)
#####################
siamcat.sputum.train <- list()
siamcat.sputum.test <- list()
siamcat.sputum.train.obj <- list()
siamcat.sputum.train.deg <- list()
siamcat.sputum.test.deg <- list()
caret.sputum <- list()
siamcat.rfe.sputum.train <- list()
siamcat.rfe.sputum.test <- list()
set.seed(1234)

for (d in names){
    all.meta.data <- subset( s.meta.data, Transfer %in% compare.groups[[d]])
  for(num in 1:10){
    i <- paste0(d,"_",num)
 
    training.samples <- all.meta.data$Transfer %>%  createDataPartition(p = 0.7, list = FALSE)
    meta.data <- all.meta.data[training.samples,]
    feat.data <- s.feat.data[, rownames(meta.data)] %>% as.data.frame()
    
    siamcat.sputum.train[[i]] <- siamcat.cross(x=feat.data, y=meta.data)
    
    siamcat.sputum.train.obj[[i]] <- siamcat.deg(x=feat.data, y=meta.data)
    temp <- associations(siamcat.sputum.train.obj[[i]]) %>% subset(p.adj <0.1)
    genus <- rownames(temp)
    deg.feat <- subset(feat.data, rownames(feat.data) %in% genus)
    siamcat.sputum.train.deg[[i]] <- siamcat.cross(x=deg.feat, y=meta.data)
    
    feat.x <- s.feat.data[,rownames(all.meta.data)] %>% t() %>% as.data.frame() 
    feat.x <- scale(feat.x[,-nearZeroVar(feat.x)])
    feat.x <- feat.x[,-findCorrelation(cor(feat.x), .8)] %>% as.data.frame()
    feat.x.data <- feat.x[training.samples,]
      caret.sputum[[i]] <-  rfe(feat.x.data, as.factor(meta.data$Transfer), sizes = c(2:50), 
                                rfeControl = rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 10, allowParallel = TRUE))
      
      rfe.genus <- caret.sputum[[i]]$optVariables
      rfe.feat <- subset(s.feat.data, rownames(s.feat.data) %in% rfe.genus)
      train.feat <- rfe.feat[,rownames(meta.data)]
      siamcat.rfe.sputum.train[[i]] <- siamcat.cross(x=train.feat, y=meta.data)
      
      test.meta <-  all.meta.data[-training.samples,]
      test.rfe.feat <- rfe.feat[,rownames(test.meta)]
      siamcat.rfe.sputum.test[[i]] <- siamcat(feat = test.rfe.feat, meta = test.meta,
                                              label='Transfer', case=compare.groups[[d]][2])
      siamcat.rfe.sputum.test[[i]] <-  make.predictions(siamcat.rfe.sputum.train[[i]], siamcat.rfe.sputum.test[[i]]) 
      siamcat.rfe.sputum.test[[i]] <-  evaluate.predictions(siamcat.rfe.sputum.test[[i]])
      
      test.feat <- s.feat.data[, rownames(test.meta)] %>% as.data.frame()
      siamcat.sputum.test[[i]] <- siamcat(feat = test.feat, meta = test.meta,
                                          label='Transfer', case=compare.groups[[d]][2])
      siamcat.sputum.test[[i]] <-  make.predictions(siamcat.sputum.train[[i]], siamcat.sputum.test[[i]]) 
      siamcat.sputum.test[[i]] <-  evaluate.predictions(siamcat.sputum.test[[i]])
      
      test.deg.feat <- subset(test.feat, rownames(test.feat) %in% genus)
      siamcat.sputum.test.deg[[i]] <- siamcat(feat = test.deg.feat, meta = test.meta,
                                              label='Transfer', case=compare.groups[[d]][2])
      siamcat.sputum.test.deg[[i]] <-  make.predictions(siamcat.sputum.train.deg[[i]], siamcat.sputum.test.deg[[i]]) 
      siamcat.sputum.test.deg[[i]] <-  evaluate.predictions(siamcat.sputum.test.deg[[i]])
    }
}


siamcat.fecal.train <- list()
siamcat.fecal.test <- list()
siamcat.fecal.train.obj <- list()
siamcat.fecal.train.deg <- list()
siamcat.fecal.test.deg <- list()
caret.fecal <- list()
siamcat.rfe.fecal.train <- list()
siamcat.rfe.fecal.test <- list()


set.seed(1234)
for (d in names){
  all.meta.data <- subset( f.meta.data, Transfer %in% compare.groups[[d]])
  for(num in 1:10){
    i <- paste0(d,"_",num)

    training.samples <- all.meta.data$Transfer %>%  createDataPartition(p = 0.7, list = FALSE)
    meta.data <- all.meta.data[training.samples,]
    feat.data <- f.feat.data[, rownames(meta.data)] %>% as.data.frame()
    
    siamcat.fecal.train[[i]] <- siamcat.cross(x=feat.data, y=meta.data)
    
    siamcat.fecal.train.obj[[i]] <- siamcat.deg(x=feat.data, y=meta.data)
    temp <- associations(siamcat.fecal.train.obj[[i]]) %>% subset(p.adj <0.1)
    genus <- rownames(temp)
    deg.feat <- subset(feat.data, rownames(feat.data) %in% genus)
    siamcat.fecal.train.deg[[i]] <- siamcat.cross(x=deg.feat, y=meta.data)
    
    feat.x <- f.feat.data[,rownames(all.meta.data)] %>% t() %>% as.data.frame() 
    feat.x <- scale(feat.x[,-nearZeroVar(feat.x)])
    feat.x <- feat.x[,-findCorrelation(cor(feat.x), .8)] %>% as.data.frame()
    
    feat.x.data <- feat.x[training.samples,]
    caret.fecal[[i]] <-  rfe(feat.x.data, as.factor(meta.data$Transfer), sizes = c(2:50), 
                             rfeControl = rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 10, allowParallel = TRUE))
    
    rfe.genus <- caret.fecal[[i]]$optVariables
    rfe.feat <- subset(f.feat.data, rownames(f.feat.data) %in% rfe.genus)
    train.feat <- rfe.feat[,rownames(meta.data)]
    siamcat.rfe.fecal.train[[i]] <- siamcat.cross(x=train.feat, y=meta.data)
    
    test.meta <-  all.meta.data[-training.samples,]
    test.rfe.feat <- rfe.feat[,rownames(test.meta)]
    siamcat.rfe.fecal.test[[i]] <- siamcat(feat = test.rfe.feat, meta = test.meta,
                                           label='Transfer', case=compare.groups[[d]][2])
    siamcat.rfe.fecal.test[[i]] <-  make.predictions(siamcat.rfe.fecal.train[[i]], siamcat.rfe.fecal.test[[i]]) 
    siamcat.rfe.fecal.test[[i]] <-  evaluate.predictions(siamcat.rfe.fecal.test[[i]])
    
    test.feat <- f.feat.data[, rownames(test.meta)] %>% as.data.frame()
    siamcat.fecal.test[[i]] <- siamcat(feat = test.feat, meta = test.meta,
                                       label='Transfer', case=compare.groups[[d]][2])
    siamcat.fecal.test[[i]] <-  make.predictions(siamcat.fecal.train[[i]], siamcat.fecal.test[[i]]) 
    siamcat.fecal.test[[i]] <-  evaluate.predictions(siamcat.fecal.test[[i]])
    
    test.deg.feat <- subset(test.feat, rownames(test.feat) %in% genus)
    siamcat.fecal.test.deg[[i]] <- siamcat(feat = test.deg.feat, meta = test.meta,
                                           label='Transfer', case=compare.groups[[d]][2])
    siamcat.fecal.test.deg[[i]] <-  make.predictions(siamcat.fecal.train.deg[[i]], siamcat.fecal.test.deg[[i]]) 
    siamcat.fecal.test.deg[[i]] <-  evaluate.predictions(siamcat.fecal.test.deg[[i]])
  }
}


siamcat.train <- list()
siamcat.test <- list()
siamcat.rfe.train <- list()
siamcat.rfe.test <- list()
caret.list <- list()
set.seed(1234)
for (d in names){
  all.meta.data <- subset( pheno, Transfer %in% compare.groups[[d]])
  f.feat.all <- f.dat[, rownames(all.meta.data)] %>% as.data.frame()
  s.feat.all <- s.dat[, rownames(all.meta.data)] %>% as.data.frame()
  feat.all <- rbind(s.feat.all,f.feat.all)
  feat.x <- t(feat.all) %>% as.data.frame() 
  feat.x <- scale(feat.x[,-nearZeroVar(feat.x)])
  feat.x <- feat.x[,-findCorrelation(cor(feat.x), .8)] %>% as.data.frame()
  for(num in 1:10){
    i <- paste0(d,"_",num)
  
  training.samples <- all.meta.data$Transfer %>%  createDataPartition(p = 0.9, list = FALSE)
  meta.data <- all.meta.data[training.samples,]
  feat.data <- feat.all[,rownames(meta.data)]
  
  siamcat.train[[i]] <- siamcat.cross(x=feat.data, y=meta.data)
  test.meta <-  all.meta.data[-training.samples,]
  test.feat <- feat.all[, rownames(test.meta)] %>% as.data.frame()
  
  siamcat.test[[i]] <- siamcat(feat = test.feat, meta = test.meta,
                               label='Transfer', case=compare.groups[[d]][2]
  )
  siamcat.test[[i]] <-  make.predictions(siamcat.train[[i]], siamcat.test[[i]]) 
  siamcat.test[[i]] <-  evaluate.predictions(siamcat.test[[i]])
  
  feat.x.data <- feat.x[training.samples,]
  caret.list[[i]] <- rfe(feat.x.data, as.factor(meta.data$Transfer), sizes = c(2:50), 
                         rfeControl = rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 10, allowParallel = TRUE))
  
  genus <- caret.list[[i]]$optVariables
  rfe.feat <- subset(feat.all, rownames(feat.all) %in% genus)
  train.feat <- rfe.feat[,rownames(meta.data)]
  siamcat.rfe.train[[i]] <- siamcat.cross(x=train.feat, y=meta.data)
 
  test.feat <-  rfe.feat[,rownames(test.meta)]
  
  siamcat.rfe.test[[i]] <- siamcat(feat = test.feat, meta = test.meta,
                                   label='Transfer', case=compare.groups[[d]][2]
  )
  siamcat.rfe.test[[i]] <-  make.predictions(siamcat.rfe.train[[i]], siamcat.rfe.test[[i]]) 
  siamcat.rfe.test[[i]] <-  evaluate.predictions(siamcat.rfe.test[[i]])
    }
}

##
# mean of AUC

library(pROC)
roc.result <- c()
for(y in 1:10)
{
  roc.result <- rbind(roc.result,auc(eval_data(siamcat.test[[y]])$roc))
}
mean(roc.result)
##########
# plot of feature importance and heatmap
feature_mean <- foreach(y = 1:10, .combine = rbind) %do% { 
  feature_weights(siamcat.train[[y]]) %>% rownames_to_column(var = "Genus") %>%
    mutate(r.mea=rank(-abs(mean.weight))) %>% arrange(r.mea) %>% 
    mutate(r.mea.top=case_when(r.mea > 20~0, TRUE~1))
}  %>% dplyr::select(Genus, mean.weight, r.mea.top) %>% group_by(Genus) %>%
  mutate(num=sum(r.mea.top), weight_mean = mean(mean.weight)) %>% 
  mutate (robustness = (num/10)*100) %>%
  arrange(robustness,weight_mean) 

feature_mean$Source <- ifelse(grepl("F:",feature_mean$Genus),"Fecal","Sputum")
feature.plot <- feature_mean %>% subset(Genus %in% tail(unique(feature_mean$Genus),n=20L) )

p <- ggplot( data = feature.plot, aes( x= Genus, y= mean.weight) )+
  geom_boxplot(aes(fill = Source))+
  coord_flip()+
  scale_x_discrete(limits=unique(feature.plot$Genus))+
  geom_text(data = feature.plot, aes( x= Genus, y= max(mean.weight)+0.01,label=paste0(robustness,"%"))) +
  labs( title = "Importance Score",  y = "Mean weights in accuracy", x=NULL, fill="Source")
genus.plot <- p + theme_set(theme_classic()) + theme(panel.grid.major=element_line(colour=NA)) +
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+ 
  theme(axis.title.x=element_text(size=14,face="bold"), 
        axis.title.y=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"))+
  theme(legend.position=c(0.7,0.15)) +
  scale_fill_manual(values=c("#6D9EC1","#E46726"))
print (genus.plot)


genus.select <- rev(unique(feature.plot$Genus))
pheno.chose <- pheno %>% subset(Transfer %in% c("I_III","IV")) %>% arrange(Transfer)
genus.abun <- all.dat[genus.select, rownames(pheno.chose)] #%>% t() %>% cbind(pheno.chose) %>% gather(key="Genus", value="abundance", -Disease, -Transfer) 

table(pheno.chose$Transfer)
anno_col = data.frame(Status = c(rep("I_III",27),rep("IV",38)))
rownames(anno_col) <- colnames(genus.abun)

ann_colors = list(
  Status = c( I_III = "#6CC24A", IV= "#007A53") 
);

heat.map <- pheatmap(genus.abun,scale = "row",cluster_row = F, cluster_col = F, border=NA, show_colnames = T,
                     fontsize_number = 12, number_color = "black",
                     color= c(colorRampPalette(colors = c("#3B6FB6","white","#D41645"))(100)),
                     annotation_col = anno_col, annotation_colors = ann_colors, 
                     main="Heatmap") %>% as.ggplot()
print(heat.map)
pdf(file = "all/feature_Mix_IIIvsDM.pdf", width = 15, height = 7);
plot_grid(genus.plot,heat.map)
dev.off();

require(ggplotify)



  
  #test.meta <-  all.meta.data[-training.samples,]
  test.feat <-  rfe.feat[,rownames(test.meta)]
  
  siamcat.rfe.test[[i]] <- siamcat(feat = test.feat, meta = test.meta,
                                   label='Disease', case=compare.groups[[d]][2]
  )
  siamcat.rfe.test[[i]] <-  make.predictions(siamcat.rfe.train[[i]], siamcat.rfe.test[[i]]) 
  siamcat.rfe.test[[i]] <-  evaluate.predictions(siamcat.rfe.test[[i]])
}
