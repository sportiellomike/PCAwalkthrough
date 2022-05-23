# this script exists to help you easily create PCA plots, loading plots, and scree plots.
# this script was written by Mike Sportiello on 5.19.22

# load in libraries
library(ggplot2)
library(ggpubr)
library(viridis)
library(matrixStats)
library(dplyr)
library(uwot)
# Let's load in some fake data and take a peak
# dat<-mtcars # this data set is all numeric
dat<-iris # this data set is numeric with one column that is character (species column)
head(dat)
str(dat)
nums <- unlist(lapply(dat, is.numeric), use.names = FALSE)  

colnamesdat<-colnames(dat)
colnamesofchar<-colnamesdat[!nums]
# first we must remove columns that are not numeric since PCA requires only numeric values. Don't worry, we are going to add metadata character column back in later
metadata<-as.data.frame(dat[ , !nums])
colnames(metadata)<-colnamesofchar
dat<-dat[ , nums]

dat<-as.matrix(dat) #to use some downstream arguments, we have to have our data in the format of matrix, not data frame as it currently is
head(dat)
str(dat)
dat<-t(dat) # we want to transpose this so that each row is thing we might care about (for example gene name), and each column is a sample (in this case car type)

# as you can see, all of the cells are numeric 
#plot pca loadings

rv <- rowVars(dat)
select <- order(rv, decreasing=TRUE)[seq_len(min(500,length(rv)))]
subset <- dat[select, ]
pca <- prcomp(t(subset))
pcadimdf<-as.data.frame(pca$x) # this is where we actually look at what the value of each part of the pca is
#start to encode umap stuff
datdim<-dim(dat)
datdim[1]
ifelse(16>datdim[1],nneighbors<-datdim[1]-1,nneighbors<-15)
umapdimdf<-as.data.frame(umap(t(dat),n_neighbors = nneighbors))
length(umapdimdf)
umapcolnames<-paste0(rep('umap',length(umapdimdf)),1:length((umapdimdf)))
colnames(umapdimdf)<-umapcolnames

# recombine original numeric dat, pca dims, and metadata
megadat<-data.frame(t(dat),pcadimdf,umapdimdf,metadata) 

ggplot(megadat)+geom_point(aes(x=PC1,y=PC2)) # general pca
ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=Species)) # general pca colored by character column we originally took out species
ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=Species,size=Sepal.Length)) # general pca colored by character column we originally took out Species, and size is Sepal.Length
# ggsave(file.path("pcaplot.", Sys.Date(), ".png")), height = 4, width = 8, units = "in",dpi = 600)

ggplot(megadat)+geom_point(aes(x=umap1,y=umap2)) # general umap
ggplot(megadat)+geom_point(aes(x=umap1,y=umap2,color=Species)) # general umap colored by character column we originally took out species
ggplot(megadat)+geom_point(aes(x=umap1,y=umap2,color=Species,size=Sepal.Length)) # general umap colored by character column we originally took out Species, and size is Sepal.Length
# ggsave(file.path("umapplot.", Sys.Date(), ".png")), height = 4, width = 8, units = "in",dpi = 600)

# pca loadings plot
loadings <- as.data.frame(pca$rotation)
ifelse(length(loadings)<10,
       loadingcount<-length(loadings),
       loadingcount<-10)
princ<-loadings[1:loadingcount,]
absloadings<-abs(loadings)
absloadings<-absloadings[order(-absloadings$PC1),]
princ1names<-row.names(absloadings)

absloadings<-absloadings[order(-absloadings$PC2),]
princ2names<-row.names(absloadings)
absloadings<-absloadings[order(-absloadings$PC3),]
princ3names<-row.names(absloadings)




pc1<-ggplot(princ,aes(x=reorder(princ1names,PC1),y=PC1,fill=sign(PC1)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(option='plasma',begin = .2,end = .65)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC1',y='PC1 Contribution')+
  theme(legend.position = 'none')

pc2<-ggplot(princ,aes(x=reorder(princ2names,PC2),y=PC2,fill=sign(PC2)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(option='plasma',begin = .2,end = .65)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC2',y='PC2 Contribution')+
  theme(legend.position = 'none')

pc3<-ggplot(princ,aes(x=reorder(princ3names,PC3),y=PC3,fill=sign(PC3)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(option='plasma',begin = .2,end = .65)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC3',y='PC3 Contribution')+
  theme(legend.position = 'none')

gridpcloadings<-ggarrange(pc1,pc2,pc3,nrow=1)
gridpcloadings
# ggsave(file.path("pcloadings.", Sys.Date(), ".png")), plot = gridpcloadings, height = 4, width = 8, units = "in",dpi = 600)

## the contribution to the total variance for each component
percentVar <- 100*(pca$sdev^2 / sum( pca$sdev^2 ))
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:loadingcount)
colnames(scree_plot)<-c("Variance","component_number")
scree_plot$component_number <- factor(scree_plot$component_number, levels = c(1:12),labels = paste0('PC', c(1:12)))
pcscree<- ggplot(scree_plot, mapping=aes(x=component_number, y=Variance))+       
  geom_bar(stat="identity")+
  labs(x='Component Number')+
  theme(
    plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
    axis.title = element_text(color='black'),
    axis.text = element_text(color='black'),
    text = element_text(face= 'bold',color='black'))
pcscree
# ggsave(file.path("./", paste0("pcscree.integrinskept.", Sys.Date(), ".png")), plot = pcscree, height = 3, width = 8, units = "in",dpi = 600)

# add umap stuff
megaumap<-as.data.frame(umap(megadat))
numericmega<-megadat[,1:4]
nummegumap<-as.data.frame(umap(numericmega))

ggplot(megaumap)+geom_point(aes(x=V1,y=V2))
ggplot(nummegumap)+geom_point(aes(x=V1,y=V2))

sessionInfo()
######################################################################################

# old stuff can delete

# loadings$genes<-row.names(loadings)

# loadings <- as.data.frame(pca$rotation)
# loadings$genes<-row.names(loadings)
# 
# loadings<-loadings[order(-loadings$PC1),]
# head(loadings)
# princ1<-loadings[1:loadingcount,]
# loadings<-loadings[order(loadings$PC1),]
# princ1b<-loadings[1:loadingcount,]
# head(princ1)
# head(princ1b)
# princ1<-rbind(princ1,princ1b)
# loadings<-loadings[order(-loadings$PC2),]
# head(loadings)
# princ2<-loadings[1:loadingcount,]
# loadings<-loadings[order(loadings$PC2),]
# princ2b<-loadings[1:loadingcount,]
# head(princ2)
# head(princ2b)
# princ2<-rbind(princ2,princ2b)
# loadings<-loadings[order(-loadings$PC3),]
# head(loadings)
# princ3<-loadings[1:loadingcount,]
# loadings<-loadings[order(loadings$PC3),]
# princ3b<-loadings[1:loadingcount,]
# head(princ3)
# head(princ3b)
# princ3<-rbind(princ3,princ3b)
#plots
# pc1<-ggplot(princ1,aes(x=reorder(princ1names,PC1),y=PC1,fill=sign(PC1)))+
#   geom_col()+coord_flip()+
#   scale_fill_viridis(option='plasma',begin = .2,end = .65)+
#   theme(axis.title.y = element_blank())+
#   theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
#         axis.title = element_text(color='black'),
#         axis.text = element_text(color='black'),
#         text = element_text(face= 'bold',color='black'))+
#   labs(title='PC1',y='PC1 Contribution')+
#   theme(legend.position = 'none')


# pc2<-ggplot(princ2,aes(x=reorder(princ2names,PC2),y=PC2,fill=sign(PC2)))+geom_col()+coord_flip()+scale_fill_viridis(option='plasma',begin = .2,end = .65)+
#   theme(axis.title.y = element_blank())+
#   theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
#         axis.title = element_text(color='black'),
#         axis.text = element_text(color='black'),
#         text = element_text(face= 'bold',color='black'))+
#   labs(title='PC2',y='PC2 Contribution')+
#   theme(legend.position = 'none')
# pc3<-ggplot(princ3,aes(x=reorder(princ3names,PC3),y=PC3,fill=sign(PC3)))+geom_col()+coord_flip()+scale_fill_viridis(option='plasma',begin = .2,end = .65)+
#   theme(axis.title.y = element_blank())+
#   theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
#         axis.title = element_text(color='black'),
#         axis.text = element_text(color='black'),
#         text = element_text(face= 'bold',color='black'))+
#   labs(title='PC3',y='PC3 Contribution')+
#   theme(legend.position = 'none')