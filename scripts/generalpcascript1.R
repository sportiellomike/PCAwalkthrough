# this script exists to help you easily create PCA plots, umap plots, pacmap plots, loading plots, and scree plots.
# this script was written by Mike Sportiello on 5.19.22
# this script was edited by Mike Sportiello on 8.18.22

# make sure that you have both python installed as well as the called imports if you want to run the pacmap stuff.
# i recommend installing python through miniconda: https://docs.conda.io/en/latest/miniconda.html 
# i think you may have to assign which python rstudio is using. To do so, go to Tools, then global options, then python, then click which one you want to use.
# if you are fine with only PCA and umap, you can comment out/ delete all the pacmap section, 
# and then make sure you aren't adding something called pacmapdat to your megadat
# alrightalrightalright let's start this script.

# load in libraries
library(ggplot2)
library(ggpubr)
library(viridis)
library(matrixStats)
library(dplyr)
library(parallel)
library(uwot)
library(factoextra)
library(NbClust)
library(reticulate)
library(readxl)
# Let's load in some fake data and take a peak
# dat<-mtcars # this data set is all numeric
dat<-iris # this data set is numeric with one column that is character (species column)
# dat<-read_xlsx('./exampledatasetforpcascript.xlsx') 
# change iris to your own data if you want to run this script for your own purposes. You may need to read.csv() or read_xlsx() it in and assing it to dat.
head(dat) # let's look at what the first few columns look like
str(dat) # let's look at the actual structure: which columns are numeric? Which are character? etc.
nums <- unlist(lapply(dat, is.numeric), use.names = FALSE)  # let's look at which columns are numeric and therefore able to be entered into downstream PCA and other dimensional reduction techniques
nums # as you can see, for the iris dataset, the first four columns are numeric, and the last is not numeric.

colnamesdat<-colnames(dat) # let's save the colnames to its own character vector so we can subset specific colnames later.
colnamesofchar<-colnamesdat[!nums] # here let's save the portion of colnamesdat that are not numeric, and save it as a thing called colnameschar (char was short for character)
colnamesofchar # as you can see, the only colname that saved here is the one that is not numeric.

# first we must remove columns that are not numeric since PCA requires only numeric values. Don't worry, we are going to add metadata character column back in later
metadata<-as.data.frame(dat[ , !nums])

colnames(metadata)<-colnamesofchar # this is going to rename the metadata colnames appropriately

dat<-dat[ , nums] # this subsets the dataset dat to only include the columns that are numeric, as saved in our nums vector earlier

dat<-as.matrix(dat) #to use some downstream arguments, we have to have our data in the format of matrix, not data frame as it currently is
head(dat) # let's take a look at the first few lines

str(dat) # let's take a look at the actual structure
dat<-t(dat) # we want to transpose this so that each row is thing we might care about (for example gene name), and each column is a sample (in this case plant type)
head(dat) # let's take a look at the first few lines

str(dat) # let's take a look at the actual structure

# as you can see, all of the cells are numeric 
#plot pca loadings

rv <- rowVars(dat) # this is looking at the overal variance for each row. Remember, PCAs or principal component analyses try to maximize the amount of variance captured within each principal commponent starting with PC1 and being able to capture slightly less variance each component thereafter.

select <- order(rv, decreasing=TRUE)[seq_len(min(500,length(rv)))] # this line is going to reorder rv from most variance to least. The second part of this line starting with 'seq_len' is basically subsetting rv to a certain size. in this case, we are subsetting rv to whichever is shorter: 500 long, or however long rv actually is. In the case of this iris dataset it is only 4 long, so it is not subsetting it at all. If our dataset was 1000 long, then we would subset rv to only the first 500 rows.

subset <- dat[select, ] # this is the actual subsetting step described above
pca <- prcomp(t(subset)) # this step actually calls prcomp which is the PCA argument. 
pca
# View(pca) # let's actually look at what pca contains. Yes it contains the dataframe of "principal component numbers" we would want for plotting, but there's other important informatio we might want for doing more in depth analyses as well.

pcadimdf<-as.data.frame(pca$x) # this is where we actually look at what the value of each part of the pca is
# so we have now finished creating and saving the prinicipal components for this dataset. Let's also run a umap on this data in case that is more relevant to your analysis.

#start to encode UMAP stuff
datdim<-dim(dat) # let's look at how many dimensions we have
datdim # so we have 4 rows and 150 colums
datdim[1] # this is another way to quickly look at the number of rows (basically it's telling the computer to show the first part of the output that we saved into dat dim which is the number of rows)

ifelse(16>datdim[1],nneighbors<-datdim[1]-1,nneighbors<-15) # this argument is important because the way UMAP works is that is requires you to input paramaters, including nneighbors. You should read and understand what this means, but I can't do it justice in this script. One thing you have to know is that you can't have more nneighbors than you do rows, so this line is an ifelse statement that says if there are fewer 16 rows in your dat, then we are going to use one fewer than however many rows you have for nneighbors. Otherwise, it's going to use exactly 15 nneighbors whether you have 17 or 1 billion rows.

umapdimdf<-as.data.frame(umap(t(dat),n_neighbors = nneighbors,n_threads = parallel::detectCores()-1)) # this line calls the umap function and saves the relevant output as umapdimdf

# umapdimdf<-as.data.frame(umap(t(dat),n_neighbors = nneighbors) # if your computer keeps crashing when running the above umap, run this one instead, which is slightly different in the way the funtion sends info to your CPU and memory and may help with that. Don't run this as it is slower unless you have to.


length(umapdimdf) # this tells us how many umapdims it made. we can set it to a higher number if we want to I think but generally you only ever want the first two which is the default
umapcolnames<-paste0(rep('umap.',length(umapdimdf)),1:length((umapdimdf))) # let's change the colnames to something more meaningful
colnames(umapdimdf)<-umapcolnames

# let's get the pacmap stuff going
# you should probably read some about it: https://github.com/YingfanWang/PaCMAP
# these next few lines are basically running R as python, and it won't make much sense to your eyes if you know R or python, so I'm not going to comment the below stuff out.

python_pandas <- import("pandas")
python_pacmap <- import("pacmap")
python_numpy <- import("numpy")
tdat<-as.data.frame(t(dat))
pydat <- reticulate::r_to_py(tdat)
nparray <- pydat$values
nparray <- nparray$astype(python_numpy$float)
embedding <- python_pacmap$PaCMAP(n_dims = as.integer(2), MN_ratio=0.5, FP_ratio=2.0)  # if you want more than 2 pacmap components you can change 2 to the number you actually want.

X_transformed <- embedding$fit_transform(nparray, init="pca")
pacmapdat <- data.frame(X_transformed)
colnames(pacmapdat)<-paste0('pacmap.',seq(1:length(colnames(pacmapdat))))
head(pacmapdat)
# let's cluster

# first let's scale the data
tdatscaled <- scale(tdat)

# now we can try to find the optimal number of clusters so we aren't just guessing. Here are presented a number of ways.
# you should probably read more about that here: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
set.seed(585414) # we set a seed to make it reproducible. You can pick any number for this.
# Elbow method
fviz_nbclust(tdatscaled, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(tdatscaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# gap method
fviz_nbclust(tdatscaled, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

resultNbClust<-NbClust(tdatscaled, diss=NULL, distance = "euclidean", min.nc=2, max.nc=15, 
        method = "kmeans", index = "all") 


tresultnbclust<-as.data.frame(t(resultNbClust$Best.nc))
ggplot(tresultnbclust)+geom_bar(aes(x=factor(Number_clusters)))+
  theme(panel.grid.major = element_line(color='light grey'))+
  xlab('Number of proposed closters')+
  ylab('Number of indices that voted for this number of clusters')


ggplot(tresultnbclust)+geom_col(aes(x=rownames(tresultnbclust),y=Number_clusters))+
  # scale_x_continuous(limits = c())+
  theme(panel.grid.major = element_line(color='light grey'),
        axis.text.x = element_text(angle=90))
# now that we picked the numberof clusters, set it to k
k<-6 # set this number to your desired number of clusters
km.res <- kmeans(tdat, centers = k, nstart = 25)
cluster<-as.factor(km.res$cluster)

# recombine original numeric dat, pca dims, pacmapdat, umapdat, and metadata. 
# This will be your main dataframe from which you plot now.
megadat<-data.frame(t(dat),pcadimdf,pacmapdat,umapdimdf,cluster,metadata) 



# let's do some plotting!

# first let's set some themes
theme_set(theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_rect(colour = 'black',fill = NA)))

# now we can plot
ggplot(megadat)+geom_point(aes(x=PC1,y=PC2)) # general pca
ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=Species)) + scale_color_viridis(discrete = T,end=.8)# general pca colored by character column we originally took out species
# if you want to change color, you can delete scale color viridis and use scale_color_manual(c('colrcodes'))

ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=cluster)) +
  scale_color_viridis(discrete = T,end=.8)+
  stat_conf_ellipse(aes(x=PC1,y=PC2,color=cluster)) # general pca colored by cluster, which correlates almost perfectly with species as you'll notice

ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=Species))+
  stat_conf_ellipse(aes(x=PC1,y=PC2,color=Species),level=.95) +
  scale_color_viridis(discrete = T,end=.8)# This is the same thing as the last plot but we added a confidence ellipse of a confidence interval of 95% as well. 

ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=Species,size=Sepal.Length)) + scale_color_viridis(discrete = T,end=.8) # general pca colored by character column we originally took out Species, and size is Sepal.Length
# ggsave(file.path("pcaplot-.", Sys.Date(), ".png")), height = 4, width = 8, units = "in",dpi = 600)

ggplot(megadat)+geom_point(aes(x=umap.1,y=umap.2)) # general umap
ggplot(megadat)+geom_point(aes(x=umap.1,y=umap.2,color=Species)) + scale_color_viridis(discrete = T,end=.8) # general umap colored by character column we originally took out species
ggplot(megadat)+geom_point(aes(x=umap.1,y=umap.2,color=cluster)) + scale_color_viridis(discrete = T,end=.8) # general umap colored by character column we originally took out species
ggplot(megadat)+geom_point(aes(x=umap.1,y=umap.2,color=Species,size=Sepal.Length)) # general umap colored by character column we originally took out Species, and size is Sepal.Length

ggplot(megadat)+geom_point(aes(x=pacmap.1,y=pacmap.2))  # general pacmap
ggplot(megadat)+geom_point(aes(x=pacmap.1,y=pacmap.2,color=Species))+stat_conf_ellipse(aes(x=pacmap.1,y=pacmap.2,color=Species),level=.95)+ scale_color_viridis(discrete = T,end=.8) # general pacmap colored by character column we originally took out species

ggplot(megadat)+geom_point(aes(x=pacmap.1,y=pacmap.2,color=cluster))+stat_conf_ellipse(aes(x=pacmap.1,y=pacmap.2,color=cluster),level=.95)+ scale_color_viridis(discrete = T,end=.8) # general pacmap colored by character column we originally took out cluster
ggplot(megadat)+geom_point(aes(x=pacmap.1,y=pacmap.2,color=Species,size=Sepal.Length)) # general pacmap colored by character column we originally took out Species, and size is Sepal.Length

# ggsave(file.path("umapplot-.", Sys.Date(), ".png")), height = 4, width = 8, units = "in",dpi = 600)





pccluster<-ggplot(megadat)+geom_point(aes(x=PC1,y=PC2,color=cluster)) +
  scale_color_viridis(discrete = T,end=.8)+
  stat_conf_ellipse(aes(x=PC1,y=PC2,color=cluster))
umapcluster<-ggplot(megadat)+geom_point(aes(x=umap.1,y=umap.2,color=cluster)) +
  stat_conf_ellipse(aes(x=umap.1,y=umap.2,color=cluster),level=.95)+
  scale_color_viridis(discrete = T,end=.8)
pacmapcluster<-ggplot(megadat)+geom_point(aes(x=pacmap.1,y=pacmap.2,color=cluster))+stat_conf_ellipse(aes(x=pacmap.1,y=pacmap.2,color=cluster),level=.95)+ scale_color_viridis(discrete = T,end=.8)


ggarrange(pccluster,umapcluster,pacmapcluster,nrow = 1,common.legend = T,legend = 'top')


# pca loadings plot
loadings <- as.data.frame(pca$rotation) # this is the part that actually tells us how much variance is captured in each principal component


ifelse(length(loadings)<10, # these arguments basically look at how big to make the loadings plot. we only have four rows so we are going to make it the same length
       loadingcount<-length(loadings),
       loadingcount<-10)
princ<-loadings[1:loadingcount,]
absloadings<-abs(loadings) # this is taking absolute value of the loadings since a very positive or very negative value both contribute highly to the PCA
absloadings<-absloadings[order(-absloadings$PC1),] # now we reorder the PC1 loadings
princ1names<-row.names(absloadings) # and we save the rownames in princ1names

absloadings<-absloadings[order(-absloadings$PC2),] # we do same thing for PC2 and PC3
princ2names<-row.names(absloadings)
absloadings<-absloadings[order(-absloadings$PC3),]
princ3names<-row.names(absloadings)

# now we plot the contributions we have for each saved in princ, but we reorder princ by princ1names
pc1<-ggplot(princ,aes(x=reorder(rownames(princ),PC1),y=PC1,fill=sign(PC1)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(begin = .1,end = .8)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC1',y='PC1 Contribution')+
  theme(legend.position = 'none')
# we make similar plots for pc2 and pc3
pc2<-ggplot(princ,aes(x=reorder(rownames(princ),PC2),y=PC2,fill=sign(PC2)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(begin = .1,end = .8)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC2',y='PC2 Contribution')+
  theme(legend.position = 'none')

pc3<-ggplot(princ,aes(x=reorder(rownames(princ),PC3),y=PC3,fill=sign(PC3)))+
  geom_col()+coord_flip()+
  scale_fill_viridis(begin = .1,end = .8)+
  theme(axis.title.y = element_blank())+
  theme(plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
        axis.title = element_text(color='black'),
        axis.text = element_text(color='black'),
        text = element_text(face= 'bold',color='black'))+
  labs(title='PC3',y='PC3 Contribution')+
  theme(legend.position = 'none')

gridpcloadings<-ggarrange(pc1,pc2,pc3,nrow=1) # we then put these plots in one plot
gridpcloadings # let's look at the plot
# ggsave(file.path("pcloadings.", Sys.Date(), ".png")), plot = gridpcloadings, height = 4, width = 8, units = "in",dpi = 600)

## the contribution to the total variance for each component
percentVar <- 100*(pca$sdev^2 / sum( pca$sdev^2 ))
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:loadingcount)
colnames(scree_plot)<-c("Variance","component_number")
scree_plot$component_number <- factor(scree_plot$component_number, levels = c(1:length(scree_plot$component_number)),labels = paste0('PC', c(1:length(scree_plot$component_number))))

scree_plot$cumulativevariance <- cumsum(scree_plot$Variance)
pcscreewithcumulline<- ggplot(scree_plot, mapping=aes(x=component_number, y=Variance))+       
  geom_bar(stat="identity")+
  geom_line(aes(x=component_number,y=cumulativevariance,group=1))+
  labs(x='Component Number')+
  theme(
    plot.title = element_text(color="black", size=14, face="bold",hjust = .5),
    axis.title = element_text(color='black'),
    axis.text = element_text(color='black'),
    text = element_text(face= 'bold',color='black'))
pcscreewithcumulline

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

sessionInfo()

### FIN ###
# Any questions or comments? head to https://github.com/sportiellomike
