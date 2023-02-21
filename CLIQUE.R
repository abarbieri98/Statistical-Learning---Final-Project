library(MASS)
library(subspace)
library(ggplot2)

##### Util functions ######

generate_dataset_gaussian<-function(n, concentration = 1, distance = 0, caos = 0){
  cov_1 <- 0
  cov_2 <- 0
  cov_3 <- 0
  if(caos != 0){
    cov_1 <- runif(1,-1,1)
    cov_2 <- runif(1,-1,1)
    cov_3 <- runif(1,-1,1)
  }
  
  sigma<-matrix(c((runif(1,1,10)*caos)+(1/concentration),cov_1,cov_1,(runif(1,1,10)*caos)+(1/concentration)),2,2)
  sigma2<-matrix(c((runif(1,1,10)*caos)+(1/concentration),cov_2,cov_2, (runif(1,1,10)*caos)+(1/concentration)),2,2)
  sigma3<-matrix(c((runif(1,1,10)*caos)+(1/concentration),cov_3,cov_3,(runif(1,1,10)*caos)+(1/concentration)),2,2)
  
  samples <- matrix(rep(0,n*2), ncol= 2)
  for (i in 1:n) {
    coin <- runif(1,0,3)
    val <- 0
    if(coin < 1){
      val <- mvrnorm(1, c(2,3+distance),sigma)
    }
    else if(coin > 2){
      val <- mvrnorm(1, c(4+distance,-3),sigma2)
    }
    else{
      val <- mvrnorm(1, mu = c(-1-distance,-2),Sigma = sigma3)
    }
    
    samples[i,] <- val
  }
  return(samples)
}

grids_2d<-function(n, data){
  ticks_x1<-seq(min(data[,1]), max(data[,1]), length.out = n)
  ticks_x2<-seq(min(data[,2]), max(data[,2]), length.out = n)
  abline(v=ticks_x1)
  abline(h=ticks_x2)
}

biggest_clust = function(model){
  common_dim <- rep(0, length(model))
  for(i in 1:length(model)){
    common_dim[i] <-sum(model[[i]]$subspace)
  }
  
  selected_clus = which(common_dim==max(common_dim))
  
  clust_size <- rep(0,length(selected_clus))
  for(i in 1:length(selected_clus)){
    clust_size[i] <-length(model[[selected_clus[i]]]$objects)
  }
  cat("Maximum dimensionality: ", max(common_dim), "\n")
  n_clust = length(selected_clus)
  cat("Clusters found: ", n_clust, "\n")
  cat("Cluster size: \n")
  for(i in 1:length(selected_clus)){
    cat(clust_size[i], " ")
  }
  cat("\nCluster number: \n")
  for(i in 1:length(selected_clus)){
    cat(selected_clus[i], " ")
  }
}

select_cluster = function(model, d, size=0){
  common_dim <- rep(0, length(model))
  for(i in 1:length(model)){
    common_dim[i] <-sum(model[[i]]$subspace)
  }
  
  selected_clus = which(common_dim==d)
  
  clust_size <- rep(0,length(selected_clus))
  for(i in 1:length(selected_clus)){
    clust_size[i] <-length(model[[selected_clus[i]]]$objects)
  }
  if(size !=0){
    selected_clus = selected_clus[which(clust_size > size)]
    clust_size <- rep(0,length(selected_clus))
    for(i in 1:length(selected_clus)){
      clust_size[i] <-length(model[[selected_clus[i]]]$objects)
    }
  }
  
  cat("Selected dimensionality: ",d , "\n")
  n_clust = length(selected_clus)
  cat("Clusters found: ", n_clust, "\n")
  cat("Cluster size: \n")
  for(i in 1:length(selected_clus)){
    cat(clust_size[i], " ")
  }
  cat("\nCluster number: \n")
  for(i in 1:length(selected_clus)){
    cat(selected_clus[i], " ")
  }
}
data <- data.frame(generate_dataset_gaussian(1000,1, distance = 10, caos = 1))
plot(data)

xi<-10
tau<- .005
model <- CLIQUE(data = data, xi = xi, tau = tau)
plot(data)
grids_2d(xi, data)
summary(model)
plot(model, data, "mix")
plot(model,data, "any")

View(model)
# For the specified parameters only two clusters are bi-dimensional

clus_1 <- model[[3]]$objects
clus_2 <- model[[4]]$objects
clus_3 <- model[[5]]$objects
clus_4 <- model[[6]]$objects
plot(data)

points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)

###### Circles #####

data <- read.csv("circles.csv")
data<-data[,2:3]
plot(data)

xi<-10
tau<- .1
model <- CLIQUE(data = data, xi = xi, tau = tau)
plot(data)
grids_2d(xi, data)
summary(model)

View(model)
# For the specified parameters only two clusters are bi-dimensional

clus_1 <- model[[2]]$objects
clus_2 <- model[[5]]$objects
clus_3 <- model[[6]]$objects
clus_4 <- model[[11]]$objects
plot(data)

points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)
points(data[clus_5,], col = 6)
points(data[clus_6,], col = 7)


##### Iris dataset ######

data <- iris
data <- data[,1:4]

model <- CLIQUE(data = data, xi = 10, tau = .01)
summary(model)
plot(model,data, "mix")
plot(model,data, "any")

library(cluster.datasets)

data("new.haven.school.scores")

matrix <-data.frame(new.haven.school.scores[,2:ncol(new.haven.school.scores)])
cluster<- CLIQUE(matrix, 5, .3)
summary(cluster)
plot(cluster,matrix)
plot(cluster, matrix, "any")

# Looking at the model information we can see that only a cluster is three-dimensional,
# hence we extract its units and plot them on the first two dimensions

obs <- as.integer( cluster[[13]]$objects)
cluster_3 <- matrix[obs,]
plot(matrix[,1],matrix[,2])
points(cluster_3[,1], cluster_3[,2], col=2)


###### UN2 dataset #####

data <- read.csv("un2.csv")
# Optimal parameters: xi = 15, tau = 0.004
data <- data[,1:2]
plot(data)

xi<-15
tau<- .004
model <- CLIQUE(data = data, xi = xi, tau = tau)
plot(data)
grids_2d(xi, data)
summary(model)

View(model)

clus_1 <- model[[2]]$objects
clus_2 <- model[[4]]$objects
clus_3 <- model[[5]]$objects
clus_4 <- model[[6]]$objects
clus_5 <- model[[7]]$objects
plot(data)

points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)
points(data[clus_5,], col = 6)


###### Spotify data #####

rm(list=ls())
data <- read.csv("spotify_data.csv", header =  T, sep = ",")
data <- data[,-1]
data_string <- data[,14:16]
data <- data[,-c(9,14:16)]
data_sd <-scale(data,center = T,scale = T)

# Dimensionality Reduction

pca = prcomp(data_sd)
summary(pca)

# Seems rough but let's take the first 5 components

data_red = pca$x[,1:5]
# Clustering

# We first check if clusters may be present
library(hopkins)
hopkins(data_red)

# We now apply the CLIQUE algorithm

xi<-5
tau<- 0.1 
model <- CLIQUE(data = data_sd, xi = xi, tau = tau)
summary(model)
View(model)

biggest_clust(model)
select_cluster(model,6,700)
clust_size

clus_1 <- model[[7]]$objects
clus_2 <- model[[45]]$objects
clus1_obs <- data_sd[clus_1,]
clus2_obs <- data_sd[clus_2,]
clus1_names <- data_string[clus_1,]
clus2_names <- data_string[clus_2,]
clus1_names
summary(clus1_obs)
summary(clus2_obs)
boxplot(clus1_obs)
diag(var(clus1_obs))/diag(var(data_sd))
diag(var(clus2_obs))

dim1 = 2
dim2 = 4
plot(data_sd[,dim1], data_sd[,dim2], xlab = colnames(data_sd)[dim1], ylab = colnames(data_sd)[dim2] )
points(clus1_obs[,dim1], clus1_obs[,dim2], col = 2)
points(clus2_obs[,dim1], clus2_obs[,dim2], col = 3)

##### Hyperplane test ######

data <- read.csv("hyperplane.csv")
# Optimal parameters: xi = 15, tau = 0.01
colors = data[,3]
data <- data[,1:2]
plot(data[,1:2], color = data$color)
points(data[color==0,], col = 2)
points(data[color==1,], col = 3)
real_obs = which(color==0)

xi<-15
tau<- .01
model <- CLIQUE(data = data, xi = xi, tau = tau)
plot(data)
biggest_clust(model)

clus_1 <- sort(model[[3]]$objects)
correct = 0

for(i in 1:length(clus_1)){
  for(j in 1:length(real_obs)){
    if(clus_1[i] == real_obs[j]){
      correct = correct + 1
      break
    }
  }
}


plot(data)

points(data[clus_1,], pch = "x")
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)
points(data[clus_5,], col = 6)

ggplot(data =  data) +
  geom_point(aes(x,y))+
  geom_point(data= data[colors==0,], aes(x,y), colour = "#ff71ce", size = 3)+
  geom_point(data= data[colors==1,], aes(x,y), colour = "#01cdfe", size = 3)+
  geom_point(data= data[clus_1,], aes(x,y, shape = "Cluster 1"), color = "black")+
  scale_shape_manual(name = "Clusters", values=c("Cluster 1" = 2))
ggsave("hyperplane_clique.png",device = "png", width = 3000, height = 2000, units = "px")

# No matter what we do, we will never divide the points in two clusters.
# However, we can clearly isolate the denser cluster on the left

##### Blob #####
data <- read.csv("blob.csv")
# Optimal parameters: xi = 25, tau = 0.001
colors = data[,3]
data <- data[,1:2]
plot(data)

# All one = xi 10 tau .001
# Four clusters = xi 10 tau .017
xi<- 30
#tau = .006
tau<- .005
model <- CLIQUE(data = data, xi = xi, tau = tau)
biggest_clust(model)

clus_1 <- model[[2]]$objects
clus_2 <- model[[5]]$objects
clus_3 <- model[[6]]$objects
clus_4 <- model[[8]]$objects
plot(data)
points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)

ggplot(data =  data) +
  geom_point(aes(x,y))+
  geom_point(data= data[colors==0,], aes(x,y), colour = "#ff71ce", size = 3)+
  geom_point(data= data[colors==1,], aes(x,y), colour = "#01cdfe", size = 3)+
  geom_point(data= data[colors==2,], aes(x,y), colour = "#05ffa1", size = 3)+
  geom_point(data= data[colors==3,], aes(x,y), colour = "#b967ff", size = 3)+
  geom_point(data= data[clus_1,], aes(x,y, shape = "Cluster 1"), color = "black")+
  scale_shape_manual(name = "Clusters", values=c("Cluster 1" = 2))
ggsave("blob_clique.png",device = "png", width = 3000, height = 2000, units = "px")

size = 2
ggplot(data =  data) +
  geom_point(aes(x,y), size = size)+
  geom_point(data= data[clus_1,], aes(x,y, color = "Cluster 1"), size = size)+
  #geom_point(data= data[clus_2,], aes(x,y, color = "Cluster 2"), size = size)+
  #geom_point(data= data[clus_3,], aes(x,y, color = "Cluster 3"), size = size)+
  #geom_point(data= data[clus_4,], aes(x,y, color = "Cluster 4"), size = size)+
  scale_color_manual(name = "Clusters", values=c("Cluster 1" = "#ff71ce"))
   # "Cluster 2" = "#01cdfe", "Cluster 3" = "#05ffa1" , "Cluster 4" =  "#b967ff"))

ggsave("blob_clique_alt.png",device = "png", width = 3000, height = 2000, units = "px")

##### Basic5 #####

data <- read.csv("basic5.csv")
# Optimal parameters: xi = 16, tau = 0.002
colors = data[,3]
data <- data[,1:2]
plot(data)

xi<- 16
tau<- .002
model <- CLIQUE(data = data, xi = xi, tau = tau)
biggest_clust(model)

clus_1 <- model[[3]]$objects
clus_2 <- model[[4]]$objects
clus_3 <- model[[5]]$objects
plot(data)
points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)

ggplot(data =  data) +
  geom_point(aes(x,y))+
  geom_point(data= data[colors==0,], aes(x,y), colour = "#ff71ce", size = 3)+
  geom_point(data= data[colors==1,], aes(x,y), colour = "#01cdfe", size = 3)+
  geom_point(data= data[colors==2,], aes(x,y), colour = "#05ffa1", size = 3)+
  geom_point(data= data[clus_1,], aes(x,y, shape = "Cluster 1"), color = "black")+
  geom_point(data= data[clus_2,], aes(x,y, shape = "Cluster 2"))+
  geom_point(data= data[clus_3,], aes(x,y, shape = "Cluster 3"))+
  scale_shape_manual(name = "Clusters", values=c("Cluster 1" = 2, "Cluster 2" = 3, "Cluster 3" = 4))

size = 2
ggplot(data =  data) +
  geom_point(aes(x,y), size = size)+
  geom_point(data= data[clus_1,], aes(x,y, color = "Cluster 1"), size = size)+
  geom_point(data= data[clus_2,], aes(x,y, color = "Cluster 2"), size = size)+
  geom_point(data= data[clus_3,], aes(x,y, color = "Cluster 3"), size = size)+
  #geom_point(data= data[clus_4,], aes(x,y, color = "Cluster 4"), size = size)+
  scale_color_manual(name = "Clusters", values=c("Cluster 1" = "#ff71ce",
    "Cluster 2" = "#01cdfe", "Cluster 3" = "#05ffa1"))


ggsave("basic5_clique.png",device = "png", width = 3000, height = 2000, units = "px")
##### Face #####

data <- read.csv("face.csv")
# Optimal parameters: xi = 15, tau = 0.003
colors = data[,4]
data <- data[,2:3]
plot(data)

xi<- 15
tau<- .003
model <- CLIQUE(data = data, xi = xi, tau = tau)
biggest_clust(model)

clus_1 <- model[[3]]$objects
clus_2 <- model[[4]]$objects
clus_3 <- model[[5]]$objects
clus_4 <- model[[6]]$objects
plot(data)
points(data[clus_1,], col = 2)
points(data[clus_2,], col = 3)
points(data[clus_3,], col = 4)
points(data[clus_4,], col = 5)

ggplot(data =  data) +
  geom_point(aes(x,y))+
  geom_point(data= data[colors==0,], aes(x,y), colour = "#ff71ce", size = 3)+
  geom_point(data= data[colors==1,], aes(x,y), colour = "#01cdfe", size = 3)+
  geom_point(data= data[colors==2,], aes(x,y), colour = "#05ffa1", size = 3)+
  geom_point(data= data[colors==3,], aes(x,y), colour = "#b967ff", size = 3)+
  geom_point(data= data[clus_1,], aes(x,y, shape = "Cluster 1"), color = "black")+
  geom_point(data= data[clus_2,], aes(x,y, shape = "Cluster 2"))+
  geom_point(data= data[clus_3,], aes(x,y, shape = "Cluster 3"))+
  geom_point(data= data[clus_4,], aes(x,y, shape = "Cluster 4"))+
  scale_shape_manual(name = "Clusters", values=c("Cluster 1" = 2, 
                                                 "Cluster 2" = 3, "Cluster 3" = 4, "Cluster 4" = 5))
ggsave("face_clique.png",device = "png", width = 3000, height = 2000, units = "px")
  
size = 2
ggplot(data =  data) +
  geom_point(aes(x,y), size = size)+
  geom_point(data= data[clus_1,], aes(x,y, color = "Cluster 1"), size = size)+
  geom_point(data= data[clus_2,], aes(x,y, color = "Cluster 2"), size = size)+
  geom_point(data= data[clus_3,], aes(x,y, color = "Cluster 3"), size = size)+
  geom_point(data= data[clus_4,], aes(x,y, color = "Cluster 4"), size = size)+
  scale_color_manual(name = "Clusters", values=c("Cluster 1" = "#ff71ce",
                                                 "Cluster 2" = "#01cdfe", "Cluster 3" = "#05ffa1" , "Cluster 4" =  "#b967ff"))


##### True cluster plot #####

blob <- read.csv("blob.csv")
# Optimal parameters: xi = 25, tau = 0.001
colors_blob = blob[,3]
blob <- blob[,1:2]

face <- read.csv("face.csv")
# Optimal parameters: xi = 15, tau = 0.003
colors_face = face[,4]
face <- face[,2:3]

basic <- read.csv("basic5.csv")
# Optimal parameters: xi = 15, tau = 0.003
colors_basic = basic[,3]
basic <- basic[,1:2]

size = 1
face_plot = ggplot(data =  face) +
  geom_point(data= face[colors_face==0,], aes(x,y), colour = "#ff71ce", size = size)+
  geom_point(data= face[colors_face==1,], aes(x,y), colour = "#01cdfe", size = size)+
  geom_point(data= face[colors_face==2,], aes(x,y), colour = "#05ffa1", size = size)+
  geom_point(data= face[colors_face==3,], aes(x,y), colour = "#b967ff", size = size)+
  labs(title = "Face")
face_plot

blob_plot = ggplot(data =  blob) +
  geom_point(data= blob[colors_blob==0,], aes(x,y), colour = "#ff71ce", size = size)+
  geom_point(data= blob[colors_blob==1,], aes(x,y), colour = "#01cdfe", size = size)+
  geom_point(data= blob[colors_blob==2,], aes(x,y), colour = "#05ffa1", size = size)+
  geom_point(data= blob[colors_blob==3,], aes(x,y), colour = "#b967ff", size = size)+
  labs(title = "Blob")
blob_plot

basic_plot = ggplot(data =  basic) +
  geom_point(data= basic[colors_basic==0,], aes(x,y), colour = "#ff71ce", size = size)+
  geom_point(data= basic[colors_basic==1,], aes(x,y), colour = "#01cdfe", size = size)+
  geom_point(data= basic[colors_basic==2,], aes(x,y), colour = "#05ffa1", size = size)+
  geom_point(data= basic[colors_basic==3,], aes(x,y), colour = "#b967ff", size = size)+
  labs(title = "Basic5")
basic_plot

library(gridExtra)
multiplot = grid.arrange(basic_plot,blob_plot, face_plot, ncol=3, nrow = 2)

ggsave("artificial_datasets.png",plot = multiplot,  device = "png",  
       width = 6000, height = 4000, units = "px", dpi = 400)
