# Packages required
  install.packages("tidyverse")   # Data manipulation
  
  install.packages("cluster")     # Clustering algorithms
  
  install.packages("factoextra")  # Clustering visualisation
  
  install.packages("dendxtend")   # For comparing two dendograms
  
  install.packages("ape")

  install.packages("NbClust")
  
  install.packages("clValid")
  
  install.packages("summarytools")

  
# Load libraries
  library(tidyverse)
  library(cluster)
  library(factoextra)
  library(dendextend)
  library(summarytools)
  library(foreign)
  library(ape)
  library(NbClust)
  library(clValid)
  
  
  
#--------------------------------------------------------------------------------------------------------------------
# Set working directory
  setwd("R:/PRISM_DataSci/Stroke Medicine/ESUS Data")

  esus <- read.dta("esus_imputed.dta", convert.factors = FALSE)
  
  # Exclude date variables and other patient IDs from the dataset for clustering
    esus_clust <- esus[c(9:99)]
  
  
#--------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVE SUMMARY
  #view(dfSummary(esus_clean))

  view(dfSummary(esus_clust))


#--------------------------------------------------------------------------------------------------------------------
# CLUSTERING ANALYSIS
  
  # Prepare data for clustering analysis
  
    # 1. Scaling/standardizing the data
         # Convert all factor variables to numeric
           df <- scale(esus_clust)
    
    

#----------------------------------------------------------------------------------------
# Agglomerative Hierarchical clustering
  # Show the first rows
    head(df)
    
# 1. Dissimilarity matrix - commute dissimilarity value
     d <- dist(df, method = "euclidean")
     
# 2. Linkage
     # Hierarchical clustering using Complete Linkage
       hc_complete <- hclust(d, method = "complete")
       
       #  Plot the obtained dendrogram
          plot(hc_complete, cex = 0.6, hang = -1)
     
          
     # Compute with agnes
       hc_agnes <- agnes(df, method = "complete")
       
       #  Agglomerative co-efficient
          # Measures the amount of clustering structure found (values closer to 1
          # suggest strong clustering structure)
          hc_agnes$ac
          #0.7970116
          
          # Methods to assess
            m <- c( "average", "single", "complete", "ward")
            names(m) <- c( "average", "single", "complete", "ward")
          
          # Function to compute coefficient
            ac <- function(x) {
              agnes(df, method = x)$ac
            }
          
            map_dbl(m, ac)
            # average    single      complete      ward 
            # 0.7514840  0.7043804   0.7970116     0.9107829 
        
          
          
# 3. Determine optimal number of clusters
     # Elbow method
       fviz_nbclust(df, FUN = hcut, method = "wss")
            
     # Average Silhouette Method
       fviz_nbclust(df, FUN = hcut, method = "silhouette")     
          
     # Gap Statistic Method
       gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 8, B = 50)
       fviz_gap_stat(gap_stat)
        
        
# 5. Cut the tree into 5 clusters        
     hc <- hclust(d, method = "ward.D2")
     
     # Cut tree into 5 groups
       sub_grp <- cutree(hc, k=5)
       
     # Number of members in each cluster
       table(sub_grp)
       # 127 544  11  38  80 
       
       # Visualise the tree
         plot(hc, cex = 0.6)
         rect.hclust(hc5, k = 5, border = 2:5)
       
       
       # Visualise final clusters
         fviz_cluster(list(cluster = sub_grp, data = df, repel = TRUE, geom = c("point"))) + 
         scale_color_brewer(palette = "Spectral")+
         scale_fill_brewer(palette = "Spectral") +
         theme_minimal()
        
       
       
  
#-------------------------------------------------------------------------------------------------
# Hierarchical K-means Clustering
  # For improving k-means results
        
          
# 1.  Estimating the optimal number of clusters
      # Simplified format:
        fviz_nbclust(df, FUNcluster, method = c("silhouette","wss","gap_stat"))    
      
      # Elbow method
        fviz_nbclust(df, kmeans, method = "wss")+ geom_vline(xintercept = 4, linetype = 2)+
            labs(subtitle = "Elbow method")
          
      # Silhouette method
        fviz_nbclust(df, kmeans, method = "silhouette")+ labs(subtitle = "Silhouette method")
          
        
      # Gap statistic
        # nboot = 50 to keep the function speedy. 
        # recommended value: nboot= 500 for your analysis.
        # Use verbose = FALSE to hide computing progression.
          set.seed(123)
          fviz_nbclust(df, kmeans, nstart = 25, method = "gap_stat", nboot = 50)+
            labs(subtitle = "Gap statistic method")
  
          

      # Combined majority clustering rule
      # Compute the number of clusters
          
          
      
        library(NbClust)
        nb <- NbClust(data = esus_clust, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "complete")
        
      
        # Visualize the result
        library(factoextra)#
        cluster_graph <- fviz_nbclust(nb) + theme_minimal()
        cluster_graph$data[3,2]=5
        cluster_graph$data[5,2]=7
        p <- cluster_graph$data
          
          
          
          barplot(p$freq, names = p$Number_clusters, xlab = "Number of clusters k", 
                  ylab = "Frequency among all 30 indices", main = "Optimal clusters k = 4", col = "green")
          
          p$Number_clusters <- as.numeric(as.character(p$Number_clusters))
          
          
          
          ggplot(p, aes(x = Number_clusters, y = freq)) + geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue4", width = 0.75) +
            xlab("Number of clusters k") + ylab("Frequency among all 30 indices") +
            ggtitle("Optimal clusters k = 4") + theme_bw() + scale_y_continuous(breaks = seq(0, 7, by = 1)) + scale_x_continuous(breaks = seq(0, 15, by = 1))
           
          
       
          
# 2.  Computing clusters - 4 clusters to match published work
      res.hk <- hkmeans(df, 4)
      
      res.pca <- prcomp(df, scale = TRUE)
      
      #2 PCA - 16% variation
      
      scree1 <- fviz_eig(res.pca, ncp = 25, main = "", xlab = "Number of Principal Components")
      scree1
      sum(scree1$data$eig)
      
      #20 PCA - 65.1% variation
      
      scree2 <- fviz_eig(res.pca, ncp = 50, main = "", xlab = "Number of Principal Components")
      scree2
      sum(scree2$data$eig)
      
      #50 PCA - 87.7% variation
      
      scree3 <- fviz_eig(res.pca, ncp = 75, main = "", xlab = "Number of Principal Components")
      scree3
      sum(scree3$data$eig)
      
      #80 PCA - 98.3% variation (thus only 16 redundant features)
      
      
      # Print the results
        res.hk
         
      # Visualise the tree
        fviz_dend(res.hk, cex=0.6, rect = TRUE, rect_fill = TRUE) + 
          scale_color_brewer(palette = "Set1")+
          scale_fill_brewer(palette = "Set1") +
          theme_minimal()
      
      # Visualise the hkmeans final clusters
        fviz_cluster(res.hk, repel = TRUE, geom = "point", xlab = "Principal component 1", ylab = "Principal component 2") + 
          scale_color_brewer(palette = "Set1")+
          scale_fill_brewer(palette = "Set1") +
          theme_minimal()

      # Cluster size
        res.hk$size
        # 177 430  44 149  
        
        
# 2b. Computing clusters -  2 clusters
      res.hk_2 <- hkmeans(df, 2)
        
      # Print the results
        res.hk_2
        
      # Visualise the tree
        fviz_dend(res.hk_2, cex=0.6, rect = TRUE, rect_fill = TRUE) + 
          scale_color_brewer(palette = "Set1")+
          scale_fill_brewer(palette = "Set1") +
          theme_minimal()
        
      # Visualise the hkmeans final clusters
        fviz_cluster(res.hk_2, repel = TRUE, geom = "point") + 
          scale_color_brewer(palette = "Set1")+
          scale_fill_brewer(palette = "Set1") +
          theme_minimal()

        
      # Cluster size
        res.hk_2$size
        # 223 577   
        

          
# 3.  To assign the point classification to the original data
      hk.clusted <- cbind(esus, cluster = res.hk$cluster)
      
      hk.clusted_2 <- cbind(esus, cluster = res.hk_2$cluster)
      
      
      
# 4.  Export the clustered dataset
      write.csv(hk.clusted, file = "esus_clustered.csv", sep = "", row.names = TRUE, 
                  col.names = TRUE)
          
      write.csv(hk.clusted_2, file = "esus_clustered_2.csv", sep = "", row.names = TRUE, 
                col.names = TRUE)    
      
          
