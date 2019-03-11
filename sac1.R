library(igraph)
library(lsa)

# read data
folder <- "/Users/xfang7/Google Drive/Courses/CSC591-603/hw/MarketSegmentation/"
graph = read_graph(paste(folder,"fb_caltech_small_edgelist.txt",sep = ''),format =c("edgelist"))
attribute_data <- read.csv(paste(folder,"fb_caltech_small_attrlist.csv", sep=''),header = T)

# funciton to compute Gain of modularity attribute
get_delta_attr <- function(attrList, h, membership, values)
{
  indices <- which(values == membership)
  similar <- 0
  for(i in indices)
  {
    similar <- similar + cosine(as.numeric(attrList[h,]), as.numeric(attrList[i,]))
  }
  similar <- similar/length(indices)
}

#implement phase 1 of Sac1 algorithm
phase1 <- function(graph,attributes, mapped_communities, alpha){
  #limit maximum of 15 iterations
  for(it in 1:15)
  {
    x <- mapped_communities
    for(i in 1:vcount(graph))
    {
      # index: store the index of  vertex with max deltaQ
      index <- 0
      # maxQ:store the maximum of DeltaW
      maxQ <- 0
      # neighbors of vertex i
      n <- neighbors(graph, i)
      
      for(j in unique(mapped_communities[n]))
      {
        tmp <- mapped_communities
        old_modularity <- modularity(graph,tmp)
        #Remove i from its community, place to j’s community
        tmp[i] <- j
        new_modularity <- modularity(graph,tmp)
        
        #Compute the composite modularity gain
        delta_Q_newman <- new_modularity - old_modularity
        delta_Q_attr <- get_delta_attr(attributes, i, j, mapped_communities)
        delta_Q <- (1-alpha)*delta_Q_attr + (alpha)*delta_Q_newman
        
        # update the vertex index and corresponding max deltaQ
        if(i!=j && delta_Q > maxQ){
          index <- j
          maxQ <- delta_Q
        }
      }
      # move i to j's comunity, where j with maximum positive gain (if exists) 
      if(index !=0){
        mapped_communities[i] <- index
      }
    }
    # if no further improvement in modularity,break 
    if(isTRUE(all.equal(x, mapped_communities)))
    {
      break
    }
    x <- mapped_communities
  }
  mapped_communities  
}

#Phase2 of sac1 algorithm
phase2 <- function(graph,attributes, mapped_communities, alpha){
  x <- mapped_communities
  for(i in 1:15)
  {
    combined_graph <- contract.vertices(graph, mapped_communities)
    new_graph <- simplify(combined_graph, remove.multiple = TRUE, remove.loops = TRUE)
    #reapply phase1
    mapped_communities <- phase1(new_graph, attributes,mapped_communities, alpha)
    # no futher improvement,break
    if(isTRUE(all.equal(x, mapped_communities)))
    {
      break
    }
    x <- mapped_communities
  }
  
  mapped_communities
}

#sac1 algorithm
sac1 <- function(alpha, attributes = attribute_data){
  r1 <- phase1(graph, attributes,alpha=alpha, mapped_communities = c(1:324))
  communities <- phase2(graph, attributes,alpha=alpha, mapped_communities = r1)
  return(communities)
}
# save result to file
save_file <- function(communities, alpha){
  if(alpha == 0.5){
    alpha = 5
  }
  
  file_name<-paste("communities",alpha,sep="_")
  file_name<-paste(file_name,"txt",sep=".")
  f<-file(file_name,"w")
  
  for(i in 1:length(unique(communities)))
  {
    community <- vector("numeric")
    for(j in 1:324)
    {
      if(communities[j]==unique(communities)[i]){
        community <- append(community, j-1, after = length(community))
      }
    }
    cat(as.character(community), file=f, sep = ",")
    cat("\n", file=f)
  }
  
  close(f)
  
}

# read parameter(alpha) from command line
args <- commandArgs(trailingOnly = TRUE)
alpha = as.numeric(args[1])
# run sac1 algorithm and save result to file
result <- sac1(alpha = alpha)
save_file(result, alpha = alpha)

