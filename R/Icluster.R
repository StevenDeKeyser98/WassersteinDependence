# This function clusters the columns of the dataset data via hierarchical agglomerative clustering
# Data should be a matrix with observations in the rows and variables in the columns
# Possibilities for est_method are
#
#      - list("Phi", "MI", "Gauss", "omegas") for Gaussian copula mutual information 
#      - list("Phi", "Hel", "Gauss", "omegas") for Gaussian copula half Hellinger distance (normalized)
#      - list("Phi", phi(t), "HAC", type, M) for general phi-dependence with specified function phi(t),
#        estimated by fitting a nested (hierarchical) Archimedean copula of given type, based on a MS sample of size M
#      - list("Phi", phi(t), "NPHAC", type) for general phi-dependence with specified function phi(t), estimated
#        via non-parametric beta kernel estimation and local bandwidth selection by using ML fitted nested Archimedean
#        copula as reference copula
#      - list("Phi", phi(t), "NP", estimator, bw_method) for general phi-dependence with specified function phi(t), estimated
#        via non-parametric beta kernel estimation and local bandwidth selection, either by using a non-parametric beta kernel
#        estimator as reference copula if bw_method = 1, or by using a big oh bandwidth rule if bw_method = 2
#      - list("Phi", phi(t), "Ellip", grid) for general phi-dependence with specified function phi(t), estimated 
#        via the improved MECIP procedure on the specified grid in (0,infty), and parameter selection done in
#        the function est_PHI using the Gaussian generator as reference generator
#      - list("OT", coef, "omegas") for optimal transport measures, either coefficient 1 or coefficient 2
#
# Max_dim indicates the maximal allowed dimension of the dependence measure used, 
# otherwise, the average link function is used
# Norm is a possible normalization function of the dependence measure


Icluster = function(data,est_method,max_dim = Inf,norm = NULL,link = "Average"){
  
  start_time = Sys.time() # Time the algorithm
  
  if(is.null(norm)){norm = function(t){t}} # Identity normalization if none is given
  
  # Function sim is the similarity measure used throughout the algorithm
  
  sim = function(sample,dim,method,est_method){
    # Computes the similarity between the columns of sample corresponding to dimensions specified in dim
    # Method will be "No link" in the algorithm
    # Est_method is as specified in the function Icluster, and will be the argument for est_PHI function
    
    if(est_method[[1]] == "Phi"){ # If (normalized) Phi-dependence measures are used
      
      if(is.character(est_method[[2]]) && est_method[[2]] == "MI"){ # If Gaussian copula mutual information
        
        if(length(est_method[[4]]) == 1){
          return(norm(MI_normal(est_R(sample,est_method[[4]]),dim)))
        } else{
          omega = CVomega(sample,est_method[[4]],5)
          return(norm(MI_normal(est_R(sample,omega),dim)))
        }
      } 
      
      else if(is.character(est_method[[2]]) && est_method[[2]] == "Hel"){ # If Gaussian copula Hellinger distance
        
        if(length(est_method[[4]]) == 1){
          return(norm(Hel_normal(est_R(sample,est_method[[4]]),dim)))
        } else{
          omega = CVomega(sample,est_method[[4]],5)
          return(norm(Hel_normal(est_R(sample,omega),dim)))
        }
      } else{ # If other method using est_PHI function
        return(norm(est_PHI(sample,dim,method,est_method[-c(1,2)],est_method[[2]])))
      }
    }
    
    else if(est_method[[1]] == "OT"){ # If optimal transport measures are used
      
      if(length(est_method[[3]]) == 1){
        
        if(est_method[[2]] == 1){
          return(norm(BWD1(est_R(sample,est_method[[3]]),dim)))
        } else{
          return(norm(BWD2(est_R(sample,est_method[[3]]),dim)))
        }
      } else{
        omega = CVomega(sample,est_method[[3]],5)
        
        if(est_method[[2]] == 1){
          return(norm(BWD1(est_R(sample,omega),dim)))
        } else{
          return(norm(BWD2(est_R(sample,omega),dim)))
        }
      }
    }
  }
  
  # Initialization 
  
  d = ncol(data)
  n = nrow(data)
  colnames(data) = integer(d)
  
  for(j in 1:d){ # Name the columns of data X1,X2,...,Xd
    name = paste("X",toString(j),sep = "")
    colnames(data)[j] = name
  }
  
  hierarchy = sets::set(as.set(colnames(data))) # Initial hierarchy is {{X1,X2,...,Xd}}
  all = hash() # Empty hash object to store all similarities that are computed
  cluster = as.set(colnames(data)) # Initial clustering (partition) of all variables {X1,...,Xd}
  STC_computed = sets::set() # Empty set for storing similarities that are already computed
  finish = 0 # Index indicating when algorithm terminates
  index = 1 # Index indicating which partition has been found 
  stepp = 1 # Index keeping track of the amount of bivariate similarities that are computed
  
  # Main algorithm 
  
  while(finish != 1){
    
    names_cluster = c() # Vector to store all names of the current partition
    
    # If for example cluster = sets::set(c("X34"),c("X2"),c("X3","X4"),c("X25","X1","X22")) = {X34,X2,c(X3,X4),c(X25,X1,X22)}
    # then names_cluster = c("X34","X2","(X3,X4)","(X1 X22 X25)"), not necessarily in this order,
    # but within a cluster, variables are sorted according to index
    
    for(i in cluster){names_cluster = c(names_cluster,set_to_string(sets::set(i)))}
    
    # Next, we order the names according to the first appearing index in each cluster, i.e. for the example
    # names_cluster = c("(X1 X22 X25)","X2","(X3,X4)","X34")
    
    names_cluster = unlist(strsplit(set_to_string(as.set(names_cluster)), "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE))
    
    # Names for all pairwise combinations of which similarities need to be computed are contained in 
    # names_combn, e.g. for the example, names_combn = c("(X1 X22 X25) X2","(X1 X22 X25) (X3 X4)","(X1 X22 X25) X34","X2 (X3 X4)","X2 X34")
    
    names_combn = apply(combn(names_cluster,2), 2, paste, collapse=" ")
    
    # All possible combinations are also contained in the set combn, e.g. for the example
    # combn = {{X2,X34},{X2,c(X3,X4)},{X2,c(X25,X1,X22)},{X34,c(X3,X4)},{X34,c(X25,X1,X22)},{c(X3,X4),c(X25,X1,X22)}}
    # Recall that set elements have no order
    
    combn = set_combn(cluster,2)
    
    # Combinations that have not yet been considered are stored in STC_new, with a total of length(STC_new)
    
    STC_new = set_complement(STC_computed,combn)
    length_STC_new = length(STC_new)
    
    count = 1 # Index keeping track of how many similarities of length(STC_new) are already computed
    
    for(set in STC_new){ # Take a couple for which the similarity is to be computed, e.g. {x2,c(X3,X4)}
      str_set = set_to_string(set) # Unique string representation, e.g. "X2 (X3 X4)"
      dim = c() # Vector for dimensions of two clusters contained in set
      sample = c() # Vector (matrix) for sample data of two clusters contained in set
      
      for(i in set){ # For each of the two clusters, e.g. X2 and c(X3,X4)
        dim = c(dim,length(i)) # Concatenate dimensions, e.g. c(1,2)
        sample = cbind(sample,data[,i,drop = F]) # Concatenate sample data (recall that colnames(data) = X1,...,Xd)
      }
      
      names_to_keep = colnames(sample) # Keep the names of the current sample in consideration
      colnames(sample) = NULL # Then removes colnames of sample
      
      # In case ncol(sample) == 2, we have two univariate random variables
      if(ncol(sample) == 2){ 
        S = sim(sample,dim,"No link",est_method) # Compute bivariate similarity 
        print(paste("Bivariate similarity ", stepp, " out of ", choose(d,2), " is found"))
        stepp = stepp + 1 # Update stepp index 
      } 
      
      # In case ncol(sample) > max_dim, use link function to compute similarity
      else if(ncol(sample) > max_dim){ 
        S = link(names_to_keep,all,dim, link = link) # Compute link similarity
      } else{ # Or compute multivariate similarity
        S = sim(sample,dim,"No link",est_method) 
        print(paste(ncol(sample), " dimensional similarity is found"))
      }
      
      all[[str_set]] = S # Update hash object, where key = unique string representation of the two clusters
      #                           values = similarity between the two clusters
      
      if(choose(d,2) != length_STC_new){ # In case we are not in the first iteration of the algorithm
        # where only bivariate similarities are computed
        print(paste("Similarity", count, " out of ", length_STC_new, " is found"))
        count = count + 1 # Update count index
      }
      
    }
    
    # Val_all is a named numerical vector with names the unique string representations and 
    # values the corresponding similarity
    val_all = unlist(as.list(all))
    val_all = val_all[names(val_all) %in% names_combn] # Only consider those computed in the current iteration
    # Take the ones with largest similarity, e.g. "X2 (X3 X4)" and convert back to a set, e.g. {X2, c(X3,X4)}
    best_set = string_to_set(names(val_all[val_all == max(val_all)])[1])
    # Put all variables included in best_set in a vector ordered according to index
    best_vct = mixedsort(set_to_vct(best_set))
    
    # A new partition is obtained by removing best_set and adding the merged two clusters
    # e.g. remove x2 and c(X3,X4) and add c(X2,X3,X4)
    cluster = set_union(set_symdiff(cluster,best_set),sets::set(best_vct))
    # Add the new partition to the hierarchy
    hierarchy = set_union(hierarchy,sets::set(cluster))
    print(paste("Partition of ", d - index , " elements is found"))
    STC_computed = set_union(STC_new,STC_computed) # Update already computed similarities 
    index = index + 1 # Update index 
    
    if(length(cluster) == 2){ # Stop merging when only two clusters remain (can be altered)
      finish = 1
    }
  }
  end_time = Sys.time()
  print(difftime(end_time, start_time, units='mins')) # Print total running time
  
  # Finalization 
  
  # Complete the hierarchy by adding the final partition (consisting of all variables in one single cluster)
  hierarchy = set_union(hierarchy,sets::set(sets::set(colnames(data))))
  
  # We summarize all partitions in a hash object with keys (aleph1 = first partition, .... alephd = d'th partition)
  # and values the unique string representation of the respective partition
  new_hierarchy = hash()
  pos = 0
  
  for(set in hierarchy){
    set = sort_set(set)
    name = paste("Aleph_",length(hierarchy) - pos, sep = "")
    new_hierarchy[[name]] = set_to_string(set)
    pos = pos + 1
  }
  
  # We compute the average diameters of all partitions
  avg_diameters = integer(length(hierarchy))
  pos = 1
  
  for(set in hierarchy){
    avg_diameters[pos] = avg_diam(set,all)
    pos = pos + 1
  }
  
  # We compute the maximum splits of all partitions
  max_splits = integer(length(hierarchy)-1)
  pos = 1
  
  for(set in hierarchy){
    
    if(length(set) != 1){
      max_splits[pos] = max_split(set,all)
      pos = pos + 1
    }
  }
  end_time = Sys.time()
  print(difftime(end_time, start_time, units='mins')) # Print total running time
  # Return a list containing the "hierarchy" (hash object), "all" similarities computed in the algorithm (hash)
  # the average diameters "diam" (vector) and the maximum splits "split" (vector)
  return(list("hierarchy" = new_hierarchy, "all" = all, "diam" = avg_diameters, "split" = max_splits))
}

sort_set = function(set){
  
  # This function sorts all vectors in set
  # Each vector of the form c("Xi1",...,"Xik") gets sorted according to i1,...,ik
  # and converted to a string of the form "(Xj1 ... Xjl)"
  
  for(i in set){
    
    if(length(i) > 1){
      set = set_union(set_symdiff(set,sets::set(i)),sets::set(paste("(",paste(i[order(as.numeric(gsub('X', '', i)))],collapse=" "),")",sep = "")))
    }
  }
  return(set)
}

set_to_vct = function(set){
  
  # Puts set elements in a vector
  
  vct = c()
  
  for(i in set){
    vct = c(vct,i)
  }
  return(vct)
}

set_to_string = function(set){
  
  # Gives unique string representation of a partition by first sorting within each group,
  # and converting to string vector notation (using sort_set) and then sorting the groups
  # according to their first X index
  # Returns one string representing the clustering 
  # Example: sets::set(c("X34"),c("X2"),c("X3","X4"),c("X25","X1","X22"))
  # gives "(X1 X22 X25) X2 (X3 X4) X34"
  
  vect = set_to_vct(sort_set(set))
  
  for(i in 1:length(vect)){
    
    if(substr(vect[i],1,1) == "("){
      vect[i] = substring(vect[i], 2)
    }
  }
  firsts = c()
  
  for(i in 1:length(vect)){
    firsts = c(firsts, sub(" .*", "", vect[i]))
  }
  order = order(as.numeric(gsub('X', '', firsts)))
  vect = vect[order]
  
  for(i in 1:length(vect)){
    
    if(gsub(" ", "",vect[i], fixed = TRUE) != vect[i]){
      vect[i] = paste("(",vect[i], sep = "")
    }
  }
  return(paste(vect, collapse = ' '))
}

string_to_set = function(string){
  
  # Reverse operation of set_to_string
  
  string = unlist(strsplit(string, "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE))
  set = as.set(NULL)
  
  for(i in 1:length(string)){
    set = set_union(set,sets::set(unlist(strsplit(gsub("\\(|\\)", "", string[i]), "\\s+"))))
  }
  return(set)
}

avg_diam = function(partition,all){
  
  # Average diameter of a partition
  
  diameters = integer(length(partition))
  index = 1
  
  for(i in partition){
    
    if(length(i) == 1){
      min_s = 1
    } else{
      STC = set_combn(i,2)
      S_values = integer(length(STC))
      pos = 1
      
      for(set in STC){
        str_set = set_to_string(set)
        S_values[pos] = all[[str_set]]
        pos = pos + 1
      }
      min_s = min(S_values)
    }
    diameters[index] = min_s
    index = index + 1
  }
  return(mean(diameters))
}

max_split = function(partition,all){
  
  # Maximum split of a partition
  
  splits = integer(length(partition))
  index = 1
  
  for(i in partition){
    rest = set_symdiff(partition,sets::set(i))
    S_values = c()
    
    for(j in i){
      
      for(k in rest){
        
        for(l in k){
          str_set = set_to_string(sets::set(j,l))
          S_values = c(S_values,all[[str_set]])
        }
      }
    }
    splits[index] = max(S_values)
    index = index + 1
  }
  return(max(splits))
}

## Link function used in Icluster algorithm 

link = function(names,all,dim,link){
  
  # This function computes the similarity between two groups by taking the average, single 
  # or complete linkage of all pairwise inter similarities already computed before and contained in all
  # Names contains the variable names of the two groups 
  
  d1 = ifelse(dim[1] == 1,1,dim[1])
  d2 = ifelse(dim[2] == 1,1,dim[2])
  dep = matrix(0,d1,d2)
  
  for(i in 1:d1){
    
    for(j in 1:d2){
      
      if(d1 == 1){
        name1 = names[1]
      } else{
        name1 = names[i]
      }
      
      if(d2 == 1){
        name2 = names[length(names)]
      } else{
        name2 = names[dim[1] + j]
      }
      dep[i,j] = all[[set_to_string(sets::set(name1,name2))]]
    }
  }
  
  if(link == "Average"){
    return(mean(dep))
  }
  
  if(link == "Single"){
    return(min(dep))
  }
  
  if(link == "Complete"){
    return(max(dep))
  }
}

convert_to_hclust = function(Iclust){
  
  # This function converts the Iclust object (output from Icluster functions) to an hclust object
  
  a = list() 
  hier = Iclust$hierarchy
  l = length(hier)
  merge = matrix(0,nrow = l-1, ncol = 2)
  clusters = list()
  
  for(i in 1:l){
    clusters[[i]] = i
  }
  
  for(k in 1:(l-1)){
    name_old = paste("Aleph_",k,sep = "")
    name_new = paste("Aleph_",k+1,sep = "")
    old = as.set(unlist(strsplit(hier[[name_old]], "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE)))
    new = as.set(unlist(strsplit(hier[[name_new]], "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE)))
    change = set_to_vct(set_complement(old,new))
    match = as.numeric(unlist(regmatches(change, gregexpr("[[:digit:]]+", change))))
    
    for(i in 1:(length(clusters)-1)){
      
      for(j in (i+1):length(clusters)){
        
        if((sum(match %in% c(clusters[[i]],clusters[[j]])) == length(match))){
          new_cluster = sort(c(clusters[[i]],clusters[[j]]))
          
          if(i <= l && j <= l){
            merge[k,] = c(-i,-j) 
          }
          if(i <= l && j > l){
            merge[k,] = c(-i,j-l)
          }
          if(i > l && j <= l){
            merge[k,] = c(i-l,-j)
          }
          if(i > l && j > l){
            merge[k,] = c(i-l,j-l)
          }
        }
      }
    }
    clusters[[l+k]] = new_cluster
  }
  
  a$merge = merge
  a$height = 1:(l-1)
  a$order = 1:l 
  names = integer(l)
  
  for(j in 1:l){
    name = paste("X",toString(j),sep = "")
    names[j] = name
  }
  a$labels = names
  class(a) = "hclust"
  return(a)
  
}

plot_indices = function(Iclust,k){
  
  # This function plots a dendrogram of an Iclust object (output from Icluster function)
  # It also plots the average diameter and maximum split as a function of the number of clusters
  # k is the number of colors used in the dendrogram 
  
  dev.new()
  l = length(Iclust$split) + 1
  diams = Iclust$diam[2:l]
  splits = Iclust$split
  d = as.dendrogram(convert_to_hclust(Iclust))
  d %>% dendextend::set("leaves_pch", 19) %>%  
    dendextend::set("leaves_cex", 0.4) %>%  
    dendextend::set("leaves_col", "black") %>%
    dendextend::set("branches_k_color","black") %>% # k = k for colors
    dendextend::set("labels_cex", 0.7) %>%
    plot(main = "Dendrogram", cex.main = 2.5, axes = F) 
  abline(h = 56, lty = 2)
  # text(-6, 90, "6", xpd = NA)
  # dev.new()
  frame = as.data.frame(cbind(seq(1,l),Iclust$diam[1:l]))
  g = ggplot(frame) +
    geom_line(aes(x = V1, y = V2), size = 1, col = "black") +
    xlab("Number of classes") + ylab("Average diameter") + 
    scale_x_continuous(breaks= pretty_breaks()) + 
    theme_bw() + theme(axis.text = element_text(size=15),
                       axis.title = element_text(size=20)) 
  # print(g)
  frame = as.data.frame(cbind(seq(2,l),splits))
  # dev.new()
  g = ggplot(frame) +
    geom_line(aes(x = V1, y = splits), size = 1, col = "black") +
    xlab("Number of classes") + ylab("Maximum split") + 
    scale_x_continuous(breaks= pretty_breaks()) + 
    theme_bw() + theme(axis.text = element_text(size=15),
                       axis.title = element_text(size=20)) 
  # print(g)
  # dev.new()
  frame = as.data.frame(cbind(seq(2,l),1-diams,splits))
  names(frame) = c("NoC","diams","splits")
  # ggplot(frame) + geom_point(aes(x = diams, y = splits), col = "black", size = 3) + 
  #  geom_text(aes(x = diams, y = splits, label = ifelse(NoC <= 10,NoC," ")), 
  #            hjust = 0.5, vjust = -0.8, size = 5) +
  #  geom_segment(aes(xend = diams, yend = splits), x = 0, y = 0, col = 'black', lty = "dotted") + 
  #  xlab("1-Average diameter") + ylab("Maximum split") + xlim(0,1) + ylim(0,1) + 
  #  theme_bw() + theme(axis.text = element_text(size=15),
  #                     axis.title = element_text(size=20)) 
}

Silhouette = function(Iclust){
  
  # Silhouette index (based on dissimilarity = 1 - similarity)
  
  q = length(Iclust$hierarchy)
  all = Iclust$all
  meanss = integer(q)
  
  for(p in 1:q){
    print(paste("p equals", p, "out of", q))
    part = Iclust$hierarchy[[paste("Aleph_",toString(p),sep = "")]]
    part = string_to_set(part)
    ais = integer(length(q))
    bis = ais
    indexais = 1
    indexbis = 1
    
    for(i in part){
      l = length(i)
      
      if(l == 1){
        ais[indexais] = 0 # zero dissimilarity
        indexais = indexais + 1
      } else{
        
        for(j in 1:l){
          Dis_values = integer(l-1)
          pos = 1
          
          for(k in setdiff(1:l,j)){
            Dis_values[pos] = 1 - all[[set_to_string(as.set(c(i[j],i[k])))]]
            pos = pos + 1
          }
          ais[indexais] = mean(Dis_values)
          indexais = indexais + 1
        }
      }
      
      if(l == q){
        bis = rep(Inf,q)
      } else{
        
        for(j in 1:l){
          Dis = integer(length(part)-1)
          rest = set_symdiff(part,sets::set(i))
          pos1 = 1
          
          for(k in rest){
            Dis_values = integer(length(k))
            pos2 = 1
            
            for(m in k){
              Dis_values[pos2] = 1 - all[[set_to_string(as.set(c(i[j],m)))]]
              pos2 = pos2 + 1
            }
            Dis[pos1] = mean(Dis_values)
            pos1 = pos1 + 1
          }
          bis[indexbis] = min(Dis)
          indexbis = indexbis + 1
        }
      }
    }
    ss = integer(q)
    
    for(i in 1:q){
      
      if(ais[i] == 0 || bis[i] == Inf){
        ss[i] = 0
      } else{
        ss[i] = (bis[i] - ais[i])/max(ais[i],bis[i])
      }
    }
    meanss[p] = mean(ss)
  }
 
   return(meanss)
  
}

Redundancy = function(data,Iclust,coef = "MI",est_method = "Gauss", omegas = NULL, phi = NULL,norm = NULL){
  
  # Multivariate redundancy 
  
  if(is.null(norm)){norm = function(t){t}} # Identity normalization if none is given
  d = ncol(data)
  colnames(data) = integer(d)
  
  for(j in 1:d){ # Name the columns of data X1,X2,...,Xd
    name = paste("X",toString(j),sep = "")
    colnames(data)[j] = name
  }
  
  if(is.null(omegas)){
    omega = 1
  } else{
    omega = CVomega(data,omegas,5)
  }
  red = integer(d-1) # Redundancy
  
  for(i in 1:(d-1)){
    print(paste("i = ", i, "out of", d-1))
    current_partition = Iclust$hierarchy[[paste("Aleph_",toString(i),sep = "")]]
    nclusters = d-i+1
    dim_clusters = integer(nclusters)
    sample = c()
    index = 1
    
    for(j in string_to_set(current_partition)){
      dim_clusters[index] = length(j)
      sample = cbind(sample,data[,j])
      index = index + 1
    }
    
    if(est_method[[1]] == "Gauss"){
      R_est = est_R(sample,omega)
      
      if(coef == "MI"){
        red[i] = norm(MI_normal(R_est,dim_clusters))
      } 
      if(coef == "Hel"){
        red[i] = norm(Hel_normal(R_est,dim_clusters))
      }
      if(coef == "BWD1"){
        red[i] = norm(BWD1(R_est,dim_clusters))
      }
      if(coef == "BWD2"){
        red[i] = norm(BWD2(R_est,dim_clusters))
      }
    } else{
      red[i] = norm(est_PHI(sample,dim_clusters,"No link",est_method,phi))
    }
  }
  dev.new()
  frame = as.data.frame(cbind(seq(2,d),red))
  g = ggplot(frame) +
    geom_line(aes(x = V1, y = rev(red)), size = 1, col = "black") +
    xlab("# clusters") + ylab("Redundancy") + 
    scale_x_continuous(breaks= pretty_breaks()) + 
    theme_bw() + theme(axis.text = element_text(size=15),
                       axis.title = element_text(size=20),
                       plot.title = element_text(hjust=0.5, size = 25)) 
  print(g)
  return(frame)
  
}

RI = function(optimal,part,adj = FALSE){
  
  # Computes the (adjusted, if adj = TRUE) rand index between the partition optimal
  # (target partition) and the partition part (computed partition)
  
  picked = integer(length(optimal))
  nclust = length(unique(optimal))
  part = unlist(strsplit(part, "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE))
  
  for(i in 1:nclust){
    ind = as.numeric(unlist(regmatches(part[i], gregexpr("[[:digit:]]+", part[i]))))
    picked[ind] = i
  }
  
  if(adj == TRUE){
    return(adj.rand.index(optimal,picked))
  } else{
    return(rand.index(optimal,picked))
  }
} 
