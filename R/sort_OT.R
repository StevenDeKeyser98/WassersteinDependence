# Sorts columns of sample such that the vector dimensions dim are in ascending order
# Sample is a matrix with the observations in the rows

sort_OT = function(sample,dim){
  
  dim_new = sort(dim)
  sample_new = sample
  order = order(dim)
  k = length(dim)
  start = 1
  
  for(i in 1:k){
    sumdim = sum(dim_new[1:i])
    pick = sum(dim[1:order[i]]) - dim[order[i]] + 1
    sample_new[,start:sumdim] = sample[,pick:(pick + dim[order[i]] - 1)]
    start = sumdim + 1
  }
  
  return(list("sample" = sample_new, "dim" = dim_new))
  
}
