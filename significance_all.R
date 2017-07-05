# function for testing ALL subsets
significance.full <- function(sets_list, my_matrix, alpha, minlen) {
  
  count = 1
  results <- matrix(0,ncol=6)
  i = 1
  names <- dimnames(my_matrix)[[2]]
  
  while (i<length(sets_list)) {

    # for each set, extract their numbers in original matrix
    ind <- sapply(sets_list[[i]],function(x) match(x,names))
    # put the "bad-community" flag down
    flag = 0
    
    # calculate co-occurence numbers: size of community, number of common sites,
    # and number of sites where neither is present
    N <- length(ind)
    common_ones <- length(which(apply(my_matrix[,ind],1,sum)==N))
    common_zeros <- length(which(apply(my_matrix[,ind],1,sum)==0))
    
    #make some space in the matrix
    S <- 2^(N-1) - 1
    results <- rbind(results,matrix(0,nrow=S,ncol=6))
    
    #start testing indepence for each subset vs all others
    for (k in 1:(as.integer(N/2))) {
      sub <- combn(1:N,k)
      R = length(sub[1,])
      #for even N we need to test half of N/2-communities
      if (k==N/2)
        R <- as.integer(length(sub[1,])/2)
      for (j in 1:R) {
        # ones for first subsets, zeros for second
        c <- length(which(apply(as.matrix(my_matrix[,ind[sub[,j]]]),1,sum)==k)) - common_ones
        # other way round
        b <- length(which(apply(as.matrix(my_matrix[,ind[-sub[,j]]]),1,sum)==N-k)) - common_ones
        fish = fisher.test(matrix(c(common_ones,b,c,common_zeros), nrow = 2), alternative="greater")
        results[count,1] <- paste(k,j)
        results[count,2] <- concatenate(sets_list[[i]])
        results[count,3] <- common_ones
        results[count,4] <- i
        results[count,5] <- N
        results[count,6] <- fish$p.value
        count <- count+1
        #if p-value is bad - just separate the set
        if (as.numeric(fish$p.value) > alpha) { 
          flag = 1
          break
        }
      }
      if (flag==1) break
    }
    #if the set was definitely non-significant, cut it in two and run them again
    if (flag==1) {
        if (length(ind[sub[,j]])>minlen) {
            tmp1 <- sets_list[[i]][sub[,j]]
            sets_list <- list(sets_list,tmp1)
        }
        if (length(ind[-sub[,j]])>minlen) {
            tmp2 <- sets_list[[i]][-sub[,j]]
            sets_list <- list(sets_list,tmp2)
        }
    }

    i <- i+1
  }
  
  #get rid of spare zeros
  results <- results[which(as.numeric(results[,3])>0),]
  
  #write all the sets
  sets <- unique(results[,2])
  
  return (list(results,sets))
}

#function to write down new subsets
concatenate <- function(..., sep=",") {
  paste(...,collapse=sep)
}
