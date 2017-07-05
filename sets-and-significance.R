library(multtest,verbose=FALSE,quietly=TRUE)
library(arules,verbose=FALSE,quietly=TRUE)

source("halt.R")
source("significance_all.R")
source("significance_new.R")

options <- commandArgs(trailingOnly = T)
WD <- getwd();
if (!is.null(WD)) setwd(WD);

input  <- options[1]
support <- as.numeric(options[2])
minlen <- as.numeric(options[3])
maxlen <- as.numeric(options[4])
alpha <- as.numeric(options[5])

if (is.na(options[1])) {
    print("No data given\n")
    halt()
}
if (length(options)<5) {
    print("Wrong input, please check ReadMe.txt\n")
    halt()
}

# read the matrix and extract only meaningful sites and species

  data <- read.table(input,sep=",")
  my_matrix <- data.matrix(data)

# take only the sites with >1 bacteria and bacteria that occur in >2 different sites

  while (length(which(apply(my_matrix,2,sum)<3))+length(which(apply(my_matrix,1,sum)<2))>0) {
    my_matrix <- my_matrix[,which(apply(my_matrix,2,sum)>2)]
    my_matrix <- my_matrix[which(apply(my_matrix,1,sum)>1),]
  }
  places <- dimnames(my_matrix)[[1]]
  names <- dimnames(my_matrix)[[2]]

# provide rows and columns names if needed

  if (length(places)==0) {
    places <- vector("character",length(my_matrix[,1]))
    for (i in 1:length(my_matrix[,1])) {
        places[i] <- paste("SITE",as.character(i),sep="_")
    }
    dimnames(my_matrix)[[1]] <- places
  }

  if (length(names)==0) {
    names <- vector("character",length(my_matrix[1,]))
    for (i in 1:length(my_matrix[1,])) {
        names[i] <- paste("SPECIE",as.character(i),sep="_")
    }
    dimnames(my_matrix)[[2]] <- names
  }
  
# threshold for highly abundant species based on support (75%-quantile of frequencies)

  freq <- apply(my_matrix,2,sum)/length(places)
  F <- as.numeric(quantile(freq)[3])
  threshold <- support / (F*length(places))
  if (!is.na(options[7])) threshold <- as.numeric(options[7])

# getting different frequent sets

  support <- support/length(places)
  res <- apriori(my_matrix,parameter=list(support=support,minlen=minlen,maxlen=maxlen,target="closed frequent itemsets"))
  sets <- labels(res)
  print("apriori candidate sets calculated")

# getting rid of repeats and inclusions - shall be rewritten and optimized

# get rid of { } which are written automatically by apriori
  for (i in 1:length(sets)) sets[i] <- substr(sets[i],2,nchar(sets[i])-1)
# transform vector of characters "sets" to a list of vectors
  sets_list <- lapply(sets, function(x) unlist(strsplit(x,",")))
  sets_list_length <- sapply(sets_list, function(x) length(x))
# for every set check all the sets of smaller length
  i = which(sets_list_length==(minlen+1))[1]
  while (!is.na(sets[i])) {
    N = length(sets_list[[i]])
    j=1
    while (j<which(sets_list_length==N)[1]) {
      if (length(which((sets_list[[j]] %in% sets_list[[i]])==FALSE))==0) {
        sets_list <- sets_list[-j]
        sets <- sets[-j]
        sets_list_length <- sets_list_length[-j]
        i = i-1
        j = j-1
      }
      j = j+1
    }
    i = i+1
  }
# save the number of sets before significance check
  length(sets) -> length.old
  print("apriori sets cleaned")

# checking the significance

  sets_list <- lapply(sets, function(x) unlist(strsplit(x,",")))
  sets_list_length <- sapply(sets_list, function(x) length(x))

  print(paste("Number of sets for significant check: ",length(sets_list),sep=""))

  if (max(sets_list_length) >= 10)
    list <- significance(sets_list, my_matrix, threshold, alpha, minlen)
  if (max(sets_list_length) < 10)
    list <- significance.full(sets_list, my_matrix, alpha, minlen)

  results <- list[[1]]
  sets <- list[[2]]
  print("significance check done")
 
  
# multiple testing correction
  raw_pvalues = as.double(results[,6])
  adjusted = mt.rawp2adjp(raw_pvalues, proc=c("BY") )
  results = cbind(results[,1:5], adjusted$adjp[order(adjusted$index),])
  print(dim(results))

  
# getting rid of non-significant communities and writing them down
  got.out.sets <- vector("numeric")
  count <- 1
  i=1
  while (!is.na(sets[i])) {
    if (!is.matrix(results)) break;
    results_rows_for_set_i <- which(as.numeric(results[,4])==i)
    if (length(results_rows_for_set_i)>0) {
      if (length(which(as.numeric(results[results_rows_for_set_i,7])>alpha))>0) {
        results <- results[-results_rows_for_set_i,]
        got.out.sets[count] <- i
        count <- count+1
      }
    }
    i <- i+1

  }

# delete non-significant sets
  if (length(got.out.sets)>0) sets <- sets[-got.out.sets]
    print("non-significant communities deleted")

  if (!is.matrix(results) || (length(sets)==0)) print("No significant results")

  if (is.matrix(results)) {
#  get rid of spare sets if there were separations in significance check - needs to be optimized as well

    if (length(sets)>length.old) {
    
      sets_list <- lapply(sets, function(x) unlist(strsplit(x,",")))
      sets_list_length <- sapply(sets_list, function(x) length(x))
  
      sets <- sets[order(sets_list_length)]
      sets_list_length <- order(sets_list_length)
  
      sets_list <- lapply(sets, function(x) unlist(strsplit(x,",")))
 
      i=2
      while (!is.na(sets[i])){
        N = length(sets_list[[i]])
        j=1
        while (j<i) {
          if (length(which((sets_list[[j]] %in% sets_list[[i]])==FALSE))==0) {
            sets_list <- sets_list[-j]
            sets <- sets[-j]
            sets_list_length <- sets_list_length[-j]
            i = i-1
            j = j-1
          }
          j = j+1
        }
        i = i+1
      }
  
    }
    print("results cleaned")


# writing down the final answer
    output_results_final <- matrix(0,nrow=length(sets),ncol=5)
    dimnames(output_results_final)[[2]] <- list("Community","Size","Number of common sites","Raw p-value","Corrected p-value")

    for (i in 1:length(sets)) {
      line <- match(sets[i],results[,2])
      if (!is.na(line)) {
        output_results_final[i,1] <- sets[i]
        output_results_final[i,2] <- results[line,5]
        output_results_final[i,3] <- results[line,3]
		output_results_final[i,4] <- results[line,6]
		output_results_final[i,5] <- results[line,7]
      }
      else {
        indices <- sapply(sets_list[i],function(x) match(x,names))
        N <- length(indices)
        common_ones <- length(which(apply(my_matrix[,indices],1,sum)==N))
        output_results_final[i,1] <- sets[i]
        output_results_final[i,2] <- N
        output_results_final[i,3] <- common_ones
	    output_results_final[i,4] <- 0
	    output_results_final[i,5] <- 0
      }
    }
    print("final results calculated")


    if (is.matrix(output_results_final)) {
      output_results_final <- output_results_final[which(as.numeric(output_results_final[,3])>0),]
      output_results_final = output_results_final[order(as.numeric(output_results_final[,2]),decreasing=FALSE),]
      write.table(output_results_final,file="cooccurrence-results.dat", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
	  print("final significant results written in cooccurrence-results.dat")
    }
    if (!is.matrix(output_results_final)) {
      if (as.numeric(output_results_final[3])>0) {
        write.table(output_results_final,file="cooccurrence-results.dat", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        print("final significant results written in cooccurrence-results.dat")
      }
	  else print("No significant communities")
    }

  }


